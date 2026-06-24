# Drift-free integer DDA + sliding-window permute kernel.
#
# Per-lane phase is the *chip index mod L*, advanced by an exact `step_num/step_den`
# chips-per-sample DDA (so it never drifts). We require oversampling (`step ≤ 1` chip
# per sample, i.e. `step_num ≤ step_den`), which guarantees a vector of W consecutive
# samples spans at most W ≤ 64 chips — small enough to read from one 64-chip vpermb
# window with relative indices 0..63. See the module docstring for the full argument.
#
# Index/phase arithmetic uses Int16 (L = 1023 ≪ typemax(Int16)); chips are Int8.

# Chip index mod L (L ≤ 10230 for the longest GNSS primary codes ≪ typemax(Int16)) is
# kept 16-bit (Int32 for the rare long tables) so each phase Vec{64} stays compact.
#
# Fixed-point DDA: the chips-per-sample rate is represented as `step_num / 2^_B` with a
# power-of-two denominator, so the per-lane index/remainder split is a shift (`>> _B`) and
# a mask (`& (step_den-1)`) — no `div`/`SignedMultiplicativeInverse` needed — and the init
# can materialise the running product `p = step_num·sample` once and derive both index and
# remainder from it in a single pass.
#
# `_RemT` is unsigned so the DDA carry `remainder + frac_step` can't go negative. With
# `_B = 30` the remainder lives in `[0, 2^30)` and any single advance adds `frac_step < 2^30`,
# so the (pre-carry) sum stays `< 2^31 < typemax(UInt32)` — no overflow in UInt32. `_B = 30`
# is the largest that keeps that carry in UInt32 (the bound is `B ≤ 31`); it pins the
# rate-quantization "drift" to ≤ N·2^-31 chips over an N-sample run (≈1e-4 chips over a
# 200k-sample epoch — ~250× finer than `_B = 24`). The scalar-tail product `step_num·(N-1)`
# (step_num < 2^30) stays within Int64 for any realistic N (< 2^33).
const _IndexT = Int16
const _RemT   = UInt32

# Fixed-point fractional bits and the (power-of-two) DDA denominator. B = 30 is the largest
# that keeps the DDA carry (rem + frac_step < 2^(B+1)) inside the UInt32 `_RemT`, so the
# steady-state hot loop is unchanged from B = 24 — only the one-time init widens.
const _B       = 30
const _STEP_DEN = Int(1) << _B          # = 2^_B

# Per-lane index vectors `[0, 1, …, W-1]`, built once at load. The AVX-512 init (`_init_rel`)
# materialises the running product `p = step_num·sample` for all 64 lanes at once with a real
# SIMD Int64 multiply (AVX-512 has `vpmullq`), then derives rel/remainder from `p` with a
# vector shift + mask. AVX2/Portable (`_init_state`) has no 64-bit SIMD multiply, so it splits
# `step_num` into two _B/2-bit halves and reconstructs `p` from Int32 partial products (see
# there), using these Int32 lanes.
const _LANES64 = Vec{64,Int64}(ntuple(j -> Int64(j - 1), Val(64)))
@inline @generated _lanes_i32(::Val{W}) where {W} =
    :(Vec{$W,Int32}($(Expr(:tuple, (Int32(j - 1) for j in 1:W)...))))

# Convert a chips-per-sample rate to the fixed-point `(step_num, step_den)` pair used by the
# DDA: `step_den = 2^_B`, `step_num = round(cps · 2^_B)`. Drops `rationalize`'s continued
# fraction (one multiply + round instead). Requires `0 < cps ≤ 1` (oversampling).
@inline function _fixed_point_step(cps::Real)
    sd = _STEP_DEN
    sn = round(Int, cps * sd)
    (sn, sd)
end

"""
    chips_per_sample(code_frequency, sampling_frequency)

Normalised chip rate `code_frequency / sampling_frequency` (chips advanced per sample).
Pass the result to [`generate_code!`](@ref) / [`generate_code`](@ref). Must be ≤ 1
(the front-end must oversample, `sampling_frequency ≥ code_frequency`).
"""
chips_per_sample(code_frequency, sampling_frequency) = code_frequency / sampling_frequency

# ===== generate_code! =====
"""
    generate_code!(out, table, step_numerator, step_denominator; phase=0, backend=…)
    generate_code!(out, table, chips_per_sample::Real;            phase=0, backend=…)
    generate_code!(out, table; code_frequency, sampling_frequency, phase=0, backend=…)

Fill `out::AbstractVector{Int8}` with `table`'s code resampled so the chip index advances
by an exact `step_numerator / step_denominator` chips per sample (drift-free integer DDA).
`phase` is an initial **integer chip offset** (default 0). Requires
`0 < step_numerator ≤ step_denominator` (oversampling) and a sane denominator
`0 < step_denominator ≤ 2^_B` (the fixed-point DDA denominator; see `_fixed_point_step`).
"""
function generate_code!(out::AbstractVector{<:Integer}, table::CodeTable,
                        step_numerator::Integer, step_denominator::Integer;
                        phase::Integer = 0, backend::Backend = default_backend(table))
    (0 < step_denominator ≤ _STEP_DEN) ||
        throw(ArgumentError("need 0 < step_denominator ≤ 2^$_B"))
    (0 < step_numerator ≤ step_denominator) ||
        throw(ArgumentError("need 0 < step_numerator ≤ step_denominator (must oversample, chips/sample ≤ 1)"))
    _generate!(out, table, Int(step_numerator), Int(step_denominator), Int(phase), backend)
    out
end
function generate_code!(out::AbstractVector{<:Integer}, table::CodeTable, cps::Real; kw...)
    sn, sd = _fixed_point_step(cps)
    generate_code!(out, table, sn, sd; kw...)
end
function generate_code!(out::AbstractVector{<:Integer}, table::CodeTable;
                        code_frequency::Real, sampling_frequency::Real, kw...)
    generate_code!(out, table, chips_per_sample(code_frequency, sampling_frequency); kw...)
end

# Initialise one DDA state of W lanes whose first lane is sample `start_sample`:
#   phase[j]     = (div(step_num*(start_sample+j), step_den) + phase_offset) mod L
#   remainder[j] =  mod(step_num*(start_sample+j), step_den)        (j = 0..W-1)
# Phase type `T` is Int16 for tables that fit (fast), Int32 for longer ones (the AVX2
# path can't use the scalar-base rel trick — it regressed — so it widens the phase vector
# instead). AVX-512 uses the rel/scalar-base path (`_init_rel`) and ignores this.
@inline _phase_type(L) = L <= typemax(Int16) ? Int16 : Int32
# Fixed-point init: chip index = `p >> _B`, remainder = `p & (step_den-1)`, p = step_num·sample.
# AVX2 has no 64-bit SIMD multiply and `p` (< 2^37 at _B=30) overflows Int32, so SPLIT
# `step_num = hi·2^H + lo` (H = _B÷2 = 15): the partial products `hi·s` and `lo·s` both fit
# Int32 (native `vpmulld`), and `p>>_B` / `p&(2^_B-1)` are reconstructed from them with shifts
# and masks — no `Vec{W,UInt64}` (which has no AVX2 instruction and wrecked the hot loop's
# register allocation). `@noinline`: runs ONCE per stream at setup, kept out of line so its
# temporaries don't bloat the steady-state loops it feeds (`_generate_simd_avx2!`, iterators).
@noinline function _init_state(::Val{W}, step_num, step_den, L, start_sample, phase_offset, ::Type{T}) where {W,T}
    H      = _B >> 1                                           # split point (15 for _B = 30)
    himask = Int32((1 << H) - 1)
    s      = _lanes_i32(Val(W)) + Int32(start_sample)          # per-lane sample
    A      = Vec{W,Int32}(Int32(step_num >> H)) * s            # hi·s  (fits Int32: hi<2^15, s≤4W-1)
    Bv     = Vec{W,Int32}(Int32(step_num & Int(himask))) * s   # lo·s  (fits Int32)
    lowB   = ((A & himask) << H) + Bv                          # low _B bits of p (+ 1 carry bit)
    remainder = convert(Vec{W,_RemT}, lowB & Int32(step_den - 1))
    chip   = (A >> H) + (lowB >> _B)                           # p >> _B  (chip index, ≤ 4W-1)
    if L < 4W   # pathological short table (chip can exceed 2L): per-lane scalar mod
        phase = convert(Vec{W,T}, Vec{W,Int32}(ntuple(j -> Int32(mod(Int64(@inbounds chip[j]) + phase_offset, L)), Val(W))))
        return (phase, remainder)
    end
    po  = Int32(mod(Int64(phase_offset), Int64(L)))
    idx = chip + po                                            # < 2L
    idx = vifelse(idx >= Int32(L), idx - Int32(L), idx)        # mod L (one conditional subtract)
    (convert(Vec{W,T}, idx), remainder)
end

# One windowed lookup: phase (chip indices mod L for W consecutive samples) -> W chips.
# base = lane-0 chip; relative index of every lane is in 0..W-1 (≤ 63) because we
# oversample, with a +L correction for the (≤1 per vector) wrap across the code boundary.
# AVX-512: one 64-chip window, one vpermb over all 64 lanes.
@inline function _window_lookup(::AVX512, padded, phase::Vec{64,T}, L) where {T}
    base = Int(@inbounds phase[1])
    rel = phase - T(base)
    rel = vifelse(rel < zero(T), rel + T(L), rel)
    index = convert(Vec{64,Int8}, rel)
    window = @inbounds padded[VecRange{64}(base + 1)]   # Vec{64,Int8}, always contiguous
    _permute(window, index)
end

# AVX2: vpshufb shuffles each 128-bit lane independently, so we give each lane its OWN
# 16-chip window — the 32 samples are two independent 16-sample halves (each spans ≤ 15
# chips when oversampled). No 4×pshufb+blend composition needed: a single vpshufb does
# both halves at once. `_AVX2_LOW16` selects per-lane base (lane-0 chip for the low half,
# lane-16 chip for the high half).
const _AVX2_LOW16 = Vec{32,Bool}(ntuple(j -> j <= 16, Val(32)))
@inline function _window_lookup(::AVX2, padded, phase::Vec{32,T}, L) where {T}
    Lc = T(L)
    base_lo = Int(@inbounds phase[1])
    base_hi = Int(@inbounds phase[17])
    basevec = vifelse(_AVX2_LOW16, Vec{32,T}(T(base_lo)), Vec{32,T}(T(base_hi)))
    rel = phase - basevec
    rel = vifelse(rel < zero(T), rel + Lc, rel)
    index = convert(Vec{32,Int8}, rel)
    lo16 = @inbounds padded[VecRange{16}(base_lo + 1)]                     # Vec{16,Int8} → low lane
    hi16 = @inbounds padded[VecRange{16}(base_hi + 1)]                     # Vec{16,Int8} → high lane
    window = shufflevector(lo16, hi16, Val(ntuple(j -> j - 1, Val(32))))   # Vec{32,Int8}
    _pshufb(window, index)
end

# ── relative-index DDA ───────────────────────────────────────────────────────────────
# The value the permute consumes is the window-relative index `rel`, which fits Int8 (it
# is ≤ 63 for AVX-512's single 64-chip window, ≤ 15 per half for AVX2's two 16-chip
# windows). So instead of carrying the Int16 chip index and subtracting a base every
# vector, we carry `rel` directly (Int8 → half the registers) plus a scalar base per
# window, advancing both with the identity  rel[j] += carry[j] − carry[base-of-j],
# base += whole + carry[base-lane]  (mod L). This drops the per-vector subtract, the
# Int16→Int8 narrow and the lane-0 extract — ~1.8× on AVX-512.

# AVX-512: one 64-chip window per stream → one scalar base.
@inline function _init_rel(::Val{64}, step_num, step_den, L, start, phase_offset)
    # Materialise `p = step_num·sample` for all 64 lanes at once (real SIMD multiply), then
    # derive rel (relative chip index, Int8) and remainder from `p` with a vector shift + mask.
    p = Vec{64,Int64}(Int64(step_num)) * (_LANES64 + Int64(start))
    d0 = (Int64(step_num) * Int64(start)) >> _B
    rel = convert(Vec{64,Int8},  (p >> _B) - Int64(d0))
    rem = convert(Vec{64,_RemT},  p & Int64(step_den - 1))
    (rel, rem, Int(mod(d0 + phase_offset, L)))
end
@inline _rel_lookup(::AVX512, padded, rel::Vec{64,Int8}, base::Int) =
    _permute((@inbounds padded[VecRange{64}(base + 1)]), rel)

# AVX2: two independent 16-chip windows (low/high 128-bit lane) → two scalar bases. `rel`
# is relative to lane 0 in the low half and lane 16 in the high half.
function _generate!(out, table::CodeTable, step_num, step_den, phase_offset, ::AVX512)
    _generate_simd_avx512!(out, table, step_num, step_den, phase_offset)
end
function _generate_simd_avx512!(out, table::CodeTable, step_num, step_den, phase_offset)
    W = 64; stride = 4W
    L = table.length; padded = table.padded
    whole_step = div(stride * step_num, step_den)            # chips per stride (base wraps mod L)
    frac_step  = _RemT(mod(stride * step_num, step_den)); modulus = _RemT(step_den)
    zero8 = zero(Vec{64,Int8}); one8 = one(Vec{64,Int8})
    rel1, rem1, b1 = _init_rel(Val(W), step_num, step_den, L, 0,  phase_offset)
    rel2, rem2, b2 = _init_rel(Val(W), step_num, step_den, L, W,  phase_offset)
    rel3, rem3, b3 = _init_rel(Val(W), step_num, step_den, L, 2W, phase_offset)
    rel4, rem4, b4 = _init_rel(Val(W), step_num, step_den, L, 3W, phase_offset)
    num = length(out); chunk_start = 0
    @inbounds while chunk_start + stride <= num
        out[VecRange{W}(chunk_start + 1)]      = _rel_lookup(AVX512(), padded, rel1, b1)
        out[VecRange{W}(chunk_start + W + 1)]  = _rel_lookup(AVX512(), padded, rel2, b2)
        out[VecRange{W}(chunk_start + 2W + 1)] = _rel_lookup(AVX512(), padded, rel3, b3)
        out[VecRange{W}(chunk_start + 3W + 1)] = _rel_lookup(AVX512(), padded, rel4, b4)
        rem1 += frac_step; rem2 += frac_step; rem3 += frac_step; rem4 += frac_step
        c1 = rem1 >= modulus; c2 = rem2 >= modulus; c3 = rem3 >= modulus; c4 = rem4 >= modulus
        rem1 = vifelse(c1, rem1 - modulus, rem1); rem2 = vifelse(c2, rem2 - modulus, rem2)
        rem3 = vifelse(c3, rem3 - modulus, rem3); rem4 = vifelse(c4, rem4 - modulus, rem4)
        h1 = Int(c1[1]); h2 = Int(c2[1]); h3 = Int(c3[1]); h4 = Int(c4[1])  # lane-0 carry
        rel1 += vifelse(c1, one8, zero8) - Int8(h1); rel2 += vifelse(c2, one8, zero8) - Int8(h2)
        rel3 += vifelse(c3, one8, zero8) - Int8(h3); rel4 += vifelse(c4, one8, zero8) - Int8(h4)
        b1 = mod(b1 + whole_step + h1, L); b2 = mod(b2 + whole_step + h2, L)
        b3 = mod(b3 + whole_step + h3, L); b4 = mod(b4 + whole_step + h4, L)
        chunk_start += stride
    end
    # Leftover < 4W samples: do up to 3 full W-blocks as single-stream SIMD before the
    # scalar tail. Streams 1–3 are already positioned at chunk_start + {0,W,2W}, so reuse them.
    @inbounds if chunk_start + W <= num
        out[VecRange{W}(chunk_start + 1)] = _rel_lookup(AVX512(), padded, rel1, b1); chunk_start += W
        if chunk_start + W <= num
            out[VecRange{W}(chunk_start + 1)] = _rel_lookup(AVX512(), padded, rel2, b2); chunk_start += W
            if chunk_start + W <= num
                out[VecRange{W}(chunk_start + 1)] = _rel_lookup(AVX512(), padded, rel3, b3); chunk_start += W
            end
        end
    end
    _generate_tail!(out, table, step_num, step_den, phase_offset, chunk_start + 1)
end

# AVX2 keeps the phase-based DDA + two-window vpshufb lookup (the AVX-512 rel/scalar-base
# trick regressed AVX2: the per-half "halfcarry" vector it needs costs more than the
# in-vector subtract + Int16→Int8 narrow it saves — measured ~95 vs ~56 ps). The phase
# type adapts to the table: Int16 (~70 ps) for length ≤ 32767, Int32 for longer tables such
# as L1C-P's 122760-entry TMBOC table — still ~20× over the scalar fallback.
#
# Four interleaved streams (stride = 4W = 128), as on AVX-512, so the DDA carry chains
# overlap for full throughput. (The UInt32 remainder for the fine `_B = 30` rate is fine
# at 4 streams — the earlier ~2400 ps/sample regression was `_init_state`'s SIMD
# widening-multiply build poisoning the loop's register allocation, not stream pressure;
# see `_init_state`.) Byte-identical to the AVX-512 / Portable output.
function _generate!(out, table::CodeTable, step_num, step_den, phase_offset, ::AVX2)
    _check_avx2_length(table)
    if table.length <= typemax(Int16)
        _generate_simd_avx2!(out, table, step_num, step_den, phase_offset, Int16)
    else
        _generate_simd_avx2!(out, table, step_num, step_den, phase_offset, Int32)
    end
end
@inline _check_avx2_length(table::CodeTable) =
    table.length <= typemax(Int32) || throw(ArgumentError(
        "AVX2 backend supports code length ≤ $(Int(typemax(Int32))); table length is " *
        "$(table.length). Use backend=AVX512() or Portable()."))
function _generate_simd_avx2!(out, table::CodeTable, step_num, step_den, phase_offset, ::Type{T}) where {T}
    W = 32; stride = 4W
    L = table.length; padded = table.padded
    Lc = T(L)
    whole_step = T(div(stride * step_num, step_den) % L)   # integer chips per stride, mod L
    frac_step  = _RemT(mod(stride * step_num, step_den))
    modulus    = _RemT(step_den)
    phase1, rem1 = _init_state(Val(W), step_num, step_den, L, 0,  phase_offset, T)
    phase2, rem2 = _init_state(Val(W), step_num, step_den, L, W,  phase_offset, T)
    phase3, rem3 = _init_state(Val(W), step_num, step_den, L, 2W, phase_offset, T)
    phase4, rem4 = _init_state(Val(W), step_num, step_den, L, 3W, phase_offset, T)
    num = length(out); chunk_start = 0
    @inbounds while chunk_start + stride <= num
        out[VecRange{W}(chunk_start + 1)]      = _window_lookup(AVX2(), padded, phase1, L)
        out[VecRange{W}(chunk_start + W + 1)]  = _window_lookup(AVX2(), padded, phase2, L)
        out[VecRange{W}(chunk_start + 2W + 1)] = _window_lookup(AVX2(), padded, phase3, L)
        out[VecRange{W}(chunk_start + 3W + 1)] = _window_lookup(AVX2(), padded, phase4, L)
        rem1 += frac_step; rem2 += frac_step; rem3 += frac_step; rem4 += frac_step
        c1 = rem1 >= modulus; c2 = rem2 >= modulus; c3 = rem3 >= modulus; c4 = rem4 >= modulus
        rem1 = vifelse(c1, rem1 - modulus, rem1); rem2 = vifelse(c2, rem2 - modulus, rem2)
        rem3 = vifelse(c3, rem3 - modulus, rem3); rem4 = vifelse(c4, rem4 - modulus, rem4)
        phase1 += whole_step; phase2 += whole_step; phase3 += whole_step; phase4 += whole_step
        phase1 = vifelse(c1, phase1 + one(T), phase1); phase2 = vifelse(c2, phase2 + one(T), phase2)
        phase3 = vifelse(c3, phase3 + one(T), phase3); phase4 = vifelse(c4, phase4 + one(T), phase4)
        # phase ∈ [0, 2L-1] after the advance → one conditional subtract restores [0, L-1]
        phase1 = vifelse(phase1 >= Lc, phase1 - Lc, phase1); phase2 = vifelse(phase2 >= Lc, phase2 - Lc, phase2)
        phase3 = vifelse(phase3 >= Lc, phase3 - Lc, phase3); phase4 = vifelse(phase4 >= Lc, phase4 - Lc, phase4)
        chunk_start += stride
    end
    # Leftover < 4W samples: up to 3 full W-blocks as single-stream SIMD before the scalar
    # tail. Streams 1–3 are already positioned at chunk_start + {0,W,2W}, so reuse them.
    @inbounds if chunk_start + W <= num
        out[VecRange{W}(chunk_start + 1)] = _window_lookup(AVX2(), padded, phase1, L); chunk_start += W
        if chunk_start + W <= num
            out[VecRange{W}(chunk_start + 1)] = _window_lookup(AVX2(), padded, phase2, L); chunk_start += W
            if chunk_start + W <= num
                out[VecRange{W}(chunk_start + 1)] = _window_lookup(AVX2(), padded, phase3, L); chunk_start += W
            end
        end
    end
    _generate_tail!(out, table, step_num, step_den, phase_offset, chunk_start + 1)
end

function _generate!(out, table::CodeTable, step_num, step_den, phase_offset, ::Portable)
    _generate_tail!(out, table, step_num, step_den, phase_offset, 1)
end

@inline function _generate_tail!(out, table::CodeTable, step_num, step_den, phase_offset, sample)
    L = table.length; chips = table.chips
    # Fixed-point: with step_den = 2^_B the per-sample index is a shift (no idiv); the
    # step_den arg is unused but kept for the call-site signature.
    @inbounds while sample <= length(out)
        index = mod((step_num * Int64(sample - 1)) >> _B + phase_offset, L)
        out[sample] = chips[index + 1]
        sample += 1
    end
    out
end
