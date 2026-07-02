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

# Wrap a scalar chip index back into [0, L) with a single conditional subtract instead of
# an integer `mod` (idiv). Valid whenever the pre-wrap value is < 2L — which holds for the
# AVX-512 base advance once the per-stride whole-chip step is pre-reduced mod L (so
# `base + whole + carry ∈ [0, 2L-1]`), exactly the assumption the AVX2/NEON windowed phase
# DDA already makes. Removes the only idiv from the AVX-512 steady-state loop.
@inline _wrapL(x::Int, L::Int) = x >= L ? x - L : x

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
                        phase::Integer = 0, rem0::Integer = 0,
                        backend::Backend = default_backend(table))
    (0 < step_denominator ≤ _STEP_DEN) ||
        throw(ArgumentError("need 0 < step_denominator ≤ 2^$_B"))
    (0 < step_numerator ≤ step_denominator) ||
        throw(ArgumentError("need 0 < step_numerator ≤ step_denominator (must oversample, chips/sample ≤ 1)"))
    (0 ≤ rem0 < step_denominator) ||
        throw(ArgumentError("need 0 ≤ rem0 < step_denominator (fractional sub-chip offset)"))
    if out isa AbstractVector{Int8} && _use_boundary(Int(step_numerator), Int(step_denominator), backend)
        _generate_boundary!(out, table, Int(step_numerator), Int(step_denominator), Int(phase), _RemT(rem0))
    else
        _generate!(out, table, Int(step_numerator), Int(step_denominator), Int(phase), backend, _RemT(rem0))
    end
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
# temporaries don't bloat the steady-state loops it feeds (`_generate_simd_windowed!`, iterators).
@noinline function _init_state(::Val{W}, step_num, step_den, L, start_sample, phase_offset, ::Type{T},
                               rem0::_RemT = _RemT(0)) where {W,T}
    H      = _B >> 1                                           # split point (15 for _B = 30)
    himask = Int32((1 << H) - 1)
    s      = _lanes_i32(Val(W)) + Int32(start_sample)          # per-lane sample
    A      = Vec{W,Int32}(Int32(step_num >> H)) * s            # hi·s  (fits Int32: hi<2^15, s≤4W-1)
    Bv     = Vec{W,Int32}(Int32(step_num & Int(himask))) * s   # lo·s  (fits Int32)
    lowB   = ((A & himask) << H) + Bv                          # low _B bits of p (+ 1 carry bit)
    remainder = convert(Vec{W,_RemT}, lowB & Int32(step_den - 1))
    chip   = (A >> H) + (lowB >> _B)                           # p >> _B  (chip index, ≤ 4W-1)
    local phase::Vec{W,T}
    if L < 4W   # pathological short table (chip can exceed 2L): per-lane scalar mod
        phase = convert(Vec{W,T}, Vec{W,Int32}(ntuple(j -> Int32(mod(Int64(@inbounds chip[j]) + phase_offset, L)), Val(W))))
    else
        po  = Int32(mod(Int64(phase_offset), Int64(L)))
        idx = chip + po                                            # < 2L
        idx = vifelse(idx >= Int32(L), idx - Int32(L), idx)        # mod L (one conditional subtract)
        phase = convert(Vec{W,T}, idx)
    end
    # Seed the fractional sub-chip offset (one DDA micro-step, frac = rem0, whole = 0). Since
    # rem0 < step_den and each lane's remainder < step_den, the sum < 2·step_den ⇒ carry ∈ {0,1}.
    # Folding rem0 into the Int32 split-multiply above would overflow `lowB` (near the Int32
    # edge), so we add it here instead — done once at setup, never in the hot loop.
    if rem0 != _RemT(0)
        phase, remainder = _seed_micro_phase(phase, remainder, rem0, _RemT(step_den), T(L))
    end
    (phase, remainder)
end

# One DDA micro-advance (frac = rem0, whole = 0) for the phase-form state. Adds the fractional
# sub-chip start offset to every lane's running remainder, carrying ≤1 chip. Used to seed a
# fractional start phase at setup; the steady-state advance is unchanged.
@inline _seed_micro_phase(phase::Vec{W,T}, rem::Vec{W,_RemT}, rem0::_RemT,
                          modulus::_RemT, Lc::T) where {W,T} =
    _advance_phase(phase, rem, rem0, zero(T), modulus, Lc)   # whole = 0: pure fractional seed

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

# NEON: a 128-bit register is 16 bytes, so `tbl1` does a single 16-chip → 16-lane lookup —
# exactly ONE of AVX2's two `vpshufb` halves. So this is the single-window analogue of the
# AVX2 two-window lookup: one base (lane-0 chip), one 16-chip window, one `tbl1`. The 16
# samples span ≤ 15 chips when oversampled, so the relative indices are 0..15 (never the
# ≥16 out-of-range that `tbl1` zeroes).
@inline function _window_lookup(::Neon, padded, phase::Vec{16,T}, L) where {T}
    base = Int(@inbounds phase[1])
    Lc = T(L)
    rel = phase - T(base)
    rel = vifelse(rel < zero(T), rel + Lc, rel)
    index = convert(Vec{16,Int8}, rel)
    window = @inbounds padded[VecRange{16}(base + 1)]   # Vec{16,Int8}, always contiguous
    _tbl1(window, index)
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
@inline function _init_rel(::Val{64}, step_num, step_den, L, start, phase_offset, rem0::_RemT = _RemT(0))
    # Materialise `p = step_num·sample + rem0` for all 64 lanes at once (real SIMD multiply),
    # then derive rel (relative chip index, Int8) and remainder from `p` with a vector shift +
    # mask. `rem0 ∈ [0, step_den)` is the fixed-point fractional sub-chip start offset, folded
    # directly into the Int64 product (exact, overflow-safe). rem0 = 0 ⇒ original integer-phase.
    p = Vec{64,Int64}(Int64(step_num)) * (_LANES64 + Int64(start)) + Int64(rem0)
    d0 = (Int64(step_num) * Int64(start) + Int64(rem0)) >> _B
    rel = convert(Vec{64,Int8},  (p >> _B) - Int64(d0))
    rem = convert(Vec{64,_RemT},  p & Int64(step_den - 1))
    (rel, rem, Int(mod(d0 + phase_offset, L)))
end
@inline _rel_lookup(::AVX512, padded, rel::Vec{64,Int8}, base::Int) =
    _permute((@inbounds padded[VecRange{64}(base + 1)]), rel)
@inline _rel_lookup(::Neon, padded, rel::Vec{16,Int8}, base::Int) =
    _tbl1((@inbounds padded[VecRange{16}(base + 1)]), rel)

# ── one-shot single-window kernel (AVX-512 W=64 `vpermb`, NEON W=16 `tbl1`) ──────────────
# Split-constant recompute. Within a W-sample block the window-relative index is, exactly,
#   rel[j] = ⌊(r + step_num·j) / 2^_B⌋ = base_lane[j] + [ r + frac_lane[j] ≥ 2^_B ]
# where base_lane[j] = (step_num·j) >> _B and frac_lane[j] = (step_num·j) & (2^_B−1) are
# COMPILE-TIME-CONSTANT per-lane vectors (r < 2^_B and frac_lane[j] < 2^_B ⇒ the sum is a
# single 0/1 carry mask, no 64-bit product). Only a scalar `(r, base)` pair advances between
# blocks, so there is no loop-carried per-lane remainder and no lane-0 vector→scalar extract —
# ~1.4× (AVX-512) / ~1.6× (AVX2) over the previous four-stream phase/rel DDA. `Wstep_whole` is
# reduced mod L so the single-subtract `_wrapL` stays valid (needs L ≥ W; every GNSS code is
# ≫ 64). Byte-identical to the Portable oracle. NEON's `tbl1` llvmcall is aarch64-only, but
# this kernel is only instantiated for `backend = Neon()`, which never happens on x86.
function _generate!(out, table::CodeTable, step_num, step_den, phase_offset, ::AVX512, rem0::_RemT = _RemT(0))
    _generate_simd_single!(out, table, step_num, step_den, phase_offset, AVX512(), Val(64), rem0)
end
function _generate!(out, table::CodeTable, step_num, step_den, phase_offset, ::Neon, rem0::_RemT = _RemT(0))
    _check_windowed_length(table, Neon())
    _generate_simd_single!(out, table, step_num, step_den, phase_offset, Neon(), Val(16), rem0)
end
function _generate_simd_single!(out, table::CodeTable, step_num, step_den, phase_offset,
                                backend, ::Val{W}, rem0::_RemT = _RemT(0)) where {W}
    L = table.length; padded = table.padded
    SN = Int64(step_num); mask = Int64(step_den - 1); sd = _RemT(step_den)
    base_lane = Vec{W,Int8}(ntuple(j -> Int8((SN * (j - 1)) >> _B), Val(W)))
    frac_lane = Vec{W,UInt32}(ntuple(j -> UInt32((SN * (j - 1)) & mask), Val(W)))
    thr   = Vec{W,UInt32}(sd); one8 = one(Vec{W,Int8}); zero8 = zero(Vec{W,Int8})
    Wstep = Int64(W) * SN
    Wstep_whole = Int(Wstep >> _B) % L             # chips per block (≤ W), reduced mod L
    Wstep_frac  = _RemT(Wstep & mask)
    r    = _RemT(Int64(rem0) & mask)               # block-start fractional phase (< 2^_B)
    base = Int(mod(Int64(rem0) >> _B + phase_offset, L))
    num  = length(out); blk = 0
    @inbounds while blk + W <= num
        rel = base_lane + vifelse(frac_lane + Vec{W,UInt32}(r) >= thr, one8, zero8)
        out[VecRange{W}(blk + 1)] = _rel_lookup(backend, padded, rel, base)
        t = r + Wstep_frac
        c = t >= sd ? 1 : 0
        r = _RemT(Int64(t) & mask)
        base = _wrapL(base + Wstep_whole + c, L)
        blk += W
    end
    _generate_tail!(out, table, step_num, step_den, phase_offset, blk + 1, rem0)
end

# ── one-shot AVX2 kernel (W=32, two `vpshufb` 16-chip windows) ────────────────────────────
# `vpshufb` shuffles each 128-bit lane independently, so the 32 lanes are TWO independent
# 16-lane single-window sub-blocks: the low half (lanes 0..15, window `base`) and the high
# half (lanes 16..31, window `base_hi`) — the latter is the same recurrence advanced 16
# samples, derived from the shared scalar `r` with one cheap scalar micro-step per block
# (`r16`, `base_hi`). One `vpshufb` over the packed `[lo;hi]` windows. This split-constant
# form does NOT carry a per-half remainder vector, so it avoids the "halfcarry" cost that made
# the earlier phase/rel DDA regress on AVX2 — it is ~1.6× faster than that kernel here.
@inline _rel_lookup_avx2(padded, rel::Vec{32,Int8}, base_lo::Int, base_hi::Int) = begin
    lo16 = @inbounds padded[VecRange{16}(base_lo + 1)]
    hi16 = @inbounds padded[VecRange{16}(base_hi + 1)]
    window = shufflevector(lo16, hi16, Val(ntuple(j -> j - 1, Val(32))))
    _pshufb(window, rel)
end
function _generate!(out, table::CodeTable, step_num, step_den, phase_offset, ::AVX2, rem0::_RemT = _RemT(0))
    _check_windowed_length(table, AVX2())
    _generate_simd_avx2!(out, table, step_num, step_den, phase_offset, rem0)
end
function _generate_simd_avx2!(out, table::CodeTable, step_num, step_den, phase_offset, rem0::_RemT = _RemT(0))
    W = 32; H = 16
    L = table.length; padded = table.padded
    SN = Int64(step_num); mask = Int64(step_den - 1); sd = _RemT(step_den)
    # 16-lane constants (each half is relative to its own lane 0)
    base16 = Vec{32,Int8}(ntuple(j -> Int8((SN * ((j - 1) % H)) >> _B), Val(32)))
    frac16 = Vec{32,UInt32}(ntuple(j -> UInt32((SN * ((j - 1) % H)) & mask), Val(32)))
    thr = Vec{32,UInt32}(sd); one8 = one(Vec{32,Int8}); zero8 = zero(Vec{32,Int8})
    Hfrac = _RemT((Int64(H) * SN) & mask); Hwhole = Int((Int64(H) * SN) >> _B) % L   # 16-sample step
    Wfrac = _RemT((Int64(W) * SN) & mask); Wwhole = Int((Int64(W) * SN) >> _B) % L   # 32-sample step
    r    = _RemT(Int64(rem0) & mask)
    base = Int(mod(Int64(rem0) >> _B + phase_offset, L))
    num  = length(out); blk = 0
    @inbounds while blk + W <= num
        # high-half window phase (scalar advance of 16 samples from the low half)
        th = r + Hfrac; ch = th >= sd ? 1 : 0
        r16 = _RemT(Int64(th) & mask)
        base_hi = _wrapL(base + Hwhole + ch, L)
        rbroad = vifelse(_AVX2_LOW16, Vec{32,UInt32}(r), Vec{32,UInt32}(r16))
        rel = base16 + vifelse(frac16 + rbroad >= thr, one8, zero8)
        out[VecRange{W}(blk + 1)] = _rel_lookup_avx2(padded, rel, base, base_hi)
        t = r + Wfrac; c = t >= sd ? 1 : 0
        r = _RemT(Int64(t) & mask)
        base = _wrapL(base + Wwhole + c, L)
        blk += W
    end
    _generate_tail!(out, table, step_num, step_den, phase_offset, blk + 1, rem0)
end
@inline _check_windowed_length(table::CodeTable, be) =
    table.length <= typemax(Int32) || throw(ArgumentError(
        "$(backend_name(be)) backend supports code length ≤ $(Int(typemax(Int32))); table " *
        "length is $(table.length). Use backend=Portable() (or AVX512() on x86)."))

function _generate!(out, table::CodeTable, step_num, step_den, phase_offset, ::Portable, rem0::_RemT = _RemT(0))
    _generate_tail!(out, table, step_num, step_den, phase_offset, 1, rem0)
end

@inline function _generate_tail!(out, table::CodeTable, step_num, step_den, phase_offset, sample,
                                 rem0::_RemT = _RemT(0))
    L = table.length; chips = table.chips
    # Fixed-point: with step_den = 2^_B the per-sample index is a shift (no idiv); the
    # step_den arg is unused but kept for the call-site signature. `rem0` is the fractional
    # sub-chip start offset added inside the floor (rem0 = 0 ⇒ original integer-phase tail).
    r0 = Int64(rem0)
    @inbounds while sample <= length(out)
        index = mod((step_num * Int64(sample - 1) + r0) >> _B + phase_offset, L)
        out[sample] = chips[index + 1]
        sample += 1
    end
    out
end

# ── boundary fill (high-oversampling fast path) ────────────────────────────────────────
# When the chip rate is far below the sample rate (`m = 2^_B ÷ step_num` samples per chip
# is large), the windowed permute is wasteful: a W-sample block covers only W/m distinct
# chips, yet the permute recomputes all W lanes. Here we iterate over CHIPS instead: within
# a segment anchored at sample `pos` (chip `c0`, fractional chip phase `r0 < 2^_B`), chip
# `c0 + i` starts at, exactly,
#     pos + ⌈(i·2^_B − r0) / step_num⌉
# and is emitted as ONE splat store of width `SW = nextpow2(m+2) ∈ {16,32,64}` at that
# boundary (a run is ≤ m+1 < SW samples, and stores are issued left-to-right, so the next
# chip's store overwrites the overhang — the same trick as the original `gen_code!`).
# Runs longer than SW (m > SW−2, only possible at SW = 64) get extra SW-strided interior
# stores (the `EXTRAS` variant). Segments re-anchor `(c0, r0, pos)` once per code period,
# so the chip-index wrap costs nothing per chip.
#
# The ceil-division is EXACT: for non-power-of-two `step_num` it uses the (64+ℓ)-bit
# round-up magic reciprocal `magic = ⌈2^(64+ℓ)/step_num⌉` (ℓ = ⌊log₂ step_num⌋), which by
# the classical Granlund–Montgomery bound reproduces ⌈d/step_num⌉ for EVERY dividend
# d < 2^64 — far above the ≤ L·2^_B < 2^50 a segment can produce. Power-of-two rates are a
# plain shift. So the output is byte-identical to the permute kernels / Portable oracle for
# any fill length and any continuation split — unlike the `freqfix` run-fill this replaces,
# whose rounded reciprocal drifted ~1–2 samples per ~10⁵ chips and needed a `_RUNFILL_MAXFILL`
# segmentation cap and a `Val{NI}` store-width ladder. One boundary calculation is one
# independent multiply-high (no loop-carried arithmetic, unrolled 4×), one chip value is one
# byte load + broadcast, so the steady state is ~1 store per chip at store throughput. There
# is no vector arithmetic at all — the splat stores are plain unaligned SIMD.jl stores — so
# this ONE kernel serves every backend (and any table length; no `_check_windowed_length`).

# Boundary-step instruction selection: x86 favours the accumulated add/adc product,
# aarch64 the independent `umulh`s (see the two branches in `_boundary_fill!`). Both are
# bit-identical; compile-time constant, so the dead branch is eliminated.
const _BOUNDARY_ACCUM = Sys.ARCH in (:x86_64, :i686)

# ⌊log₂ SN⌋ and the round-up reciprocal; magic == 0 marks the power-of-two shift path.
@inline _boundary_lg(SN::Int64) = 63 - leading_zeros(SN)
@inline _boundary_magic(SN::Int64, lg::Int) =
    SN == Int64(1) << lg ? UInt64(0) :
        UInt64((((UInt128(1) << (64 + lg)) + UInt128(SN) - 1) ÷ UInt128(SN)) % UInt128)

# Branch-free exact ⌈d/SN⌉ for d ≥ 1 via ⌊(d−1)·magic / 2^(64+lg)⌋ + 1. `% UInt64` avoids
# the checked-conversion branch and `lg & 31` makes the shift range provable (lg ≤ _B = 30),
# so the whole thing is one `mulx` + one shift in the hot loop.
@inline _boundary_ceildiv(d::Int64, lg::Int, magic::UInt64) =
    magic == 0 ? Int((d - 1) >>> (lg & 31)) + 1 :
                 Int((((widemul((d - 1) % UInt64, magic) >> 64) % UInt64) >>> (lg & 31)) % Int64) + 1

# Scalar leftover (< SW samples): plain per-sample DDA walk from `(c0, r0)` at `pos`.
@inline function _boundary_scalar_tail!(out, padded, L::Int, SN::Int64, c0::Int, r0::Int64, pos::Int)
    mask = Int64(_STEP_DEN - 1)
    @inbounds while pos < length(out)
        out[pos + 1] = padded[c0 + 1]
        r = r0 + SN
        c0 += Int(r >> _B)
        r0 = r & mask
        c0 >= L && (c0 -= L)
        pos += 1
    end
    nothing
end

# Core boundary fill from continuation state `(c0, r0)` (chip index, fractional chip phase).
# Requires a pointer-compatible contiguous Int8 `out` (the same contract as the permute
# kernels' `VecRange` stores). Writes exactly `length(out)` samples; the caller derives the
# advanced state arithmetically (see `fill_continue!`).
function _boundary_fill!(out, padded, L::Int, SN::Int64, c0::Int, r0::Int64,
                         ::Val{SW}, ::Val{EXTRAS}) where {SW,EXTRAS}
    lg = _boundary_lg(SN)
    magic = _boundary_magic(SN, lg)
    num = length(out)
    pos = 0
    GC.@preserve out padded begin
        po = pointer(out); pp = pointer(padded)
        @inbounds while pos + SW <= num
            nc = L - c0                                     # chips before the code wraps
            # last chip i whose boundary jn keeps the SW-wide store in bounds (jn + SW ≤ num)
            ibound = (Int64(num - SW - pos) * SN + r0) >> _B
            nlim = ibound < nc ? Int(ibound) : nc
            vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + 1)), po + pos, nothing)
            jn = pos
            i = 1
            if !EXTRAS
                # a run is ≤ m+1 < SW: exactly one store per chip, at its own boundary
                if magic == 0                               # power-of-two rate: shift only
                    d = (Int64(1) << _B) - r0 - 1
                    pos1 = pos + 1
                    while i <= nlim
                        jn = pos1 + Int(d >>> (lg & 31))
                        vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + i + 1)), po + jn, nothing)
                        d += Int64(1) << _B
                        i += 1
                    end
                elseif _BOUNDARY_ACCUM
                    # x86: ACCUMULATED 128-bit product. P_i = (d_i−1)·magic advances by the
                    # constant 2^_B·magic, kept as a manual (hi, lo) pair so it lowers to
                    # add/adc (LLVM's Int128 lowering is worse) — the boundary is hi >> lg,
                    # so each chip costs one add/adc + shift instead of a multiply-high.
                    # Bit-identical to `_boundary_ceildiv`; measured ~1.15–1.2× over the
                    # mulx form on Zen 5 and it recovered the CI AVX2 GPSL1CA @ 5 MHz row
                    # (0.68 → 0.95 vs run-fill).
                    Clo = (magic << _B) % UInt64
                    Chi = (magic >>> (64 - _B)) % UInt64
                    P = widemul((((Int64(1) << _B) - r0 - 1)) % UInt64, magic)
                    lo = P % UInt64
                    hi = (P >> 64) % UInt64
                    pos1 = pos + 1
                    while i <= nlim
                        jn = pos1 + Int((hi >>> (lg & 31)) % Int64)
                        vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + i + 1)), po + jn, nothing)
                        lo, c = Base.add_with_overflow(lo, Clo)
                        hi = hi + Chi + c
                        i += 1
                    end
                else
                    # aarch64: four INDEPENDENT multiply-highs per iteration (`umulh` is
                    # cheap on Apple Silicon; the add/adc carry chain is what regressed the
                    # macos benchmark rows from 1.2–1.4× back to ~0.8×, so the ISA picks
                    # the form — output is bit-identical either way).
                    d = (Int64(1) << _B) - r0
                    while i + 3 <= nlim
                        j1 = pos + _boundary_ceildiv(d, lg, magic)
                        j2 = pos + _boundary_ceildiv(d + (Int64(1) << _B), lg, magic)
                        j3 = pos + _boundary_ceildiv(d + (Int64(2) << _B), lg, magic)
                        j4 = pos + _boundary_ceildiv(d + (Int64(3) << _B), lg, magic)
                        vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + i + 1)), po + j1, nothing)
                        vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + i + 2)), po + j2, nothing)
                        vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + i + 3)), po + j3, nothing)
                        vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + i + 4)), po + j4, nothing)
                        jn = j4
                        d += Int64(4) << _B
                        i += 4
                    end
                    while i <= nlim
                        jn = pos + _boundary_ceildiv(d, lg, magic)
                        vstore(Vec{SW,Int8}(unsafe_load(pp, c0 + i + 1)), po + jn, nothing)
                        d += Int64(1) << _B
                        i += 1
                    end
                end
            else
                # runs can exceed SW (m > SW−2): carry the current chip and add SW-strided
                # interior stores up to the next boundary
                vcur = Vec{SW,Int8}(unsafe_load(pp, c0 + 1))
                jcur = pos
                while i <= nlim
                    jn = pos + _boundary_ceildiv((Int64(i) << _B) - r0, lg, magic)
                    t = jcur + SW
                    while t < jn
                        vstore(vcur, po + t, nothing)
                        t += SW
                    end
                    vcur = Vec{SW,Int8}(unsafe_load(pp, c0 + i + 1))
                    vstore(vcur, po + jn, nothing)
                    jcur = jn
                    i += 1
                end
            end
            if nlim == nc                                   # clean code wrap: re-anchor
                r0 = r0 + Int64(jn - pos) * SN - (Int64(nc) << _B)
                c0 = 0
                pos = jn
                continue
            end
            # wind-down: the remaining chips' stores must clamp at `num`; at most SW−1
            # samples (the final partial store) fall back to scalar byte stores.
            jcur = jn
            vcur = Vec{SW,Int8}(unsafe_load(pp, c0 + i))
            while true
                jnx = i <= nc ? pos + _boundary_ceildiv((Int64(i) << _B) - r0, lg, magic) :
                                pos + _boundary_ceildiv((Int64(nc) << _B) - r0, lg, magic)
                stop = jnx < num ? jnx : num
                t = jcur
                while t + SW <= num && t < stop
                    vstore(vcur, po + t, nothing)
                    t += SW
                end
                if t < stop
                    x = vcur[1]
                    while t < stop
                        unsafe_store!(po, x, t + 1)
                        t += 1
                    end
                end
                stop >= num && return out
                if i >= nc                                  # code wrapped inside wind-down
                    r0 = r0 + Int64(jnx - pos) * SN - (Int64(nc) << _B)
                    c0 = 0
                    pos = jnx
                    break                                   # rare; outer loop re-enters
                end
                vcur = Vec{SW,Int8}(unsafe_load(pp, c0 + i + 1))
                jcur = jnx
                i += 1
            end
        end
        pos < num && _boundary_scalar_tail!(out, padded, L, SN, c0, r0, pos)
    end
    out
end

# Store width: smallest power of two ≥ m+2 (a run is m or m+1 samples, +1 slack so one
# store always covers a run), clamped to [8, 64]. Beyond m = 62 the `EXTRAS` variant adds
# the strided interior stores. Purely a store-width choice — every variant is exact for
# every rate, so this ladder (unlike the old `Val{NI}` one) cannot change the output. The
# SW = 8 rung matters on store-bandwidth-limited cores (CI runners): at m ≈ 5 a 16-wide
# store is 3.3× write amplification vs 8-wide's 1.6× (measured 1.4–1.6× on the ubuntu
# benchmark's GPSL1CA @ 5 MHz row).
@inline function _boundary_dispatch!(out, padded, L::Int, SN::Int64, c0::Int, r0::Int64)
    m = _STEP_DEN ÷ Int(SN)
    if m <= 6
        _boundary_fill!(out, padded, L, SN, c0, r0, Val(8), Val(false))
    elseif m <= 14
        _boundary_fill!(out, padded, L, SN, c0, r0, Val(16), Val(false))
    elseif m <= 30
        _boundary_fill!(out, padded, L, SN, c0, r0, Val(32), Val(false))
    elseif m <= 62
        _boundary_fill!(out, padded, L, SN, c0, r0, Val(64), Val(false))
    else
        _boundary_fill!(out, padded, L, SN, c0, r0, Val(64), Val(true))
    end
    nothing
end

# One-shot boundary fill from an integer chip phase + fractional sub-chip offset `rem0`
# (the same seeding as the permute kernels; byte-identical output).
function _generate_boundary!(out, table::CodeTable, step_num::Int, step_den::Int,
                             phase_offset::Int, rem0::_RemT = _RemT(0))
    c0 = Int(mod(Int64(rem0) >> _B + phase_offset, table.length))
    r0 = Int64(rem0) & Int64(step_den - 1)
    _boundary_dispatch!(out, table.padded, table.length, Int64(step_num), c0, r0)
    out
end

# Pick the boundary fill when there are at least `_boundary_min_m` samples per chip. The
# crossover is where one W-lane permute block (cost ≈ flat per block) equals ~W/m boundary
# stores (cost ≈ flat per chip), so it scales with the backend's block width — wider
# vectors amortise the permute over more samples and push the switch to higher m. Values
# are the measured curve crossings (Zen 5 AVX-512: boundary overtakes the split-constant
# permute at m ≈ 10; AVX2/NEON = 3, matching the old run-fill values — on the CI AVX2
# runner the m ≈ 3.9 rows (GPSL5 @ 40 MHz) lose ~20 % when sent to the permute instead). Near the crossover both kernels are within ~15 % of each other, so the
# exact values are uncritical — and since both kernels are EXACT and N-independent, the
# choice affects only speed, never output (the old N-aware short-fill overload is gone:
# neither kernel has a meaningful setup cost).
@inline _boundary_min_m(::AVX512)   = 10
@inline _boundary_min_m(::AVX2)     = 3
@inline _boundary_min_m(::Neon)     = 3
@inline _boundary_min_m(::Portable) = 2
@inline _use_boundary(step_num::Int, step_den::Int, backend::Backend) =
    step_den ÷ step_num >= _boundary_min_m(backend)
