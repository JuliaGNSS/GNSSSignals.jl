# Stateful, *continuing* code generator (module-internal engine for CodeGeneratorLUT).
#
# The one-shot `generate_code!` pays a full DDA setup (fixed-point step + four `_init_rel`
# streams) on every call — overhead that dwarfs the ~120 ns generation at 1 ms sizes. This
# engine pays that setup ONCE (at construction) and then *continues* from the carried DDA
# state across arbitrarily-sized fills, so concatenating consecutive fills equals one big
# generation. It also vectorises the non-stride-aligned tail (full 4-way strides → W-wide
# single-window steps → only the final < W samples scalar).
#
# State is a single canonical "stream 0" `(rel, rem, base)` at the current absolute sample
# offset (AVX-512 rel/scalar-base form), or the phase-based equivalent for AVX2/Portable.
# The 4-way streams used for throughput are *spawned* from stream 0 by three cheap W-window
# advances (≈20 ns, byte-identical to `_init_rel`) — not re-initialised.

# ---- per-W single-stream advance (AVX-512 rel/scalar-base form) ----
# Advance one W=64 stream by exactly W samples. Mirrors `CodeState512`'s state update
# but with per-W deltas (frac_W / whole_W) passed in, so it works for any window position.
@inline function _advance_window_512(rel::Vec{64,Int8}, rem::Vec{64,_RemT}, base::Int,
                                     frac_W::_RemT, whole_W::Int, modulus::_RemT, L::Int)
    rem2 = rem + frac_W
    carry = rem2 >= modulus
    rem2 = vifelse(carry, rem2 - modulus, rem2)
    h = Int(carry[1])
    rel2 = rel + vifelse(carry, one(Vec{64,Int8}), zero(Vec{64,Int8})) - Int8(h)
    # `whole_W` is pre-reduced mod L by the caller (CodeGenerator512 / _make_engine), so
    # `base + whole_W + h ∈ [0, 2L-1]` and a single conditional subtract replaces the idiv.
    base2 = _wrapL(base + whole_W + h, L)
    (rel2, rem2, base2)
end

# ---- mutable, continuing generator over the AVX-512 rel/scalar-base DDA ----
# Carries the canonical stream-0 state at the current absolute offset, plus per-W and
# per-stride DDA deltas (computed once at construction; no rationalize, no idiv per fill).
# Carries the FOUR interleaved streams (W apart) plus per-W / per-stride DDA deltas
# (computed once at construction; no rationalize, no idiv per fill). The bulk loop advances
# all four in lockstep by a full stride (the exact one-shot kernel loop). After a fill of M
# samples the four streams are repositioned to absolute offsets M+{0,W,2W,3W} so the next
# call resumes seamlessly.
mutable struct CodeGenerator512
    const padded::Vector{Int8}
    const step_num::Int
    const step_den::Int
    const L::Int
    const frac_W::_RemT        # remainder advance over W samples
    const whole_W::Int         # whole chips advanced over W samples
    const frac_stride::_RemT   # remainder advance over a 4W stride
    const whole_stride::Int    # whole chips advanced over a 4W stride
    const modulus::_RemT
    rel1::Vec{64,Int8}; rem1::Vec{64,_RemT}; b1::Int
    rel2::Vec{64,Int8}; rem2::Vec{64,_RemT}; b2::Int
    rel3::Vec{64,Int8}; rem3::Vec{64,_RemT}; b3::Int
    rel4::Vec{64,Int8}; rem4::Vec{64,_RemT}; b4::Int
end

function CodeGenerator512(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int)
    L = table.length; W = 64
    rel1, rem1, b1 = _init_rel(Val(W), step_num, step_den, L, 0,  phase_offset)
    rel2, rem2, b2 = _init_rel(Val(W), step_num, step_den, L, W,  phase_offset)
    rel3, rem3, b3 = _init_rel(Val(W), step_num, step_den, L, 2W, phase_offset)
    rel4, rem4, b4 = _init_rel(Val(W), step_num, step_den, L, 3W, phase_offset)
    CodeGenerator512(table.padded, step_num, step_den, L,
        _RemT(mod(W * step_num, step_den)), div(W * step_num, step_den) % L,
        _RemT(mod(4W * step_num, step_den)), div(4W * step_num, step_den) % L, _RemT(step_den),
        rel1, rem1, b1, rel2, rem2, b2, rel3, rem3, b3, rel4, rem4, b4)
end

# Advance one stream by a full 4W stride (used for the lockstep bulk advance).
@inline _advance_stride_512(rel, rem, base, fs::_RemT, ws::Int, modulus::_RemT, L::Int) =
    _advance_window_512(rel, rem, base, fs, ws, modulus, L)

# Fill out[1:end] with the next length(out) samples, continuing from `g`'s state and saving
# the state advanced by exactly length(out) back into `g`. `out` is Int8.
function fill_continue!(out::AbstractVector{<:Integer}, g::CodeGenerator512)
    W = 64; stride = 4W
    padded = g.padded; L = g.L
    frac_W = g.frac_W; whole_W = g.whole_W; modulus = g.modulus
    fstr = g.frac_stride; wstr = g.whole_stride
    rel1 = g.rel1; rem1 = g.rem1; b1 = g.b1
    rel2 = g.rel2; rem2 = g.rem2; b2 = g.b2
    rel3 = g.rel3; rem3 = g.rem3; b3 = g.b3
    rel4 = g.rel4; rem4 = g.rem4; b4 = g.b4
    num = length(out)
    nfull = num ÷ stride                       # full 4-way strides
    rem_s = num - nfull * stride               # leftover samples (0..255)
    nwin  = rem_s ÷ W                           # leftover full W-windows (0..3)
    ntail = rem_s - nwin * W                     # final < W samples

    pos = 0
    @inbounds for _ in 1:nfull
        out[VecRange{W}(pos + 1)]      = _rel_lookup(AVX512(), padded, rel1, b1)
        out[VecRange{W}(pos + W + 1)]  = _rel_lookup(AVX512(), padded, rel2, b2)
        out[VecRange{W}(pos + 2W + 1)] = _rel_lookup(AVX512(), padded, rel3, b3)
        out[VecRange{W}(pos + 3W + 1)] = _rel_lookup(AVX512(), padded, rel4, b4)
        rel1, rem1, b1 = _advance_stride_512(rel1, rem1, b1, fstr, wstr, modulus, L)
        rel2, rem2, b2 = _advance_stride_512(rel2, rem2, b2, fstr, wstr, modulus, L)
        rel3, rem3, b3 = _advance_stride_512(rel3, rem3, b3, fstr, wstr, modulus, L)
        rel4, rem4, b4 = _advance_stride_512(rel4, rem4, b4, fstr, wstr, modulus, L)
        pos += stride
    end
    # Leftover (< stride): emit from the streams in order (stream k covers samples [k·W,..));
    # then reposition all four streams to absolute offsets +rem_s so the next call resumes.
    if rem_s > 0
        streams = ((rel1, b1), (rel2, b2), (rel3, b3), (rel4, b4))
        # full W-windows of the leftover
        @inbounds for k in 1:nwin
            rk, bk = streams[k]
            out[VecRange{W}(pos + 1)] = _rel_lookup(AVX512(), padded, rk, bk)
            pos += W
        end
        if ntail > 0
            rk, bk = streams[nwin + 1]
            win = _rel_lookup(AVX512(), padded, rk, bk)
            @inbounds for j in 1:ntail
                out[pos + j] = win[j]
            end
            pos += ntail
        end
        # reposition: advance stream 1 by rem_s = nwin·W + ntail samples (the only scalar
        # advance), then respawn streams 2..4 as stream 1 + W, 2W, 3W (cheap W-advances).
        rel1, rem1, b1 = _advance_by_512(rel1, rem1, b1, nwin, ntail, frac_W, whole_W, g.step_num, g.step_den, modulus, L)
        rel2, rem2, b2 = _advance_window_512(rel1, rem1, b1, frac_W, whole_W, modulus, L)
        rel3, rem3, b3 = _advance_window_512(rel2, rem2, b2, frac_W, whole_W, modulus, L)
        rel4, rem4, b4 = _advance_window_512(rel3, rem3, b3, frac_W, whole_W, modulus, L)
    end
    g.rel1 = rel1; g.rem1 = rem1; g.b1 = b1
    g.rel2 = rel2; g.rem2 = rem2; g.b2 = b2
    g.rel3 = rel3; g.rem3 = rem3; g.b3 = b3
    g.rel4 = rel4; g.rem4 = rem4; g.b4 = b4
    out
end

# Advance one stream by nwin full W-windows + ntail (< W) scalar samples.
@inline function _advance_by_512(rel, rem, base, nwin::Int, ntail::Int, frac_W::_RemT,
                                 whole_W::Int, step_num::Int, step_den::Int, modulus::_RemT, L::Int)
    rl = rel; r = rem; b = base
    @inbounds for _ in 1:nwin
        rl, r, b = _advance_window_512(rl, r, b, frac_W, whole_W, modulus, L)
    end
    ntail > 0 && ((rl, r, b) = _advance_scalar_512(rl, r, b, ntail, step_num, step_den, modulus, L))
    (rl, r, b)
end

# Advance one W=64 stream by `m` (0 ≤ m < W) samples, sample by sample. Used only for the
# final sub-window tail, so the scalar cost is bounded (< 64 iterations) and off the hot path.
@inline function _advance_scalar_512(rel::Vec{64,Int8}, rem::Vec{64,_RemT}, base::Int,
                                     m::Int, step_num::Int, step_den::Int, modulus::_RemT, L::Int)
    frac1 = _RemT(step_num % step_den); whole1 = step_num ÷ step_den
    r = rem; rl = rel; b = base
    @inbounds for _ in 1:m
        r2 = r + frac1
        carry = r2 >= modulus
        r = vifelse(carry, r2 - modulus, r2)
        h = Int(carry[1])
        rl = rl + vifelse(carry, one(Vec{64,Int8}), zero(Vec{64,Int8})) - Int8(h)
        b = _wrapL(b + whole1 + h, L)   # whole1 ∈ {0,1} < L ⇒ sum < 2L
    end
    (rl, r, b)
end

# ---- AVX2 / Portable continuing generator (phase-based, single stream, W at a time) ----
# Keeps the phase-based DDA (the rel/scalar-base trick regressed AVX2). Correct continuation
# is the same scheme: carry one canonical phase stream, emit W-wide windows, scalar tail.
mutable struct CodeGeneratorPhase{W,T,Prep}
    const prepared::Prep
    const step_num::Int
    const step_den::Int
    const whole_W::T
    const frac_W::_RemT
    const modulus::_RemT
    const code_length::T
    phase::Vec{W,T}
    rem::Vec{W,_RemT}
end

# Branch on the table length so each arm passes a *literal* phase type (Int16/Int32) into the
# `_build_phase` barrier. `_phase_type(L)` returns a runtime type *value*; using it directly
# made `_init_state` a dynamic dispatch and boxed the whole construction (the one-shot /
# threaded allocation hot spot). The two arms give inference a small 2-way Union instead.
function CodeGeneratorPhase(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int,
                            backend::Backend, vw::Val{W}) where {W}
    table.length <= typemax(Int16) ?
        _build_phase(table, step_num, step_den, phase_offset, backend, vw, Int16) :
        _build_phase(table, step_num, step_den, phase_offset, backend, vw, Int32)
end

@inline function _build_phase(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int,
                              backend::Backend, ::Val{W}, ::Type{T}) where {W,T}
    L = table.length
    phase, rem = _init_state(Val(W), step_num, step_den, L, 0, phase_offset, T)
    prepared = prepare_code(table; backend = backend)
    CodeGeneratorPhase{W,T,typeof(prepared)}(prepared, step_num, step_den,
        T(div(W * step_num, step_den) % L), _RemT(mod(W * step_num, step_den)),
        _RemT(step_den), T(L), phase, rem)
end

@inline function _advance_window_phase(phase::Vec{W,T}, rem::Vec{W,_RemT}, whole_W::T,
                                       frac_W::_RemT, modulus::_RemT, Lc::T) where {W,T}
    rem2 = rem + frac_W
    carry = rem2 >= modulus
    rem2 = vifelse(carry, rem2 - modulus, rem2)
    phase2 = vifelse(carry, phase + whole_W + one(T), phase + whole_W)
    phase2 = vifelse(phase2 >= Lc, phase2 - Lc, phase2)
    (phase2, rem2)
end

function fill_continue!(out::AbstractVector{<:Integer}, g::CodeGeneratorPhase{W,T}) where {W,T}
    p = g.prepared; whole_W = g.whole_W; frac_W = g.frac_W; modulus = g.modulus; Lc = g.code_length
    phase = g.phase; rem = g.rem
    num = length(out)
    nwin = num ÷ W; ntail = num - nwin * W
    pos = 0
    @inbounds for _ in 1:nwin
        out[VecRange{W}(pos + 1)] = p(phase)
        phase, rem = _advance_window_phase(phase, rem, whole_W, frac_W, modulus, Lc)
        pos += W
    end
    if ntail > 0
        win = p(phase)
        @inbounds for j in 1:ntail
            out[pos + j] = win[j]
        end
        phase, rem = _advance_scalar_phase(phase, rem, ntail, g.step_num, g.step_den, modulus, Lc)
        pos += ntail
    end
    g.phase = phase; g.rem = rem
    out
end

@inline function _advance_scalar_phase(phase::Vec{W,T}, rem::Vec{W,_RemT}, m::Int,
                                       step_num::Int, step_den::Int, modulus::_RemT, Lc::T) where {W,T}
    frac1 = _RemT(step_num % step_den); whole1 = T(step_num ÷ step_den)
    ph = phase; r = rem
    @inbounds for _ in 1:m
        r2 = r + frac1
        carry = r2 >= modulus
        r = vifelse(carry, r2 - modulus, r2)
        ph = vifelse(carry, ph + whole1 + one(T), ph + whole1)
        ph = vifelse(ph >= Lc, ph - Lc, ph)
    end
    (ph, r)
end

# ---- continuing run-fill generator (high-oversampling fast path) ----
# The continuing counterpart of `_generate_runfill!`: at high oversampling broadcast-fills
# runs of identical baked chips (store-bandwidth bound) instead of one permute per W samples.
# Effectively byte-identical to the permute engines (the approximate samples-per-chip DDA in
# `_runfill_core!` matches the exact `floor(step_num·n/2^_B)` boundaries for any single fill,
# drifting ≤ a couple of samples only over ~10⁵ continued chips). State carried across fills
# is `(c, acc, pos)`: the table chip index, the chip-start fractional phase `acc`, and `pos`
# = samples of the current chip already emitted (so a fill may resume mid-run).
#
# `NI` (the padded fixed inner-store count) is a *type* parameter so `fill_continue!` calls
# the `Val{NI}`-specialised kernel statically — a runtime `Val(ni)` would box and allocate
# every call, breaking the 0-alloc steady-state guarantee. `NI == 0` is the sentinel for the
# generic (`ni > _RUNFILL_MAX_NI`) kernel, whose runtime count lives in the `ni` field.
mutable struct CodeGeneratorRunFill{NI}
    const chips::Vector{Int8}
    const L::Int
    const ff::Int             # samples-per-chip in _RUNFILL_FP fixed point
    const ni::Int             # padded inner count (used by the NI==0 generic path)
    c::Int                    # current table chip index
    acc::Int                  # chip-start fractional phase (delta & mask)
    pos::Int                  # samples of the current chip already emitted
end

function CodeGeneratorRunFill(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int)
    ni = _runfill_ni(step_num)
    NI = ni <= _RUNFILL_MAX_NI ? ni : 0      # 0 ⇒ generic kernel (runtime `ni`)
    CodeGeneratorRunFill{NI}(table.chips, table.length, _runfill_freqfix(step_num),
                            ni, mod(phase_offset, table.length), _RUNFILL_MASK, 0)
end

function fill_continue!(out::AbstractVector{<:Integer}, g::CodeGeneratorRunFill{NI}) where {NI}
    g.c, g.acc, g.pos = NI == 0 ?
        _runfill_seg_generic!(out, g.chips, g.L, g.ff, g.ni, g.c, g.acc, g.pos) :
        _runfill_seg!(out, g.chips, g.L, g.ff, g.c, g.acc, g.pos, Val(NI))
    out
end

# Portable (W=1): degenerate phase generator. Window emit + advance are scalar but correct.
const CodeGeneratorAny = Union{CodeGenerator512,CodeGeneratorPhase,CodeGeneratorRunFill}

# ---- factory: build the right backend generator from a step ratio + phase ----
function make_generator(table::CodeTable, step_num::Integer, step_den::Integer;
                        phase::Integer = 0, backend::Backend = default_backend(table))
    (0 < step_den ≤ _STEP_DEN) ||
        throw(ArgumentError("need 0 < step_denominator ≤ 2^$_B"))
    (0 < step_num ≤ step_den) ||
        throw(ArgumentError("need 0 < step_numerator ≤ step_denominator (must oversample)"))
    sn = Int(step_num); sd = Int(step_den); ph = Int(phase)
    if _use_runfill(sn, sd, backend)
        # High oversampling: broadcast-fill runs (see `_generate_runfill!`) instead of paying
        # a permute per window. Same windowed-length requirement as the AVX2/NEON permute.
        backend isa Union{AVX2,Neon} && _check_windowed_length(table, backend)
        CodeGeneratorRunFill(table, sn, sd, ph)
    elseif backend isa AVX512
        CodeGenerator512(table, sn, sd, ph)
    elseif backend isa Union{AVX2,Neon}
        _check_windowed_length(table, backend)
        CodeGeneratorPhase(table, sn, sd, ph, backend, _vwidth(backend))
    else
        CodeGeneratorPhase(table, sn, sd, ph, backend, Val(1))
    end
end
function make_generator(table::CodeTable, cps::Real; kw...)
    sn, sd = _fixed_point_step(cps)
    make_generator(table, sn, sd; kw...)
end
function make_generator(table::CodeTable; code_frequency::Real, sampling_frequency::Real, kw...)
    make_generator(table, chips_per_sample(code_frequency, sampling_frequency); kw...)
end
