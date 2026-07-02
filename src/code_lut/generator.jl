# Value-based, *continuing* code generator (module-internal engine for the fill `gen_code!`).
#
# The one-shot `generate_code!` pays a full DDA setup (fixed-point step + four `_init_rel`
# streams) on every call — overhead that dwarfs the ~120 ns generation at 1 ms sizes. A
# `FillEngine` pays that setup ONCE (at construction) and then *continues* from an explicit,
# isbits state across arbitrarily-sized fills, so concatenating consecutive fills equals one
# big generation. It also vectorises the non-stride-aligned tail (full 4-way strides → W-wide
# single-window steps → only the final < W samples scalar).
#
# State is explicit and immutable: a single canonical "stream 0" `(rel, rem, base)` at the
# current absolute sample offset (AVX-512 rel/scalar-base form, reusing `CodeState512`), or
# the phase-based equivalent (`CodeStatePhase`) for AVX2/Portable, or `(c, r)`
# (`BoundaryState`) for the high-oversampling boundary-fill path. `fill_continue!(out, eng, st)`
# returns the state advanced by exactly `length(out)`, so the caller threads it — no mutation,
# no hidden state, and an isbits state means the steady-state fill is allocation-free.
#
# The 4-way streams used for throughput (AVX-512) are *spawned* from the carried stream 0 by
# three cheap W-window advances (≈20 ns, byte-identical to `_init_rel`) at the start of each
# fill — not stored in the state.

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
    # `whole_W` is pre-reduced mod L by the caller (FillEngine512 / _make_engine), so
    # `base + whole_W + h ∈ [0, 2L-1]` and a single conditional subtract replaces the idiv.
    base2 = _wrapL(base + whole_W + h, L)
    (rel2, rem2, base2)
end

# ---- value-based, continuing fill engine over the AVX-512 rel/scalar-base DDA ----
# Immutable loop-invariant config: the padded table, the fixed-point step, per-W and
# per-stride DDA deltas (computed once at construction; no rationalize, no idiv per fill), and
# the initial phase so `fill_state` can seed the canonical stream-0 state. The four
# interleaved streams (W apart) used for throughput are spawned from the carried stream-0 at
# the start of each fill. After a fill of M samples the returned state is stream-0 at absolute
# offset M, so the next call resumes seamlessly.
struct FillEngine512
    padded::Vector{Int8}
    step_num::Int
    step_den::Int
    L::Int
    frac_W::_RemT          # remainder advance over W samples
    whole_W::Int           # whole chips advanced over W samples
    frac_stride::_RemT     # remainder advance over a 4W stride
    whole_stride::Int      # whole chips advanced over a 4W stride
    modulus::_RemT
    phase_offset::Int      # initial integer chip phase (for fill_state)
    rem0::_RemT            # initial fractional sub-chip offset (for fill_state)
end

function FillEngine512(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int,
                       rem0::_RemT = _RemT(0))
    L = table.length; W = 64
    FillEngine512(table.padded, step_num, step_den, L,
        _RemT(mod(W * step_num, step_den)), div(W * step_num, step_den) % L,
        _RemT(mod(4W * step_num, step_den)), div(4W * step_num, step_den) % L, _RemT(step_den),
        phase_offset, rem0)
end

# Initial canonical stream-0 state (reuses `CodeState512` from iterate.jl).
@inline function fill_state(eng::FillEngine512)
    rel, rem, base = _init_rel(Val(64), eng.step_num, eng.step_den, eng.L, 0, eng.phase_offset, eng.rem0)
    CodeState512(rel, rem, base)
end

# Advance one stream by a full 4W stride (used for the lockstep bulk advance).
@inline _advance_stride_512(rel, rem, base, fs::_RemT, ws::Int, modulus::_RemT, L::Int) =
    _advance_window_512(rel, rem, base, fs, ws, modulus, L)

# Fill out[1:end] with the next length(out) samples, continuing from `st` (canonical stream-0)
# and returning the state advanced by exactly length(out). `out` is Int8.
function fill_continue!(out::AbstractVector{<:Integer}, eng::FillEngine512, st::CodeState512)
    W = 64; stride = 4W
    padded = eng.padded; L = eng.L
    frac_W = eng.frac_W; whole_W = eng.whole_W; modulus = eng.modulus
    fstr = eng.frac_stride; wstr = eng.whole_stride
    # Spawn the four interleaved streams (W apart) from the carried stream-0.
    rel1 = st.rel; rem1 = st.rem; b1 = st.base
    rel2, rem2, b2 = _advance_window_512(rel1, rem1, b1, frac_W, whole_W, modulus, L)
    rel3, rem3, b3 = _advance_window_512(rel2, rem2, b2, frac_W, whole_W, modulus, L)
    rel4, rem4, b4 = _advance_window_512(rel3, rem3, b3, frac_W, whole_W, modulus, L)
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
    # then advance stream-0 by rem_s so the returned state is stream-0 at absolute offset num.
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
        # advance stream-0 by rem_s = nwin·W + ntail samples (the only scalar advance).
        rel1, rem1, b1 = _advance_by_512(rel1, rem1, b1, nwin, ntail, frac_W, whole_W, eng.step_num, eng.step_den, modulus, L)
    end
    # stream-0 is now at absolute offset num (after both the bulk and leftover advances).
    return CodeState512(rel1, rem1, b1)
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

# ---- AVX2 / Portable continuing fill engine (phase-based, single stream, W at a time) ----
# Keeps the phase-based DDA (the rel/scalar-base trick regressed AVX2). Correct continuation
# is the same scheme: carry one canonical phase stream (`CodeStatePhase`), emit W-wide
# windows, scalar tail. The engine holds the loop-invariant config + initial phase.
struct FillEnginePhase{W,T,Prep}
    prepared::Prep
    step_num::Int
    step_den::Int
    whole_W::T
    frac_W::_RemT
    modulus::_RemT
    code_length::T
    phase_offset::Int      # initial integer chip phase (for fill_state)
    rem0::_RemT            # initial fractional sub-chip offset (for fill_state)
end

# Branch on the table length so each arm passes a *literal* phase type (Int16/Int32) into the
# `_build_phase_engine` barrier. `_phase_type(L)` returns a runtime type *value*; using it
# directly made construction a dynamic dispatch and boxed the whole engine (the one-shot /
# threaded allocation hot spot). The two arms give inference a small 2-way Union instead.
function FillEnginePhase(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int,
                         backend::Backend, vw::Val{W}, rem0::_RemT = _RemT(0)) where {W}
    table.length <= typemax(Int16) ?
        _build_phase_engine(table, step_num, step_den, phase_offset, backend, vw, Int16, rem0) :
        _build_phase_engine(table, step_num, step_den, phase_offset, backend, vw, Int32, rem0)
end

@inline function _build_phase_engine(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int,
                                     backend::Backend, ::Val{W}, ::Type{T}, rem0::_RemT = _RemT(0)) where {W,T}
    L = table.length
    prepared = prepare_code(table; backend = backend)
    FillEnginePhase{W,T,typeof(prepared)}(prepared, step_num, step_den,
        T(div(W * step_num, step_den) % L), _RemT(mod(W * step_num, step_den)),
        _RemT(step_den), T(L), phase_offset, rem0)
end

@inline function fill_state(eng::FillEnginePhase{W,T}) where {W,T}
    phase, rem = _init_state(Val(W), eng.step_num, eng.step_den, Int(eng.code_length), 0, eng.phase_offset, T, eng.rem0)
    CodeStatePhase{W,T}(phase, rem)
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

function fill_continue!(out::AbstractVector{<:Integer}, eng::FillEnginePhase{W,T}, st::CodeStatePhase{W,T}) where {W,T}
    p = eng.prepared; whole_W = eng.whole_W; frac_W = eng.frac_W; modulus = eng.modulus; Lc = eng.code_length
    phase = st.phase; rem = st.rem
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
        phase, rem = _advance_scalar_phase(phase, rem, ntail, eng.step_num, eng.step_den, modulus, Lc)
        pos += ntail
    end
    return CodeStatePhase{W,T}(phase, rem)
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

# ---- continuing boundary fill engine (high-oversampling fast path) ----
# The continuing counterpart of `_generate_boundary!`: at high oversampling splat-fills one
# chip run per store instead of paying a permute per W samples (see kernel.jl). Because the
# boundary arithmetic is EXACT, the carried state is just `(c, r)` — the table chip index
# and the fractional chip phase — and continuation is byte-identical to one big fill for any
# chunking (the old run-fill's `(c, acc, pos)` approximate-DDA state and its documented
# multi-fill drift are gone). The store width `SW` / long-run `EXTRAS` flag are type
# parameters so `fill_continue!` calls the specialised kernel statically (0-alloc); like the
# old `Val{NI}` they are fixed at construction, but they are a pure store-width choice — any
# variant is exact for any rate.
struct BoundaryState
    c::Int                    # current table chip index
    r::_RemT                  # fractional chip phase (< 2^_B)
end

struct FillEngineBoundary{SW,EXTRAS}
    padded::Vector{Int8}
    L::Int
    step_num::Int
    c0::Int             # initial table chip index (for fill_state)
    r0::_RemT           # initial fractional chip phase (for fill_state)
end

function FillEngineBoundary(table::CodeTable, step_num::Int, step_den::Int, phase_offset::Int,
                            rem0::_RemT = _RemT(0))
    c0 = Int(mod(Int64(rem0) >> _B + phase_offset, table.length))
    r0 = _RemT(Int64(rem0) & Int64(step_den - 1))
    m = _STEP_DEN ÷ step_num
    if m <= 14
        FillEngineBoundary{16,false}(table.padded, table.length, step_num, c0, r0)
    elseif m <= 30
        FillEngineBoundary{32,false}(table.padded, table.length, step_num, c0, r0)
    elseif m <= 62
        FillEngineBoundary{64,false}(table.padded, table.length, step_num, c0, r0)
    else
        FillEngineBoundary{64,true}(table.padded, table.length, step_num, c0, r0)
    end
end

@inline fill_state(eng::FillEngineBoundary) = BoundaryState(eng.c0, eng.r0)

function fill_continue!(out::AbstractVector{<:Integer}, eng::FillEngineBoundary{SW,EXTRAS},
                        st::BoundaryState) where {SW,EXTRAS}
    SN = Int64(eng.step_num)
    _boundary_fill!(out, eng.padded, eng.L, SN, st.c, Int64(st.r), Val(SW), Val(EXTRAS))
    # advance the state by exactly length(out) samples, arithmetically (exact)
    tot = Int64(st.r) + Int64(length(out)) * SN
    BoundaryState(Int(mod(Int64(st.c) + (tot >> _B), eng.L)), _RemT(tot & Int64(_STEP_DEN - 1)))
end

# Union over the backend fill engines (Portable is the W=1 degenerate phase engine: window
# emit + advance are scalar but correct).
const FillEngineAny = Union{FillEngine512,FillEnginePhase,FillEngineBoundary}

# ---- factory: build the right backend fill engine from a step ratio + phase ----
function make_fill_engine(table::CodeTable, step_num::Integer, step_den::Integer;
                          phase::Integer = 0, rem0::Integer = 0,
                          backend::Backend = default_backend(table))
    (0 < step_den ≤ _STEP_DEN) ||
        throw(ArgumentError("need 0 < step_denominator ≤ 2^$_B"))
    (0 < step_num ≤ step_den) ||
        throw(ArgumentError("need 0 < step_numerator ≤ step_denominator (must oversample)"))
    (0 ≤ rem0 < step_den) ||
        throw(ArgumentError("need 0 ≤ rem0 < step_denominator (fractional sub-chip offset)"))
    sn = Int(step_num); sd = Int(step_den); ph = Int(phase); r0 = _RemT(rem0)
    if _use_boundary(sn, sd, backend)
        # High oversampling: splat-fill chip runs (see `_generate_boundary!`) instead of
        # paying a permute per window. No windowed-length requirement (no in-register window).
        FillEngineBoundary(table, sn, sd, ph, r0)
    elseif backend isa AVX512
        FillEngine512(table, sn, sd, ph, r0)
    elseif backend isa Union{AVX2,Neon}
        _check_windowed_length(table, backend)
        FillEnginePhase(table, sn, sd, ph, backend, _vwidth(backend), r0)
    else
        FillEnginePhase(table, sn, sd, ph, backend, Val(1), r0)
    end
end
function make_fill_engine(table::CodeTable, cps::Real; kw...)
    sn, sd = _fixed_point_step(cps)
    make_fill_engine(table, sn, sd; kw...)
end
function make_fill_engine(table::CodeTable; code_frequency::Real, sampling_frequency::Real, kw...)
    make_fill_engine(table, chips_per_sample(code_frequency, sampling_frequency); kw...)
end
