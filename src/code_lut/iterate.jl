# Allocation-free, value-based code resampling: a loop-invariant `CodeEngine` plus a
# drift-free DDA carried in an isbits `CodeState`, yielding `Vec{W,Int8}` chunks of
# resampled code. Renew the state by value each iteration and fuse it straight into your
# own loop — e.g. pair it with SinCosLUT's value-based carrier engine to build a local
# replica (code · carrier) without ever writing the code to memory:
#
#     ceng = code_engine(ct, sn, sd, Val(1)); cst = code_state(ceng)
#     for _ in 1:nchunks
#         code_vec = code_lookup(ceng, cst)   # Vec{W,Int8} lanes
#         ...                                 # mix with the carrier chunk
#         cst = code_advance(ceng, cst)
#     end

# ---- stateless primitive: hold the padded code + backend, map a phase Vec -> code Vec ----
struct PreparedCode{Backend}
    padded::Vector{Int8}
    backend::Backend
    length::Int
end

"""
    prepare_code(table; backend=default_backend(table)) -> callable

Return a callable `p` such that `p(phase_index::Vec{64,Int32}) -> Vec{64,Int8}`, looking
up the code at each (already mod-`L`) chip index via the sliding 64-chip permute window.
"""
function prepare_code(table::CodeTable; backend::Backend = default_backend(table))
    backend isa Union{AVX2,Neon} && _check_windowed_length(table, backend)
    PreparedCode(table.padded, backend, table.length)
end

@inline (p::PreparedCode{AVX512})(phase::Vec{64,T}) where {T} =
    _window_lookup(AVX512(), p.padded, phase, p.length)
@inline (p::PreparedCode{AVX2})(phase::Vec{32,T}) where {T} =
    _window_lookup(AVX2(), p.padded, phase, p.length)
@inline (p::PreparedCode{Neon})(phase::Vec{16,T}) where {T} =
    _window_lookup(Neon(), p.padded, phase, p.length)
@inline (p::PreparedCode{Portable})(phase::Vec{1,T}) where {T} =
    Vec{1,Int8}(@inbounds p.padded[Int(phase[1]) + 1])

# ─────────────────────────────────────────────────────────────────────────────
# Value-based code engine: a loop-invariant `CodeEngine` (table + baked DDA deltas) plus an
# isbits `CodeState` (the per-stream DDA position) renewed by value every step. The state
# holds no table data, so it stays in registers; fuse it straight into a correlation loop.
# Mirrors SinCosLUT's CarrierEngine/CarrierState. Build one engine per interleave factor `K`
# (`Val(K)`); hold K states `code_state(eng, k-1)` (W apart) and `code_advance` each per step.
#
# Per-backend by necessity — and the gap is large:
#   • AVX-512 (`CodeState512`): the incremental rel/scalar-base form. One `vpermb` over a
#     64-chip window per chunk; the advance only shifts an Int8 `rel` + a scalar `base`.
#     ~50 ps/sample.
#   • AVX2/NEON/portable (`CodeStatePhase`): `vpermb` does not exist, so the chip index is
#     carried as a full per-lane phase vector and the lookup is a windowed `vpshufb`/`tbl`.
#     The DDA advance (not the lookup) dominates and runs ~8× slower (~420 ps/sample) — this
#     gap is fundamental to the ISA (no 64-wide cross-lane byte permute), not a tuning miss.
# ─────────────────────────────────────────────────────────────────────────────

# ---- AVX-512 engine + state (incremental rel / scalar base) ----
struct CodeEngine512
    padded::Vector{Int8}
    step_num::Int
    step_den::Int
    phase_offset::Int
    modulus::_RemT
    frac_stride::_RemT     # remainder advance over the baked K·W stride
    whole_stride::Int      # whole chips advanced over the stride (base does mod L)
    L::Int
end
struct CodeState512
    rel::Vec{64,Int8}
    rem::Vec{64,_RemT}
    base::Int
end

# ---- phase engine + state (AVX2 / NEON / portable) ----
struct CodeEnginePhase{W,T,Prep}
    prepared::Prep
    step_num::Int
    step_den::Int
    phase_offset::Int
    modulus::_RemT
    frac_stride::_RemT
    whole_stride::T        # whole chips per stride, already reduced mod L (< L)
    L::Int
end
struct CodeStatePhase{W,T}
    phase::Vec{W,T}
    rem::Vec{W,_RemT}
end

"""
    code_engine(table, step_num, step_den, Val(K); phase=0, backend=…) -> CodeEngine
    code_engine(table, chips_per_sample::Real, Val(K); …)
    code_engine(table, Val(K); code_frequency, sampling_frequency, …)

Build the loop-invariant code engine for `table`, baking the drift-free fixed-point DDA
deltas for a `K`-way interleaved loop (stride `K·W`, `W` = backend SIMD width). Pair it with
`K` states `code_state(eng, k)` (`k = 0..K-1`, `W` samples apart) and drive each with
[`code_lookup`](@ref) / [`code_advance`](@ref). `phase` is an integer chip phase offset.
Same oversampling / denominator requirements as `generate_code!`.
"""
function code_engine(table::CodeTable, step_num::Integer, step_den::Integer, ::Val{K};
                     phase::Integer = 0, backend::Backend = default_backend(table)) where {K}
    (0 < step_den ≤ _STEP_DEN) ||
        throw(ArgumentError("need 0 < step_denominator ≤ 2^$_B"))
    (0 < step_num ≤ step_den) ||
        throw(ArgumentError("need 0 < step_numerator ≤ step_denominator (must oversample)"))
    _make_engine(table, Int(step_num), Int(step_den), Int(phase), Val(K), backend, _vwidth(backend))
end
function code_engine(table::CodeTable, cps::Real, vk::Val; kw...)
    sn, sd = _fixed_point_step(cps)
    code_engine(table, sn, sd, vk; kw...)
end
function code_engine(table::CodeTable, vk::Val; code_frequency::Real, sampling_frequency::Real, kw...)
    code_engine(table, chips_per_sample(code_frequency, sampling_frequency), vk; kw...)
end

function _make_engine(table::CodeTable, sn, sd, ph, ::Val{K}, ::AVX512, ::Val{64}) where {K}
    W = 64; stride = K * W
    CodeEngine512(table.padded, sn, sd, ph, _RemT(sd),
        _RemT(mod(stride * sn, sd)), div(stride * sn, sd) % table.length, table.length)
end
function _make_engine(table::CodeTable, sn, sd, ph, ::Val{K}, backend, ::Val{W}) where {K,W}
    backend isa Union{AVX2,Neon} && _check_windowed_length(table, backend)
    stride = K * W; L = table.length; T = _phase_type(L)
    prepared = prepare_code(table; backend = backend)
    CodeEnginePhase{W,T,typeof(prepared)}(prepared, sn, sd, ph, _RemT(sd),
        _RemT(mod(stride * sn, sd)), T(div(stride * sn, sd) % L), L)
end

"""
    code_state(eng::CodeEngine, stream::Integer = 0) -> CodeState

Initial DDA state for the `stream`-th interleaved lane (first sample at `stream·W`).
"""
@inline function code_state(eng::CodeEngine512, stream::Integer = 0)
    rel, rem, base = _init_rel(Val(64), eng.step_num, eng.step_den, eng.L, stream * 64, eng.phase_offset)
    CodeState512(rel, rem, base)
end
@inline function code_state(eng::CodeEnginePhase{W,T}, stream::Integer = 0) where {W,T}
    p, r = _init_state(Val(W), eng.step_num, eng.step_den, eng.L, stream * W, eng.phase_offset, T)
    CodeStatePhase{W,T}(p, r)
end

"""
    code_lookup(eng::CodeEngine, st::CodeState) -> Vec{W,Int8}

The resampled ±1 code chunk at `st`'s current chip position. Pure read; does not advance.
"""
@inline code_lookup(eng::CodeEngine512, st::CodeState512) =
    _rel_lookup(AVX512(), eng.padded, st.rel, st.base)
@inline code_lookup(eng::CodeEnginePhase{W}, st::CodeStatePhase{W}) where {W} =
    eng.prepared(st.phase)

# phase DDA step (carry-propagating fixed-point advance), returning the new (phase, rem)
@inline function _advance_phase(phase::Vec{W,T}, rem::Vec{W,_RemT}, frac::_RemT,
                                whole::T, modulus::_RemT, Lc::T) where {W,T}
    a = rem + frac
    carry = a >= modulus
    rem2 = vifelse(carry, a - modulus, a)
    p = vifelse(carry, phase + whole + one(T), phase + whole)
    p = vifelse(p >= Lc, p - Lc, p)
    (p, rem2)
end

"""
    code_advance(eng::CodeEngine, st::CodeState) -> CodeState

Advance the stream by one engine stride (`K·W` samples), returning a new immutable state.
"""
@inline function code_advance(eng::CodeEngine512, st::CodeState512)
    rel, rem, base = _advance_window_512(st.rel, st.rem, st.base,
                                         eng.frac_stride, eng.whole_stride, eng.modulus, eng.L)
    CodeState512(rel, rem, base)
end
@inline function code_advance(eng::CodeEnginePhase{W,T}, st::CodeStatePhase{W,T}) where {W,T}
    p, r = _advance_phase(st.phase, st.rem, eng.frac_stride, eng.whole_stride, eng.modulus, T(eng.L))
    CodeStatePhase{W,T}(p, r)
end

"""
    code_width(eng::CodeEngine) -> Int

SIMD lane count `W` of the engine (samples per chunk).
"""
@inline code_width(::CodeEngine512) = 64
@inline code_width(::CodeEnginePhase{W}) where {W} = W
