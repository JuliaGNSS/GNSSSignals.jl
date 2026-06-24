# Allocation-free iterator: a drift-free DDA carried in the (isbits) iteration state,
# yielding `Vec{W,Int8}` chunks of resampled code. Fuse it straight into your own loop —
# e.g. zip it with SinCosLUT.generate_carrier to build a local replica (code · carrier)
# without ever writing the code to memory:
#
#     for (code_vec, (sin_vec, cos_vec)) in zip(generate_code(ct, cps_code, n),
#                                                generate_carrier(sct, cps_carrier, n))
#         replica_i = code_vec * sin_vec      # Vec{W,Int8} lanes
#         ...
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

# ---- stateful iterator: drift-free phase, yields Vec{W,Int8} ----
# `T` is the phase type (Int16 for short tables, Int32 for long ones — see _phase_type).
struct CodeIterator{W,T,Prep}
    prepared::Prep
    phase_init::Vec{W,T}
    remainder_init::Vec{W,_RemT}
    whole_step::T
    frac_step::_RemT
    modulus::_RemT
    code_length::T
    num_chunks::Int
end

"""
    generate_code(table, step_numerator, step_denominator, num_samples; phase=0, backend=…)
    generate_code(table, chips_per_sample::Real, num_samples;            phase=0, backend=…)
    generate_code(table, num_samples; code_frequency, sampling_frequency, phase=0, backend=…)

Iterator over `num_samples ÷ W` chunks (W = SIMD width), each yielding a `Vec{W,Int8}`
of resampled code chips. The chip index advances by an exact `step_numerator /
step_denominator` chips per sample via a drift-free DDA in the iteration state — no code
array is allocated. Same oversampling / denominator requirements as [`generate_code!`](@ref).
Any `num_samples % W` tail is not produced (handle it yourself if needed).
"""
function generate_code(table::CodeTable, step_numerator::Integer, step_denominator::Integer,
                       num_samples::Integer; phase::Integer = 0,
                       backend::Backend = default_backend(table))
    (0 < step_denominator ≤ _STEP_DEN) ||
        throw(ArgumentError("need 0 < step_denominator ≤ 2^$_B"))
    (0 < step_numerator ≤ step_denominator) ||
        throw(ArgumentError("need 0 < step_numerator ≤ step_denominator (must oversample)"))
    _make_code(table, Int(step_numerator), Int(step_denominator), Int(num_samples),
               Int(phase), backend, _vwidth(backend))
end
function generate_code(table::CodeTable, cps::Real, num_samples::Integer; kw...)
    sn, sd = _fixed_point_step(cps)
    generate_code(table, sn, sd, num_samples; kw...)
end
function generate_code(table::CodeTable, num_samples::Integer;
                       code_frequency::Real, sampling_frequency::Real, kw...)
    generate_code(table, chips_per_sample(code_frequency, sampling_frequency), num_samples; kw...)
end

function _make_code(table::CodeTable, step_num, step_den, num_samples, phase_offset,
                    backend, ::Val{W}) where {W}
    L = table.length; T = _phase_type(L)
    phase_init, rem_init = _init_state(Val(W), step_num, step_den, L, 0, phase_offset, T)
    prepared = prepare_code(table; backend = backend)
    CodeIterator{W,T,typeof(prepared)}(prepared, phase_init, rem_init,
        T(div(W * step_num, step_den) % L), _RemT(mod(W * step_num, step_den)),
        _RemT(step_den), T(L), num_samples ÷ W)
end

Base.length(it::CodeIterator) = it.num_chunks
Base.IteratorSize(::Type{<:CodeIterator}) = Base.HasLength()
Base.eltype(::Type{<:CodeIterator{W}}) where {W} = Vec{W,Int8}

@inline function Base.iterate(it::CodeIterator{W,T},
                              state = (it.phase_init, it.remainder_init, 0)) where {W,T}
    phase, remainder, chunk = state
    chunk >= it.num_chunks && return nothing
    result = it.prepared(phase)
    remainder += it.frac_step
    carry = remainder >= it.modulus
    remainder = vifelse(carry, remainder - it.modulus, remainder)
    phase = vifelse(carry, phase + it.whole_step + one(T), phase + it.whole_step)
    phase = vifelse(phase >= it.code_length, phase - it.code_length, phase)
    (result, (phase, remainder, chunk + 1))
end

# ---- AVX-512 iterator: relative-index DDA (see kernel.jl), yields Vec{64,Int8} ----
struct CodeIterator512
    padded::Vector{Int8}
    rel_init::Vec{64,Int8}
    rem_init::Vec{64,_RemT}
    base_init::Int
    whole_step::Int
    frac_step::_RemT
    modulus::_RemT
    code_length::Int
    num_chunks::Int
end

function _make_code(table::CodeTable, step_num, step_den, num_samples, phase_offset,
                    ::AVX512, ::Val{64})
    L = table.length
    rel_init, rem_init, base_init = _init_rel(Val(64), step_num, step_den, L, 0, phase_offset)
    CodeIterator512(table.padded, rel_init, rem_init, base_init,
        div(64 * step_num, step_den), _RemT(mod(64 * step_num, step_den)),
        _RemT(step_den), L, num_samples ÷ 64)
end

Base.length(it::CodeIterator512) = it.num_chunks
Base.IteratorSize(::Type{CodeIterator512}) = Base.HasLength()
Base.eltype(::Type{CodeIterator512}) = Vec{64,Int8}

@inline function Base.iterate(it::CodeIterator512,
                              state = (it.rel_init, it.rem_init, it.base_init, 0))
    rel, remainder, base, chunk = state
    chunk >= it.num_chunks && return nothing
    result = _rel_lookup(AVX512(), it.padded, rel, base)
    remainder += it.frac_step
    carry = remainder >= it.modulus
    remainder = vifelse(carry, remainder - it.modulus, remainder)
    h = Int(carry[1])
    rel = rel + vifelse(carry, one(Vec{64,Int8}), zero(Vec{64,Int8})) - Int8(h)
    base = mod(base + it.whole_step + h, it.code_length)
    (result, (rel, remainder, base, chunk + 1))
end

# ===== generate_code4: 4-way interleaved iterator =====
# Yields an `NTuple{4,Vec{W,Int8}}` per step (4·W samples), running four interleaved DDA
# states so their carry chains overlap — reaching full loop throughput even for trivial
# consumers (array fill, sum). Use plain `generate_code` when fusing into nontrivial work
# (that work supplies its own instruction-level parallelism). Mirrors
# `SinCosLUT.generate_carrier4`; destructure the 4-tuple in the loop header.

"""
    generate_code4(table, step_numerator, step_denominator, num_samples; phase=0, backend=…)
    generate_code4(table, chips_per_sample::Real, num_samples;            phase=0, backend=…)
    generate_code4(table, num_samples; code_frequency, sampling_frequency, phase=0, backend=…)

Like [`generate_code`](@ref) but yields four `Vec{W,Int8}` chunks per step (`4·W`
samples), running four interleaved DDA states so the carry chains overlap — full
throughput even for trivial consumers. **Destructure the 4-tuple in the loop header**
(`for (a,b,c,d) in generate_code4(...)`). Produces `num_samples ÷ (4W)` steps; handle any
tail yourself.
"""
function generate_code4(table::CodeTable, step_numerator::Integer, step_denominator::Integer,
                        num_samples::Integer; phase::Integer = 0,
                        backend::Backend = default_backend(table))
    (0 < step_denominator ≤ _STEP_DEN) ||
        throw(ArgumentError("need 0 < step_denominator ≤ 2^$_B"))
    (0 < step_numerator ≤ step_denominator) ||
        throw(ArgumentError("need 0 < step_numerator ≤ step_denominator (must oversample)"))
    _make_code4(table, Int(step_numerator), Int(step_denominator), Int(num_samples),
                Int(phase), backend, _vwidth(backend))
end
function generate_code4(table::CodeTable, cps::Real, num_samples::Integer; kw...)
    sn, sd = _fixed_point_step(cps)
    generate_code4(table, sn, sd, num_samples; kw...)
end
function generate_code4(table::CodeTable, num_samples::Integer;
                        code_frequency::Real, sampling_frequency::Real, kw...)
    generate_code4(table, chips_per_sample(code_frequency, sampling_frequency), num_samples; kw...)
end

# ---- AVX-512 4-way: relative-index DDA (== the generate_code! AVX-512 loop, yielded) ----
struct CodeIterator4_512
    padded::Vector{Int8}
    rel1::Vec{64,Int8}; rel2::Vec{64,Int8}; rel3::Vec{64,Int8}; rel4::Vec{64,Int8}
    rem1::Vec{64,_RemT}; rem2::Vec{64,_RemT}; rem3::Vec{64,_RemT}; rem4::Vec{64,_RemT}
    b1::Int; b2::Int; b3::Int; b4::Int
    whole_step::Int
    frac_step::_RemT
    modulus::_RemT
    code_length::Int
    num_steps::Int
end

function _make_code4(table::CodeTable, step_num, step_den, num_samples, phase_offset,
                     ::AVX512, ::Val{64})
    W = 64; L = table.length
    rel1, rem1, b1 = _init_rel(Val(W), step_num, step_den, L, 0,  phase_offset)
    rel2, rem2, b2 = _init_rel(Val(W), step_num, step_den, L, W,  phase_offset)
    rel3, rem3, b3 = _init_rel(Val(W), step_num, step_den, L, 2W, phase_offset)
    rel4, rem4, b4 = _init_rel(Val(W), step_num, step_den, L, 3W, phase_offset)
    CodeIterator4_512(table.padded, rel1, rel2, rel3, rel4, rem1, rem2, rem3, rem4,
        b1, b2, b3, b4, div(4W * step_num, step_den), _RemT(mod(4W * step_num, step_den)),
        _RemT(step_den), L, num_samples ÷ (4W))
end

Base.length(it::CodeIterator4_512) = it.num_steps
Base.IteratorSize(::Type{CodeIterator4_512}) = Base.HasLength()
Base.eltype(::Type{CodeIterator4_512}) = NTuple{4,Vec{64,Int8}}

@inline function Base.iterate(it::CodeIterator4_512,
                              state = (it.rel1, it.rel2, it.rel3, it.rel4,
                                       it.rem1, it.rem2, it.rem3, it.rem4,
                                       it.b1, it.b2, it.b3, it.b4, 0))
    rel1, rel2, rel3, rel4, rem1, rem2, rem3, rem4, b1, b2, b3, b4, step = state
    step >= it.num_steps && return nothing
    pd = it.padded
    result = (_rel_lookup(AVX512(), pd, rel1, b1), _rel_lookup(AVX512(), pd, rel2, b2),
              _rel_lookup(AVX512(), pd, rel3, b3), _rel_lookup(AVX512(), pd, rel4, b4))
    frac = it.frac_step; modulus = it.modulus; whole = it.whole_step; L = it.code_length
    z = zero(Vec{64,Int8}); o = one(Vec{64,Int8})
    a1 = rem1 + frac; c1 = a1 >= modulus; nrem1 = vifelse(c1, a1 - modulus, a1); h1 = Int(c1[1])
    a2 = rem2 + frac; c2 = a2 >= modulus; nrem2 = vifelse(c2, a2 - modulus, a2); h2 = Int(c2[1])
    a3 = rem3 + frac; c3 = a3 >= modulus; nrem3 = vifelse(c3, a3 - modulus, a3); h3 = Int(c3[1])
    a4 = rem4 + frac; c4 = a4 >= modulus; nrem4 = vifelse(c4, a4 - modulus, a4); h4 = Int(c4[1])
    nrel1 = rel1 + vifelse(c1, o, z) - Int8(h1); nrel2 = rel2 + vifelse(c2, o, z) - Int8(h2)
    nrel3 = rel3 + vifelse(c3, o, z) - Int8(h3); nrel4 = rel4 + vifelse(c4, o, z) - Int8(h4)
    nb1 = mod(b1 + whole + h1, L); nb2 = mod(b2 + whole + h2, L)
    nb3 = mod(b3 + whole + h3, L); nb4 = mod(b4 + whole + h4, L)
    (result, (nrel1, nrel2, nrel3, nrel4, nrem1, nrem2, nrem3, nrem4, nb1, nb2, nb3, nb4, step + 1))
end

# ---- AVX2 / Portable 4-way: phase-based DDA via the prepared callable ----
struct CodeIterator4{W,T,Prep}
    prepared::Prep
    phase1::Vec{W,T}; phase2::Vec{W,T}; phase3::Vec{W,T}; phase4::Vec{W,T}
    rem1::Vec{W,_RemT}; rem2::Vec{W,_RemT}; rem3::Vec{W,_RemT}; rem4::Vec{W,_RemT}
    whole_step::T
    frac_step::_RemT
    modulus::_RemT
    code_length::T
    num_steps::Int
end

function _make_code4(table::CodeTable, step_num, step_den, num_samples, phase_offset,
                     backend, ::Val{W}) where {W}
    L = table.length; T = _phase_type(L)
    p1, r1 = _init_state(Val(W), step_num, step_den, L, 0,  phase_offset, T)
    p2, r2 = _init_state(Val(W), step_num, step_den, L, W,  phase_offset, T)
    p3, r3 = _init_state(Val(W), step_num, step_den, L, 2W, phase_offset, T)
    p4, r4 = _init_state(Val(W), step_num, step_den, L, 3W, phase_offset, T)
    prepared = prepare_code(table; backend = backend)
    CodeIterator4{W,T,typeof(prepared)}(prepared, p1, p2, p3, p4, r1, r2, r3, r4,
        T(div(4W * step_num, step_den) % L), _RemT(mod(4W * step_num, step_den)),
        _RemT(step_den), T(L), num_samples ÷ (4W))
end

Base.length(it::CodeIterator4) = it.num_steps
Base.IteratorSize(::Type{<:CodeIterator4}) = Base.HasLength()
Base.eltype(::Type{<:CodeIterator4{W}}) where {W} = NTuple{4,Vec{W,Int8}}

@inline function Base.iterate(it::CodeIterator4{W,T},
                              state = (it.phase1, it.phase2, it.phase3, it.phase4,
                                       it.rem1, it.rem2, it.rem3, it.rem4, 0)) where {W,T}
    phase1, phase2, phase3, phase4, rem1, rem2, rem3, rem4, step = state
    step >= it.num_steps && return nothing
    p = it.prepared
    result = (p(phase1), p(phase2), p(phase3), p(phase4))
    frac = it.frac_step; modulus = it.modulus; whole = it.whole_step; Lc = it.code_length
    a1 = rem1 + frac; c1 = a1 >= modulus; nrem1 = vifelse(c1, a1 - modulus, a1)
    a2 = rem2 + frac; c2 = a2 >= modulus; nrem2 = vifelse(c2, a2 - modulus, a2)
    a3 = rem3 + frac; c3 = a3 >= modulus; nrem3 = vifelse(c3, a3 - modulus, a3)
    a4 = rem4 + frac; c4 = a4 >= modulus; nrem4 = vifelse(c4, a4 - modulus, a4)
    np1 = vifelse(c1, phase1 + whole + one(T), phase1 + whole)
    np2 = vifelse(c2, phase2 + whole + one(T), phase2 + whole)
    np3 = vifelse(c3, phase3 + whole + one(T), phase3 + whole)
    np4 = vifelse(c4, phase4 + whole + one(T), phase4 + whole)
    np1 = vifelse(np1 >= Lc, np1 - Lc, np1); np2 = vifelse(np2 >= Lc, np2 - Lc, np2)
    np3 = vifelse(np3 >= Lc, np3 - Lc, np3); np4 = vifelse(np4 >= Lc, np4 - Lc, np4)
    (result, (np1, np2, np3, np4, nrem1, nrem2, nrem3, nrem4, step + 1))
end
