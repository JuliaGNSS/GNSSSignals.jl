abstract type Modulation end
abstract type BOC <: Modulation end

struct BOCsin{M <: Union{AbstractFloat, Integer}, N <: Union{AbstractFloat, Integer}} <: BOC
    m::M
    n::N
    BOCsin(m, n) =
        m >= 1 && n >= 1 ?
        new{typeof(m), typeof(n)}(m, n) :
        error("m and n must be >= 1")
end

struct BOCcos{M <: Union{AbstractFloat, Integer}, N <: Union{AbstractFloat, Integer}} <: BOC
    m::M
    n::N
    BOCcos(m, n) =
        m >= 1 && n >= 1 ?
        new{typeof(m), typeof(n)}(m, n) :
        error("m and n must be >= 1")
end

struct LOC <: Modulation end

struct CBOC{B1 <: BOC, B2 <: BOC} <: BOC
    boc1::B1
    boc2::B2
    boc1_power::Float64
    CBOC(boc1, boc2, boc1_power) =
        0 < boc1_power < 1 && boc1.n == boc2.n ?
        new{typeof(boc1), typeof(boc2)}(boc1, boc2, boc1_power) :
        error("Power of BOC1 must be between 0 and 1 and n of both BOCs must match")
end

get_code_type(system::T) where T <: AbstractGNSS = get_code_type(get_modulation(T))

get_code_type(system::LOC) = Int16
get_code_type(system::BOC) = Int16
get_code_type(system::CBOC) = Float32

function get_subcarrier_code(modulation::BOCsin, phase)
    floored_subcarrier_phase = floor(Int, phase * 2 * modulation.m)
    iseven(floored_subcarrier_phase) << 1 - 1
end

function get_subcarrier_code(modulation::BOCcos, phase)
    floored_subcarrier_phase = floor(Int, (phase + 0.25) * 2 * modulation.m)
    iseven(floored_subcarrier_phase) << 1 - 1
end

function get_subcarrier_code(modulation::CBOC, phase)
    get_subcarrier_code(modulation.boc1, phase) * sqrt(modulation.boc1_power) +
        get_subcarrier_code(modulation.boc2, phase) * sqrt(1 - modulation.boc1_power)
end

get_floored_phase(modulation::BOC, phase) = floor(Int, phase * modulation.n)
get_floored_phase(modulation::CBOC, phase) = floor(Int, phase * modulation.boc1.n)

"""
$(SIGNATURES)

Get code of type <: `AbstractGNSS` at phase `phase` of PRN `prn`.
```julia-repl
julia> get_code(GPSL1(), 1200.3, 1)
```
"""
function get_code(
    system::T,
    phase,
    prn::Integer
) where T <: AbstractGNSS
    get_code(get_modulation(T), system, phase, prn)
end

"""
$(SIGNATURES)

Get code of BOC at
phase `phase` of PRN `prn`.
"""
function get_code(
    modulation::BOC,
    system::AbstractGNSS,
    phase,
    prn::Integer
)
    floored_phase = get_floored_phase(modulation, phase)
    modded_floored_phase = mod(
        floored_phase,
        get_code_length(system) * get_secondary_code_length(system)
    )
    get_code_at_index(system, modded_floored_phase, prn) *
        get_subcarrier_code(modulation, phase)
end

"""
$(SIGNATURES)

Get code of LOC at
phase `phase` of PRN `prn`.
"""
function get_code(
    modulation::LOC,
    system::AbstractGNSS,
    phase,
    prn::Integer
)
    floored_phase = floor(Int, phase)
    modded_floored_phase = mod(
        floored_phase,
        get_code_length(system) * get_secondary_code_length(system)
    )
    get_code_at_index(system, modded_floored_phase, prn)
end

Base.@propagate_inbounds function get_code_at_index(
    gnss::AbstractGNSS,
    phase::Integer,
    prn::Integer
)
    gnss.codes[phase + 1, prn]
end

"""
$(SIGNATURES)

Get code.
It is unsafe because it omits the modding.
The phase will not be wrapped by the code length. The phase has to be smaller
than the code length incl. secondary code.
```julia-repl
julia> get_code(GPSL1(), 1200.3, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(
    system::T,
    phase,
    prn::Integer
) where T <: AbstractGNSS
    get_code_unsafe(get_modulation(T), system, phase, prn)
end

"""
$(SIGNATURES)

Get code of BOC at
phase `phase` of PRN `prn`.
The phase will not be wrapped by the code length. The phase has to be smaller
than the code length incl. secondary code.
"""
Base.@propagate_inbounds function get_code_unsafe(
    modulation::BOC,
    system::AbstractGNSS,
    phase,
    prn::Integer
)
    floored_phase = get_floored_phase(modulation, phase)
    get_code_at_index(system, floored_phase, prn) *
        get_subcarrier_code(modulation, phase)
end

"""
$(SIGNATURES)

Get code of LOC at
phase `phase` of PRN `prn`.
The phase will not be wrapped by the code length. The phase has to be smaller
than the code length incl. secondary code.
"""
Base.@propagate_inbounds function get_code_unsafe(
    modulation::LOC,
    system::AbstractGNSS,
    phase,
    prn::Integer
)
    floored_phase = floor(Int, phase)
    get_code_at_index(system, floored_phase, prn)
end