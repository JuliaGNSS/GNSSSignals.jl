"""
Abstract supertype for all modulation types.
"""
abstract type Modulation end

"""
Abstract supertype for Binary Offset Carrier (BOC) modulation types.
"""
abstract type BOC <: Modulation end

"""
    BOCsin(m, n)

Sine-phased Binary Offset Carrier modulation.

BOC(m,n) uses a subcarrier frequency of `m * 1.023 MHz` and a code rate of `n * 1.023 Mcps`.

# Arguments
- `m`: Subcarrier frequency multiplier (must be ≥ 1)
- `n`: Code rate multiplier (must be ≥ 1)

# Example
```julia
boc11 = BOCsin(1, 1)  # BOC(1,1)
```
"""
struct BOCsin{M<:Union{AbstractFloat,Integer},N<:Union{AbstractFloat,Integer}} <: BOC
    m::M
    n::N
    BOCsin(m, n) =
        m >= 1 && n >= 1 ? new{typeof(m),typeof(n)}(m, n) : error("m and n must be >= 1")
end

"""
    BOCcos(m, n)

Cosine-phased Binary Offset Carrier modulation.

BOC(m,n) uses a subcarrier frequency of `m * 1.023 MHz` and a code rate of `n * 1.023 Mcps`.

# Arguments
- `m`: Subcarrier frequency multiplier (must be ≥ 1)
- `n`: Code rate multiplier (must be ≥ 1)

# Example
```julia
boc11 = BOCcos(1, 1)  # BOC(1,1) with cosine phase
```
"""
struct BOCcos{M<:Union{AbstractFloat,Integer},N<:Union{AbstractFloat,Integer}} <: BOC
    m::M
    n::N
    BOCcos(m, n) =
        m >= 1 && n >= 1 ? new{typeof(m),typeof(n)}(m, n) : error("m and n must be >= 1")
end

"""
    LOC()

Legacy/BPSK modulation (no subcarrier).

Used for GPS L1 C/A and GPS L5 signals.
"""
struct LOC <: Modulation end

"""
    CBOC(boc1, boc2, boc1_power)

Composite Binary Offset Carrier modulation.

CBOC combines two BOC modulations with specified power distribution.
Used for Galileo E1B signals as CBOC(6,1,1/11).

# Arguments
- `boc1`: First BOC component
- `boc2`: Second BOC component
- `boc1_power`: Power fraction allocated to first BOC (0 < power < 1)

# Example
```julia
cboc = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10/11)  # CBOC(6,1,1/11)
```
"""
struct CBOC{B1<:BOC,B2<:BOC} <: BOC
    boc1::B1
    boc2::B2
    boc1_power::Float32
    CBOC(boc1, boc2, boc1_power) =
        0 < boc1_power < 1 && boc1.n == boc2.n ?
        new{typeof(boc1),typeof(boc2)}(boc1, boc2, boc1_power) :
        error("Power of BOC1 must be between 0 and 1 and n of both BOCs must match")
end

"""
$(SIGNATURES)

Get the element type for code values of a GNSS system.

Returns the numeric type used to represent code values. For BPSK (LOC) signals,
this is typically `Int16`. For CBOC signals, this is a floating-point type.

# Arguments
- `system`: A GNSS system instance

# Returns
- `Type`: The element type for code values

# Examples
```julia-repl
julia> get_code_type(GPSL1())
Int16
julia> get_code_type(GalileoE1B())
Float32
```
"""
get_code_type(system::T) where {T<:AbstractGNSS} = get_code_type(system, get_modulation(T))

get_code_type(system::AbstractGNSS{<:AbstractMatrix{T}}, modulation) where {T} = T
get_code_type(system::AbstractGNSS{<:AbstractMatrix{T}}, modulation::CBOC) where {T} =
    promote_type(T, typeof(modulation.boc1_power))

get_code_factor(system::T) where {T<:AbstractGNSS} = get_code_factor(get_modulation(T))
get_code_factor(modulation::LOC) = 1
get_code_factor(modulation::BOC) = modulation.n
get_code_factor(modulation::CBOC) = modulation.boc1.n

"""
$(SIGNATURES)

Get the modulation type for a GNSS system.

# Arguments
- `system`: A GNSS system instance or type

# Returns
- `Modulation`: The modulation type (`LOC`, `BOCsin`, `BOCcos`, or `CBOC`)

# Examples
```julia-repl
julia> get_modulation(GPSL1())
LOC()
julia> get_modulation(GalileoE1B())
CBOC{BOCsin{Int64, Int64}, BOCsin{Int64, Int64}}(BOCsin{Int64, Int64}(1, 1), BOCsin{Int64, Int64}(6, 1), 0.90909094f0)
```
"""
get_modulation(s::T) where {T<:AbstractGNSS} = get_modulation(T)

"""
$(SIGNATURES)

Get the spectral power density of a GNSS signal at a given frequency.

Computes the power spectral density based on the system's modulation type.

# Arguments
- `system`: A GNSS system instance
- `f`: Baseband frequency at which to evaluate the spectrum

# Returns
- Spectral power density value

# Examples
```julia-repl
julia> using Unitful: kHz
julia> get_code_spectrum(GPSL1(), 0kHz)
9.775171065493646e-7
```
"""
get_code_spectrum(system, f) = get_code_spectrum(get_modulation(system), system, f)
get_code_spectrum(modulation::LOC, system, f) =
    get_code_spectrum_BPSK(get_code_frequency(system), f)
get_code_spectrum(modulation::BOCsin, system, f) = get_code_spectrum_BOCsin(
    modulation.n * get_code_frequency(system),
    modulation.m * get_code_frequency(system),
    f,
)
get_code_spectrum(modulation::BOCcos, system, f) = get_code_spectrum_BOCcos(
    modulation.n * get_code_frequency(system),
    modulation.m * get_code_frequency(system),
    f,
)
function get_code_spectrum(modulation::CBOC, system, f)
    get_code_spectrum(modulation.boc1, system, f) * modulation.boc1_power +
    get_code_spectrum(modulation.boc2, system, f) * (1 - modulation.boc1_power)
end

function get_subcarrier_code(modulation::BOCsin, phase::T) where {T<:Real}
    floored_subcarrier_phase = floor(Int, phase * 2 * modulation.m)
    iseven(floored_subcarrier_phase) * 2 - 1
end

function get_subcarrier_code(modulation::BOCcos, phase::T) where {T<:Real}
    get_subcarrier_code(BOCsin(modulation.m, modulation.n), phase + T(0.25))
end

# The amplitude is the sqrt of the power see
# https://galileognss.eu/wp-content/uploads/2015/12/Galileo_OS_SIS_ICD_v1.2.pdf
# Chapter 2.3.3. E1 Signal
function get_subcarrier_code(modulation::CBOC, phase::T) where {T<:Real}
    get_subcarrier_code(modulation.boc1, phase) * sqrt(modulation.boc1_power) +
    get_subcarrier_code(modulation.boc2, phase) * sqrt(1 - modulation.boc1_power)
end

get_floored_phase(modulation::LOC, phase) = floor(Int, phase)
get_floored_phase(modulation::BOC, phase) = floor(Int, phase * modulation.n)
get_floored_phase(modulation::CBOC, phase) = floor(Int, phase * modulation.boc1.n)

"""
$(SIGNATURES)

Get the code value at a given phase for a specific PRN.

Returns the spreading code value (including subcarrier modulation for BOC signals)
at the specified code phase. The phase is automatically wrapped to the code length.

# Arguments
- `system`: A GNSS system instance (e.g., `GPSL1()`, `GPSL5()`, `GalileoE1B()`)
- `phase`: Code phase in chips
- `prn`: PRN number of the satellite

# Returns
- Code value (typically `Int8` for BPSK, `Float32` for CBOC)

# Examples
```julia-repl
julia> get_code(GPSL1(), 0.0, 1)
1
julia> get_code(GPSL1(), 1200.3, 1)
-1
julia> get_code.(GPSL1(), 0:1022, 1)  # Full code period
```
"""
function get_code(system::T, phase, prn::Integer) where {T<:AbstractGNSS}
    get_code(get_modulation(T), system, phase, prn)
end

"""
$(SIGNATURES)

Get code value for BOC-modulated signals at a given phase.

Internal method that handles BOC modulation (sine or cosine phased).
"""
function get_code(modulation::BOC, system::AbstractGNSS, phase, prn::Integer)
    floored_phase = get_floored_phase(modulation, phase)
    modded_floored_phase = mod(floored_phase, size(system.codes, 1))
    get_code_at_index(system, modded_floored_phase, prn) *
    get_subcarrier_code(modulation, phase)
end

"""
$(SIGNATURES)

Get code value for LOC (BPSK) signals at a given phase.

Internal method that handles legacy/BPSK modulation without subcarrier.
"""
function get_code(modulation::LOC, system::AbstractGNSS, phase, prn::Integer)
    floored_phase = get_floored_phase(modulation, phase)
    modded_floored_phase = mod(floored_phase, size(system.codes, 1))
    get_code_at_index(system, modded_floored_phase, prn)
end

Base.@propagate_inbounds function get_code_at_index(
    gnss::AbstractGNSS,
    phase::Integer,
    prn::Integer,
)
    gnss.codes[phase+1, prn]
end

"""
$(SIGNATURES)

Get the code value at a given phase without bounds checking.

This is a faster version of [`get_code`](@ref) that skips the modulo operation.
The phase must be within `[0, code_length * secondary_code_length)`.

!!! warning
    Using phases outside the valid range results in undefined behavior.

# Arguments
- `system`: A GNSS system instance
- `phase`: Code phase in chips (must be within valid range)
- `prn`: PRN number of the satellite

# Returns
- Code value at the given phase

# Examples
```julia-repl
julia> get_code_unsafe(GPSL1(), 500.0, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(
    system::S,
    phase,
    prn::Integer,
) where {S<:AbstractGNSS}
    get_code_unsafe(get_modulation(S), system, phase, prn)
end

"""
$(SIGNATURES)

Get code value for BOC signals without bounds checking.

Internal method for BOC modulation without phase wrapping.
"""
Base.@propagate_inbounds function get_code_unsafe(
    modulation::BOC,
    system::AbstractGNSS,
    phase,
    prn::Integer,
)
    floored_phase = get_floored_phase(modulation, phase)
    get_code_at_index(system, floored_phase, prn) * get_subcarrier_code(modulation, phase)
end

"""
$(SIGNATURES)

Get code value for LOC (BPSK) signals without bounds checking.

Internal method for BPSK modulation without phase wrapping.
"""
Base.@propagate_inbounds function get_code_unsafe(
    modulation::LOC,
    system::AbstractGNSS,
    phase,
    prn::Integer,
)
    floored_phase = get_floored_phase(modulation, phase)
    get_code_at_index(system, floored_phase, prn)
end
