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

Linear Offset Carrier modulation — the BPSK-like baseline with no
subcarrier, named in contrast to [`BOC`](@ref) (Binary Offset Carrier).

Used for GPS L1 C/A and GPS L5-I.
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
    TMBOC(boc1, boc2, pattern)

Time-Multiplexed Binary Offset Carrier modulation.

Used for the GPS L1C pilot (L1C-P) component. Within each repeating
block of `length(pattern)` primary chips, `pattern[k]` selects which
BOC variant to apply at chip-position `k - 1`: `true` for `boc2`,
`false` for `boc1`.

For GPS L1C-P specifically, the pattern is TMBOC(6,1,4/33): every
33 primary chips, the four positions `{0, 4, 6, 29}` use BOC(6,1) and
the remaining 29 use BOC(1,1). See IS-GPS-800G §3.3.

# Arguments
- `boc1`: BOC variant for the majority positions (`pattern[k] == false`)
- `boc2`: BOC variant for the minority positions (`pattern[k] == true`)
- `pattern`: `NTuple{N, Bool}` selecting BOC variant per chip-position
  within one repeat block. `N` must match the signal's TMBOC period.

The two BOC components must share the same code-rate multiplier `n`.

# Example
```julia
# GPS L1C-P: BOC(6,1) at positions 0, 4, 6, 29 within every 33 chips
pattern = ntuple(k -> (k - 1) ∈ (0, 4, 6, 29), 33)
tmboc = TMBOC(BOCsin(1, 1), BOCsin(6, 1), pattern)
```
"""
struct TMBOC{B1<:BOC,B2<:BOC,N} <: BOC
    boc1::B1
    boc2::B2
    pattern::NTuple{N, Bool}
    function TMBOC(boc1::B1, boc2::B2, pattern::NTuple{N, Bool}) where {B1<:BOC, B2<:BOC, N}
        boc1.n == boc2.n || error("n of both BOCs must match")
        N >= 1 || error("TMBOC pattern must be non-empty")
        new{B1, B2, N}(boc1, boc2, pattern)
    end
end

"""
$(SIGNATURES)

Get the element type for code values of a GNSS signal.

Returns the numeric type used to represent code values. For BPSK (LOC) signals,
this is typically `Int16`. For CBOC signals, this is a floating-point type.

# Arguments
- `signal`: A GNSS signal instance

# Returns
- `Type`: The element type for code values

# Examples
```julia-repl
julia> get_code_type(GPSL1CA())
Int16
julia> get_code_type(GalileoE1B())
Float32
```
"""
get_code_type(signal::T) where {T<:AbstractGNSSSignal} = get_code_type(signal, get_modulation(T))

get_code_type(signal::AbstractGNSSSignal{<:AbstractMatrix{T}}, modulation) where {T} = T
get_code_type(signal::AbstractGNSSSignal{<:AbstractMatrix{T}}, modulation::CBOC) where {T} =
    promote_type(T, typeof(modulation.boc1_power))

get_code_factor(signal::T) where {T<:AbstractGNSSSignal} = get_code_factor(get_modulation(T))
get_code_factor(modulation::LOC) = 1
get_code_factor(modulation::BOC) = modulation.n
get_code_factor(modulation::CBOC) = modulation.boc1.n
get_code_factor(modulation::TMBOC) = modulation.boc1.n

"""
$(SIGNATURES)

Get the modulation type for a GNSS signal.

# Arguments
- `signal`: A GNSS signal instance or type

# Returns
- `Modulation`: The modulation type (`LOC`, `BOCsin`, `BOCcos`, or `CBOC`)

# Examples
```julia-repl
julia> get_modulation(GPSL1CA())
LOC()
julia> get_modulation(GalileoE1B())
CBOC{BOCsin{Int64, Int64}, BOCsin{Int64, Int64}}(BOCsin{Int64, Int64}(1, 1), BOCsin{Int64, Int64}(6, 1), 0.90909094f0)
```
"""
get_modulation(signal::T) where {T<:AbstractGNSSSignal} = get_modulation(T)

"""
$(SIGNATURES)

Get the spectral power density of a GNSS signal at a given frequency.

Computes the power spectral density based on the signal's modulation type.

# Arguments
- `signal`: A GNSS signal instance
- `f`: Baseband frequency at which to evaluate the spectrum

# Returns
- Spectral power density value

# Examples
```julia-repl
julia> using Unitful: kHz
julia> get_code_spectrum(GPSL1CA(), 0kHz)
9.775171065493646e-7
```
"""
get_code_spectrum(signal, f) = get_code_spectrum(get_modulation(signal), signal, f)
get_code_spectrum(modulation::LOC, signal, f) =
    get_code_spectrum_BPSK(get_code_frequency(signal), f)
get_code_spectrum(modulation::BOCsin, signal, f) = get_code_spectrum_BOCsin(
    modulation.n * get_code_frequency(signal),
    modulation.m * get_code_frequency(signal),
    f,
)
get_code_spectrum(modulation::BOCcos, signal, f) = get_code_spectrum_BOCcos(
    modulation.n * get_code_frequency(signal),
    modulation.m * get_code_frequency(signal),
    f,
)
function get_code_spectrum(modulation::CBOC, signal, f)
    get_code_spectrum(modulation.boc1, signal, f) * modulation.boc1_power +
    get_code_spectrum(modulation.boc2, signal, f) * (1 - modulation.boc1_power)
end

# For TMBOC, weight by the fraction of chips that use each BOC variant.
function get_code_spectrum(modulation::TMBOC{B1, B2, N}, signal, f) where {B1, B2, N}
    n_boc2 = count(modulation.pattern)
    boc2_frac = n_boc2 / N
    get_code_spectrum(modulation.boc1, signal, f) * (1 - boc2_frac) +
    get_code_spectrum(modulation.boc2, signal, f) * boc2_frac
end

function get_subcarrier_code(modulation::BOCsin, phase::T) where {T<:Real}
    floored_subcarrier_phase = floor(Int, phase * 2 * modulation.m)
    iseven(floored_subcarrier_phase) * 2 - 1
end

function get_subcarrier_code(modulation::BOCcos, phase::T) where {T<:Real}
    get_subcarrier_code(BOCsin(modulation.m, modulation.n), phase + T(0.25))
end

# The amplitude is the sqrt of the power see
# https://galileosignal.eu/wp-content/uploads/2015/12/Galileo_OS_SIS_ICD_v1.2.pdf
# Chapter 2.3.3. E1 Signal
function get_subcarrier_code(modulation::CBOC, phase::T) where {T<:Real}
    get_subcarrier_code(modulation.boc1, phase) * sqrt(modulation.boc1_power) +
    get_subcarrier_code(modulation.boc2, phase) * sqrt(1 - modulation.boc1_power)
end

# TMBOC's subcarrier value at a given phase: pick the BOC variant for
# the current primary-chip position (mod pattern length), then evaluate
# that BOC's subcarrier at the same phase.
function get_subcarrier_code(modulation::TMBOC{B1, B2, N}, phase::T) where {B1, B2, N, T<:Real}
    chip_pos = mod(floor(Int, phase * modulation.boc1.n), N)
    @inbounds use_boc2 = modulation.pattern[chip_pos + 1]
    use_boc2 ? get_subcarrier_code(modulation.boc2, phase) :
               get_subcarrier_code(modulation.boc1, phase)
end

get_floored_phase(modulation::LOC, phase) = floor(Int, phase)
get_floored_phase(modulation::BOC, phase) = floor(Int, phase * modulation.n)
get_floored_phase(modulation::CBOC, phase) = floor(Int, phase * modulation.boc1.n)
get_floored_phase(modulation::TMBOC, phase) = floor(Int, phase * modulation.boc1.n)

"""
$(SIGNATURES)

Get the code value at a given phase for a specific PRN.

Returns the spreading code value (including subcarrier modulation for BOC signals)
at the specified code phase. The phase is automatically wrapped to the code length.

# Arguments
- `signal`: A GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B()`)
- `phase`: Code phase in chips
- `prn`: PRN number of the satellite

# Returns
- Code value (typically `Int8` for BPSK, `Float32` for CBOC)

# Examples
```julia-repl
julia> get_code(GPSL1CA(), 0.0, 1)
1
julia> get_code(GPSL1CA(), 1200.3, 1)
-1
julia> get_code.(GPSL1CA(), 0:1022, 1)  # Full code period
```
"""
function get_code(signal::T, phase, prn::Integer) where {T<:AbstractGNSSSignal}
    get_code(get_modulation(T), signal, phase, prn)
end

"""
$(SIGNATURES)

Get code value for BOC-modulated signals at a given phase.

Internal method that handles BOC modulation (sine or cosine phased).
"""
function get_code(modulation::BOC, signal::AbstractGNSSSignal, phase, prn::Integer)
    floored_phase = get_floored_phase(modulation, phase)
    primary_length = size(signal.codes, 1)
    sec = get_secondary_code(signal)
    sec_len = secondary_code_length(sec)
    absolute_chip = mod(floored_phase, primary_length * sec_len)
    chip_idx = mod(absolute_chip, primary_length)
    sec_idx = div(absolute_chip, primary_length)
    get_code_at_index(signal, chip_idx, prn) *
    secondary_value(sec, prn, sec_idx) *
    get_subcarrier_code(modulation, phase)
end

"""
$(SIGNATURES)

Get code value for LOC (BPSK) signals at a given phase.

Internal method that handles legacy/BPSK modulation without subcarrier.
"""
function get_code(modulation::LOC, signal::AbstractGNSSSignal, phase, prn::Integer)
    floored_phase = get_floored_phase(modulation, phase)
    primary_length = size(signal.codes, 1)
    sec = get_secondary_code(signal)
    sec_len = secondary_code_length(sec)
    absolute_chip = mod(floored_phase, primary_length * sec_len)
    chip_idx = mod(absolute_chip, primary_length)
    sec_idx = div(absolute_chip, primary_length)
    get_code_at_index(signal, chip_idx, prn) * secondary_value(sec, prn, sec_idx)
end

Base.@propagate_inbounds function get_code_at_index(
    signal::AbstractGNSSSignal,
    phase::Integer,
    prn::Integer,
)
    signal.codes[phase+1, prn]
end

"""
$(SIGNATURES)

Get the code value at a given phase without bounds checking.

This is a faster version of [`get_code`](@ref) that skips the modulo operation.
The phase must be within `[0, code_length * secondary_code_length)`.

!!! warning
    Using phases outside the valid range results in undefined behavior.

# Arguments
- `signal`: A GNSS signal instance
- `phase`: Code phase in chips (must be within valid range)
- `prn`: PRN number of the satellite

# Returns
- Code value at the given phase

# Examples
```julia-repl
julia> get_code_unsafe(GPSL1CA(), 500.0, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(
    signal::S,
    phase,
    prn::Integer,
) where {S<:AbstractGNSSSignal}
    get_code_unsafe(get_modulation(S), signal, phase, prn)
end

"""
$(SIGNATURES)

Get code value for BOC signals without bounds checking.

Internal method for BOC modulation without phase wrapping. Assumes
`phase` is within `[0, primary_length * secondary_length)`.
"""
Base.@propagate_inbounds function get_code_unsafe(
    modulation::BOC,
    signal::AbstractGNSSSignal,
    phase,
    prn::Integer,
)
    floored_phase = get_floored_phase(modulation, phase)
    primary_length = size(signal.codes, 1)
    sec = get_secondary_code(signal)
    chip_idx = mod(floored_phase, primary_length)
    sec_idx = div(floored_phase, primary_length)
    get_code_at_index(signal, chip_idx, prn) *
    secondary_value(sec, prn, sec_idx) *
    get_subcarrier_code(modulation, phase)
end

"""
$(SIGNATURES)

Get code value for LOC (BPSK) signals without bounds checking.

Internal method for BPSK modulation without phase wrapping. Assumes
`phase` is within `[0, primary_length * secondary_length)`.
"""
Base.@propagate_inbounds function get_code_unsafe(
    modulation::LOC,
    signal::AbstractGNSSSignal,
    phase,
    prn::Integer,
)
    floored_phase = get_floored_phase(modulation, phase)
    primary_length = size(signal.codes, 1)
    sec = get_secondary_code(signal)
    chip_idx = mod(floored_phase, primary_length)
    sec_idx = div(floored_phase, primary_length)
    get_code_at_index(signal, chip_idx, prn) * secondary_value(sec, prn, sec_idx)
end
