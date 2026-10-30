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

BOC(m,n) uses a subcarrier frequency of `m * 1.023 MHz` and a chip rate of `n * 1.023 Mcps`,
so a code chip spans `2m/n` subcarrier half-periods. A code phase counts chips, as everywhere
in the package. The code lookup table behind [`gen_code!`](@ref) needs `2m/n` to be a whole
number, and a code period to span whole subcarrier periods.

# Arguments
- `m`: Subcarrier frequency multiplier (must be ≥ 1)
- `n`: Chip rate multiplier (must be ≥ 1)

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

BOC(m,n) uses a subcarrier frequency of `m * 1.023 MHz` and a chip rate of `n * 1.023 Mcps`.
The subcarrier is the sine-phased one ([`BOCsin`](@ref)) shifted by a quarter of its period,
`n/(4m)` chips. The code lookup table behind [`gen_code!`](@ref) bakes it at `4m/n` quarter
periods per chip, which must be a whole number.

# Arguments
- `m`: Subcarrier frequency multiplier (must be ≥ 1)
- `n`: Chip rate multiplier (must be ≥ 1)

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
    BPSK()

Binary Phase Shift Keying modulation, with no subcarrier.

Used for GPS L1 C/A and GPS L5-I.
"""
struct BPSK <: Modulation end

"""
    CBOC(boc1, boc2, boc1_power, boc2_sign = 1)

Composite Binary Offset Carrier modulation.

CBOC combines two BOC modulations with specified power distribution.
Used for the Galileo E1 OS signals as CBOC(6,1,1/11).

The two CBOC components are added as
`sqrt(boc1_power) * boc1 + boc2_sign * sqrt(1 - boc1_power) * boc2`. The
`boc2_sign` selects the relative phase of the BOC(6,1) component: the
Galileo OS SIS ICD (§2.3.3) puts E1-B in phase (`+1`, "CBOC(+)") and E1-C
in anti-phase (`-1`, "CBOC(−)"), so the BOC(6,1) contributions of the two
components cancel in the composite. The sign does not affect the power
spectral density (see [`get_code_spectrum`](@ref)).

# Arguments
- `boc1`: First BOC component
- `boc2`: Second BOC component
- `boc1_power`: Power fraction allocated to first BOC (0 < power < 1)
- `boc2_sign`: Relative sign of the second BOC component (`+1` or `-1`, default `+1`)

# Example
```julia
cboc_b = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10/11)      # CBOC(6,1,1/11), E1-B
cboc_c = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10/11, -1)  # CBOC(6,1,1/11), E1-C
```
"""
struct CBOC{B1<:BOC,B2<:BOC} <: BOC
    boc1::B1
    boc2::B2
    boc1_power::Float32
    boc2_sign::Float32
    CBOC(boc1, boc2, boc1_power, boc2_sign = 1) =
        0 < boc1_power < 1 && boc1.n == boc2.n ?
        new{typeof(boc1),typeof(boc2)}(boc1, boc2, boc1_power, boc2_sign) :
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

Returns the type of the values [`get_code`](@ref) returns: the element type of the
code matrix (`Int8` for built-in signals), or a floating-point type for CBOC signals, whose two
subcarriers are weighted by their power split. [`gen_code!`](@ref) always fills
`Int8`.

# Arguments
- `signal`: A GNSS signal instance

# Returns
- `Type`: The element type for code values

# Examples
```julia-repl
julia> get_code_type(GPSL1CA())
Int8
julia> get_code_type(GalileoE1B_CBOC())
Float32
```
"""
get_code_type(signal::T) where {T<:AbstractGNSSSignal} = get_code_type(signal, get_modulation(T))

# `eltype(C)` rather than matching `AbstractMatrix{T}`: `C` is a direct parameter, so it
# always binds, whereas `T` would be unbound for a UnionAll code type like `GPSL1CA{Matrix}`
# (`Matrix <: AbstractMatrix` holds without fixing `T` — Aqua's unbound-parameter check).
get_code_type(signal::AbstractGNSSSignal{C}, modulation) where {C} = eltype(C)
get_code_type(signal::AbstractGNSSSignal{C}, modulation::CBOC) where {C} =
    promote_type(eltype(C), typeof(modulation.boc1_power))

"""
$(SIGNATURES)

Get the modulation type for a GNSS signal.

# Arguments
- `signal`: A GNSS signal instance or type

# Returns
- `Modulation`: The modulation type (`BPSK`, `BOCsin`, `BOCcos`, `CBOC` or `TMBOC`)

Each concrete signal states its modulation once, on its own type
(`get_modulation(::Type{<:GPSL1CA}) = BPSK()`, in that signal's file) — the form
the code LUT is baked from; this method forwards an instance to that type, so
the two can never disagree.

# Examples
```julia-repl
julia> get_modulation(GPSL1CA())
BPSK()
julia> get_modulation(GalileoE1B_CBOC())
CBOC{BOCsin{Int64, Int64}, BOCsin{Int64, Int64}}(BOCsin{Int64, Int64}(1, 1), BOCsin{Int64, Int64}(6, 1), 0.90909094f0)
```
"""
@inline get_modulation(signal::T) where {T<:AbstractGNSSSignal} = get_modulation(T)

# The phase is in code chips, as everywhere in the package. A chip of BOC(m, n) spans
# `2m/n` subcarrier half-periods: the subcarrier runs at m · 1.023 MHz, the code at n · 1.023
# MHz. The phase is multiplied by that factor rather than divided by `n`: the factor does not
# depend on the phase, so a loop computes it once, and in the phase's own float type it gives
# for n = 1 exactly the `phase · 2m` this was before n was carried. An `Int8` sign, so it keeps
# the code matrix's element type (`get_code_type`) rather than widening the chip to `Int`.
function get_subcarrier_code(modulation::BOCsin, phase::Real)
    half_periods_per_chip = oftype(float(phase), 2 * modulation.m / modulation.n)
    floored_subcarrier_phase = floor(Int, phase * half_periods_per_chip)
    iseven(floored_subcarrier_phase) ? Int8(1) : Int8(-1)
end

# The cosine subcarrier is the sine one a quarter subcarrier period later: sign(cos x) =
# sign(sin(x + π/2)). A subcarrier period is `n/m` chips, so the quarter is `n/(4m)` chips.
function get_subcarrier_code(modulation::BOCcos, phase::Real)
    quarter_period = oftype(float(phase), modulation.n / (4 * modulation.m))
    get_subcarrier_code(BOCsin(modulation.m, modulation.n), phase + quarter_period)
end

# The amplitude is the sqrt of the power see
# https://galileosignal.eu/wp-content/uploads/2015/12/Galileo_OS_SIS_ICD_v1.2.pdf
# Chapter 2.3.3. E1 Signal
function get_subcarrier_code(modulation::CBOC, phase::T) where {T<:Real}
    get_subcarrier_code(modulation.boc1, phase) * sqrt(modulation.boc1_power) +
    modulation.boc2_sign *
    get_subcarrier_code(modulation.boc2, phase) * sqrt(1 - modulation.boc1_power)
end

# TMBOC's subcarrier value at a given phase: pick the BOC variant for
# the current primary-chip position (mod pattern length), then evaluate
# that BOC's subcarrier at the same phase.
function get_subcarrier_code(modulation::TMBOC{B1, B2, N}, phase::T) where {B1, B2, N, T<:Real}
    chip_pos = mod(floor(Int, phase), N)
    @inbounds use_boc2 = modulation.pattern[chip_pos + 1]
    use_boc2 ? get_subcarrier_code(modulation.boc2, phase) :
               get_subcarrier_code(modulation.boc1, phase)
end

# The chip a phase falls in: the phase counts chips for every modulation.
get_floored_phase(::Modulation, phase) = floor(Int, phase)

"""
$(SIGNATURES)

Get the code value at a given phase for a specific PRN.

Returns the spreading code value (including subcarrier modulation for BOC signals)
at the specified code phase. The phase is automatically wrapped to the code length.

# Arguments
- `signal`: A GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B_CBOC()`)
- `phase`: Code phase in chips
- `prn`: PRN number of the satellite

# Returns
- Code value, of type [`get_code_type`](@ref)`(signal)`: `Int8`, or `Float32` for CBOC

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

# Multiply the primary chip by its secondary-code chip. The absolute chip is
# either wrapped by get_code or supplied in range by get_code_unsafe.
Base.@propagate_inbounds function _primary_secondary_code(signal, absolute_chip, prn)
    primary_length = size(signal.codes, 1)
    sec = get_secondary_code(signal)
    chip_idx = mod(absolute_chip, primary_length)
    sec_idx = div(absolute_chip, primary_length)
    get_code_at_index(signal, chip_idx, prn) * secondary_value(sec, prn, sec_idx)
end

@inline _modulate_code(::BPSK, code, phase) = code
@inline _modulate_code(modulation::BOC, code, phase) =
    code * get_subcarrier_code(modulation, phase)

function get_code(modulation::Union{BOC,BPSK}, signal::AbstractGNSSSignal, phase, prn::Integer)
    floored_phase = get_floored_phase(modulation, phase)
    primary_length = size(signal.codes, 1)
    sec_len = secondary_code_length(get_secondary_code(signal))
    absolute_chip = mod(floored_phase, primary_length * sec_len)
    _modulate_code(modulation, _primary_secondary_code(signal, absolute_chip, prn), phase)
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

# The caller supplies a phase within the primary × secondary code length.
Base.@propagate_inbounds function get_code_unsafe(
    modulation::Union{BOC,BPSK},
    signal::AbstractGNSSSignal,
    phase,
    prn::Integer,
)
    floored_phase = get_floored_phase(modulation, phase)
    _modulate_code(modulation, _primary_secondary_code(signal, floored_phase, prn), phase)
end
