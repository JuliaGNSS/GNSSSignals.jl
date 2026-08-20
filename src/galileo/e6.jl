"""
    GalileoE6B{C} <: AbstractGalileoSignal{C}

Galileo E6-B signal (the data-carrying component of Galileo E6).

BPSK(5) modulation ([`LOC`](@ref)) of a 5115-chip memory code at 5.115 Mcps on
the E6 band at 1278.75 MHz, so [`get_band`](@ref) returns [`E6`](@ref). The
1 ms primary code carries no secondary code
([`get_secondary_code`](@ref) is [`NoSecondaryCode`](@ref)), and the component
carries the C/NAV message at 1000 symbols/s — the channel that broadcasts the
Galileo High Accuracy Service (HAS) corrections.

The primary codes are optimised pseudo-random memory codes, not register
sequences, so they are stored rather than generated (Galileo E6-B/C Codes
Technical Note, Issue 1, §2.3.1 and Annex A.3). Codes 1-50 are defined; SVID
`n` uses code `n`.

# Example
```julia
e6b = GalileoE6B()
get_code_length(e6b)      # 5115
get_band(e6b)             # E6()
get_data_frequency(e6b)   # 1000 Hz (C/NAV)
```
"""
struct GalileoE6B{C<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GalileoE6C{C, M} <: AbstractGalileoSignal{C}

Galileo E6-C signal (the dataless pilot component of Galileo E6).

BPSK(5) modulation ([`LOC`](@ref)) of a 5115-chip memory code at 5.115 Mcps on
the E6 band (1278.75 MHz, reported as [`E6`](@ref)), overlaid with the 100-chip
per-SVID CS100 secondary code for a 100 ms tiered code. As the pilot component
it carries no navigation data, so [`get_data_frequency`](@ref) returns 0 Hz.

The primary codes are memory codes (Galileo E6-B/C Codes Technical Note, Issue
1, §2.3.1 and Annex A.4); the secondary codes are the same CS100_1-50 the
[`GalileoE5aQ`](@ref) pilot uses, assigned CS100ₙ to SVID `n` (Technical Note
§2.4). Codes 1-50 are defined.

!!! note
    The E6-C ranging codes are currently transmitted unencrypted, but the
    Technical Note states encryption is planned for future use, at which point
    the codes here stop matching the signal in space. E6-C availability is not
    guaranteed the way an Open Service component's is.

# Example
```julia
e6c = GalileoE6C()
get_code_length(e6c)            # 5115
get_secondary_code_length(e6c)  # 100
get_data_frequency(e6c)         # 0 Hz
```
"""
struct GalileoE6C{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    secondary_codes::M    # 100 × 50 Int8 ±1 CS100 matrix, exposed via PerPRNSecondaryCode
    lut::SignalLUT        # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

const E6_CODE_LENGTH = 5115
const E6_NUM_PRNS = 50

#= The E6-B/C primary codes are memory codes, so they ship as `Int8` ±1 blobs
(`data/codes_galileo_e6b.bin`, `data/codes_galileo_e6c.bin`) exactly as the E1
codes do. Both were decoded from the hexadecimal listings embedded in the
Galileo E6-B/C Codes Technical Note (Annex A, files `A3_E6B.TXT` /
`A4_E6C.TXT`), MSB-first, four chips per hex symbol, dropping the single
zero the ICD pads the last symbol with (5115 is not divisible by 4). =#
read_galileo_e6b_codes() = read_in_codes(
    Int8,
    joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_galileo_e6b.bin"),
    E6_NUM_PRNS,
    E6_CODE_LENGTH,
)

read_galileo_e6c_codes() = read_in_codes(
    Int8,
    joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_galileo_e6c.bin"),
    E6_NUM_PRNS,
    E6_CODE_LENGTH,
)

function GalileoE6B()
    codes = widen_codes_to_storage(read_galileo_e6b_codes())
    lut = build_signal_lut(get_modulation(GalileoE6B), codes, NoSecondaryCode())
    GalileoE6B(codes, lut)
end

function GalileoE6C()
    codes = widen_codes_to_storage(read_galileo_e6c_codes())
    secondary = _build_galileo_cs100_secondary(1:50)
    # The 100-chip per-SVID CS100 overlay is too long to bake (100·5115·1 > typemax(Int16)),
    # so it stays residual in the SignalLUT and is applied per primary period at gen time.
    lut = build_signal_lut(get_modulation(GalileoE6C), codes, PerPRNSecondaryCode(secondary))
    GalileoE6C(codes, secondary, lut)
end

# Shared interface (band, modulation, frequencies).

get_modulation(::Type{<:GalileoE6B}) = LOC()
get_modulation(::Type{<:GalileoE6C}) = LOC()

@inline get_band(::Type{<:GalileoE6B}) = E6()
@inline get_band(::Type{<:GalileoE6C}) = E6()

# E6-B and E6-C share one carrier component, and the OS SIS ICD v2.2 Eq. 10
# subtracts the pilot: s_E6 = (e_E6-B − e_E6-C)/√2. That 180° is the pilot's
# phase against the E6-B reference, and — unlike E1, where the anti-phase sits in
# the CBOC(−) code — nothing in the plain-BPSK E6-C definition carries it, so it
# lives here. See [`get_carrier_phase_offset`](@ref).
@inline get_carrier_phase_offset(::Type{<:GalileoE6C}) = Float64(π)

# 50/50 E6-B/E6-C power sharing (Galileo OS SIS ICD v2.2 §2.3.2), of the E6
# composite whose −155.25 dBW total (Table 13) is the unit here. See
# [`get_relative_power`](@ref).
@inline get_relative_power(::Type{<:GalileoE6B}) = 0.5
@inline get_relative_power(::Type{<:GalileoE6C}) = 0.5

@inline get_signal_name(::Type{<:GalileoE6B}) = "Galileo E6-B"
@inline get_signal_name(::Type{<:GalileoE6C}) = "Galileo E6-C"

"""
$(SIGNATURES)

Get the code length for Galileo E6 (5115 chips, both components).
"""
@inline get_code_length(::Type{<:GalileoE6B}) = E6_CODE_LENGTH
@inline get_code_length(::Type{<:GalileoE6C}) = E6_CODE_LENGTH

"""
$(SIGNATURES)

Get the code chipping rate for Galileo E6 (5.115 MHz, both components — Galileo
OS SIS ICD v2.2, Table 9).
"""
@inline get_code_frequency(::Type{<:GalileoE6B}) = 5_115_000Hz
@inline get_code_frequency(::Type{<:GalileoE6C}) = 5_115_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for Galileo E6-B.

The E6-B channel carries the C/NAV message at 1000 symbols/s (Galileo OS SIS
ICD v2.2, Table 9), the rate the HAS SIS ICD's one page per second builds on.

# Returns
- `Frequency`: 1000 Hz
"""
@inline get_data_frequency(::Type{<:GalileoE6B}) = 1000Hz

"""
$(SIGNATURES)

Get the data symbol rate for Galileo E6-C.

The E6-C component is a dataless pilot.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_data_frequency(::Type{<:GalileoE6C}) = 0Hz

"""
$(SIGNATURES)

Get the secondary code for Galileo E6-B.

Galileo E6-B has no secondary code — its 1 ms primary code repeats untiered
(Galileo E6-B/C Codes Technical Note, Table 1 and Table 3); the overlay lives
on the E6-C pilot.

# Returns
- [`NoSecondaryCode`](@ref)
"""
@inline get_secondary_code(::GalileoE6B) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the secondary code for Galileo E6-C.

Each SVID overlays its primary code with its 100-chip CS100 sequence (Galileo
E6-B/C Codes Technical Note §2.4, the OS SIS ICD's CS100_1-50), giving a 100 ms
tiered code.

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 100 × 50 CS100 matrix
"""
@inline get_secondary_code(s::GalileoE6C) = PerPRNSecondaryCode(s.secondary_codes)
