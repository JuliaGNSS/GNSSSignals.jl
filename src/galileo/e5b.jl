"""
    GalileoE5bI{C} <: AbstractGalileoSignal{C}

Galileo E5b-I signal (the in-phase, data-carrying component of Galileo E5b).

10230-chip primary code at 10.23 Mcps on the E5b sideband at 1207.14 MHz, so
[`get_band`](@ref) returns [`E5b`](@ref) — a carrier Galileo shares exactly with
BeiDou B2b. A 4-bit
secondary code (CS4, the same for every SVID) overlays the data channel, giving
a 4 ms tiered code, and the component carries the I/NAV message at 250
symbols/s.

Galileo E5b is one half of the wideband E5 AltBOC(15,10) signal. The ICD notes
that "E5a and E5b signals can be processed independently by the user receiver as
though they were two separate QPSK signals" (OS SIS ICD v2.2 §2.3.1.2), which is
what this implementation models: the E5b sideband on its own as BPSK(10)
(modulation [`LOC`](@ref)), like [`GalileoE5aI`](@ref) does for E5a.

The primary code is generated from two 14-stage shift registers per the Galileo
OS SIS ICD (§3.4.1): a common base register 1 (all-ones start, feedback
`64021₈`) XOR'd with a per-SVID base register 2 (feedback `51445₈`, start values
`E5B_I_X2_INIT`), truncated to 10230 chips. PRNs 1-50 are supported.

# Example
```julia
e5b_i = GalileoE5bI()
get_code_length(e5b_i)            # 10230
get_secondary_code_length(e5b_i)  # 4
get_band(e5b_i)                   # E5b()
```
"""
struct GalileoE5bI{C<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GalileoE5bQ{C, M} <: AbstractGalileoSignal{C}

Galileo E5b-Q signal (the quadrature, dataless pilot component of Galileo E5b).

10230-chip primary code at 10.23 Mcps on the E5b sideband (1207.14 MHz, reported
as [`E5b`](@ref)), overlaid with a 100-bit per-SVID secondary code for a 100 ms
tiered code. As the pilot component it carries no navigation data, so
[`get_data_frequency`](@ref) returns 0 Hz.

Like [`GalileoE5bI`](@ref) the E5b sideband is modelled on its own as BPSK(10)
([`LOC`](@ref)). The primary code uses the E5b base register 1 (feedback
`64021₈`) with the E5b-Q register 2 (feedback `43143₈`) and the per-SVID start
values `E5B_Q_X2_INIT`. The secondary codes are the CS100 codes 51-100 —
CS100₍ₙ₊₅₀₎ to SVID `n` (OS SIS ICD v2.2 §3.5.2) — the other half of the table
[`GalileoE5aQ`](@ref) draws from. PRNs 1-50 are supported.

# Example
```julia
e5b_q = GalileoE5bQ()
get_code_length(e5b_q)            # 10230
get_secondary_code_length(e5b_q)  # 100
get_data_frequency(e5b_q)         # 0 Hz
```
"""
struct GalileoE5bQ{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    secondary_codes::M    # 100 × 50 Int8 ±1 matrix, exposed via PerPRNSecondaryCode
    lut::SignalLUT        # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

#= E5b primary code generation (Galileo OS SIS ICD v2.2 §3.4.1, Table 16).

The generator is the one shared with E5a in `galileo/codes.jl`; E5b differs in
both feedback polynomials, and its two components differ from each other in the
register-2 polynomial as well as the start values. =#

const E5B_X1_TAP = 0o64021   # base register 1 feedback polynomial (both components)
const E5B_I_X2_TAP = 0o51445 # base register 2 feedback polynomial, E5b-I
const E5B_Q_X2_TAP = 0o43143 # base register 2 feedback polynomial, E5b-Q

# E5b-I base register 2 start values (octal), PRN 1-50 (ICD Table 19).
const E5B_I_X2_INIT = [
    0o07220, 0o26047, 0o00252, 0o17166, 0o14161, 0o02540, 0o01537, 0o26023, 0o01725, 0o20637,
    0o02364, 0o27731, 0o30640, 0o34174, 0o06464, 0o07676, 0o32231, 0o10353, 0o00755, 0o26077,
    0o11644, 0o11537, 0o35115, 0o20452, 0o34645, 0o25664, 0o21403, 0o32253, 0o02337, 0o30777,
    0o27122, 0o22377, 0o36175, 0o33075, 0o33151, 0o13134, 0o07433, 0o10216, 0o35466, 0o02533,
    0o05351, 0o30121, 0o14010, 0o32576, 0o30326, 0o37433, 0o26022, 0o35770, 0o06670, 0o12017,
]

# E5b-Q base register 2 start values (octal), PRN 1-50 (ICD Table 20).
const E5B_Q_X2_INIT = [
    0o03331, 0o06143, 0o25322, 0o23371, 0o00413, 0o36235, 0o17750, 0o04745, 0o13005, 0o37140,
    0o30155, 0o20237, 0o03461, 0o31662, 0o27146, 0o05547, 0o02456, 0o30013, 0o00322, 0o10761,
    0o26767, 0o36004, 0o30713, 0o07662, 0o21610, 0o20134, 0o11262, 0o10706, 0o34143, 0o11051,
    0o25460, 0o17665, 0o32354, 0o21230, 0o20146, 0o11362, 0o37246, 0o16344, 0o15034, 0o25471,
    0o25646, 0o22157, 0o04336, 0o16356, 0o04075, 0o02626, 0o11706, 0o37011, 0o27041, 0o31024,
]

function GalileoE5bI()
    codes = widen_codes_to_storage(read_galileo_e5_codes(E5B_X1_TAP, E5B_I_X2_TAP, E5B_I_X2_INIT))
    lut = build_signal_lut(get_modulation(GalileoE5bI), codes, _galileo_e5b_i_secondary_code())
    GalileoE5bI(codes, lut)
end

function GalileoE5bQ()
    codes = widen_codes_to_storage(read_galileo_e5_codes(E5B_X1_TAP, E5B_Q_X2_TAP, E5B_Q_X2_INIT))
    # CS100_(n+50) to SVID n (Galileo OS SIS ICD v2.2 §3.5.2) — the upper half of
    # the shared CS100 table in `galileo/codes.jl`.
    secondary = _build_galileo_cs100_secondary(51:100)
    # The 100-chip per-SVID overlay is too long to bake (100·10230·1 > typemax(Int16)),
    # so it stays residual in the SignalLUT and is applied per primary period at gen time.
    lut = build_signal_lut(get_modulation(GalileoE5bQ), codes, PerPRNSecondaryCode(secondary))
    GalileoE5bQ(codes, secondary, lut)
end

# Shared interface (band, modulation, frequencies).

get_modulation(::Type{<:GalileoE5bI}) = LOC()
get_modulation(::Type{<:GalileoE5bQ}) = LOC()

# E5b-I data rides the in-phase carrier (default 0); the E5b-Q pilot is in
# quadrature and LEADS it by 90°, the same `(I + jQ)` convention E5a follows
# (OS SIS ICD v2.2 Eq. 2). See [`get_carrier_phase_offset`](@ref).
@inline get_carrier_phase_offset(::Type{<:GalileoE5bQ}) = π / 2

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

Galileo E5b sits at 1207.14 MHz, which BeiDou B2b shares exactly, so both report
[`E5b`](@ref) — band identity is by RF frequency, not by constellation, the same
reason Galileo E1 reports [`L1`](@ref).
"""
@inline get_band(::Type{<:GalileoE5bI}) = E5b()
@inline get_band(::Type{<:GalileoE5bQ}) = E5b()

# 50/50 I/Q power sharing of the E5b composite, whose −155.25 dBW total (Galileo
# OS SIS ICD v2.2, Table 13) is the unit here. See [`get_relative_power`](@ref).
@inline get_relative_power(::Type{<:GalileoE5bI}) = 0.5
@inline get_relative_power(::Type{<:GalileoE5bQ}) = 0.5

@inline get_signal_name(::Type{<:GalileoE5bI}) = "Galileo E5b-I"
@inline get_signal_name(::Type{<:GalileoE5bQ}) = "Galileo E5b-Q"

"""
$(SIGNATURES)

Get the code length for Galileo E5b (10230 chips, both components).
"""
@inline get_code_length(::Type{<:GalileoE5bI}) = E5_CODE_LENGTH
@inline get_code_length(::Type{<:GalileoE5bQ}) = E5_CODE_LENGTH

"""
$(SIGNATURES)

Get the code chipping rate for Galileo E5b (10.23 MHz, both components).
"""
@inline get_code_frequency(::Type{<:GalileoE5bI}) = 10_230_000Hz
@inline get_code_frequency(::Type{<:GalileoE5bQ}) = 10_230_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for Galileo E5b-I.

The E5b-I channel carries the Galileo I/NAV navigation message at 250 symbols/s
(Galileo OS SIS ICD v2.2, Table 5) — the same message E1-B carries, at the same
rate.

# Returns
- `Frequency`: 250 Hz
"""
@inline get_data_frequency(::Type{<:GalileoE5bI}) = 250Hz

"""
$(SIGNATURES)

Get the data symbol rate for Galileo E5b-Q.

The E5b-Q component is a dataless pilot.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_data_frequency(::Type{<:GalileoE5bQ}) = 0Hz

"""
$(SIGNATURES)

Get the secondary code for Galileo E5b-I.

CS4 is shared across all SVIDs: every 1 ms primary period is overlaid with one
chip of the 4-bit sequence `E` → `1110` (Galileo OS SIS ICD v2.2 §3.5.1),
giving a 4 ms tiered code.

# Returns
- [`SharedSecondaryCode`](@ref) of length 4
"""
@inline get_secondary_code(::GalileoE5bI) = _galileo_e5b_i_secondary_code()

# CS4 secondary, shared across SVIDs: hex E -> 1110 MSB-first, mapped 0 -> -1, 1 -> +1.
# Factored out so the `GalileoE5bI` constructor can build the embedded `SignalLUT`
# (which needs the secondary) before an instance exists.
@inline _galileo_e5b_i_secondary_code() =
    SharedSecondaryCode(Int8(1), Int8(1), Int8(1), Int8(-1))

"""
$(SIGNATURES)

Get the secondary code for Galileo E5b-Q.

Each SVID overlays its primary code with the 100-chip CS100₍ₙ₊₅₀₎ sequence
(Galileo OS SIS ICD v2.2 §3.5.2), giving a 100 ms tiered code.

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 100 × 50 secondary matrix
"""
@inline get_secondary_code(s::GalileoE5bQ) = PerPRNSecondaryCode(s.secondary_codes)
