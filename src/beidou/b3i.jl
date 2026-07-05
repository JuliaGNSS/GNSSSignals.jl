include("b3i_constants.jl")

"""
    BeiDouB3I{C, M} <: AbstractBeiDouSignal{C}

BeiDou B3I signal — the open-service signal on the BeiDou B3 frequency
(1268.52 MHz, reported as the [`B3I`](@ref) band).

BPSK-modulated 10230-chip ranging code at 10.23 Mcps (a Gold code). As on
[`BeiDouB1I`](@ref), the 20-bit Neuman-Hoffman secondary code (NH20) overlays
the ranging code only on the MEO/IGSO satellites (PRN 6-58, D1 message, 1 kHz),
giving a 20 ms tiered code; the GEO satellites (PRN 1-5 and 59-63, D2 message)
carry no secondary code (BDS-SIS-ICD-B3I-1.0 §5.2.1).

The ranging code is the modulo-2 sum of two 13-stage shift registers G1 and G2
(BDS-SIS-ICD-B3I-1.0 §4.3). G1 (`X¹³+X⁴+X³+X+1`) starts all-ones and is reset to
that state after 8190 chips (short-cycling it to a 8190-chip period); G2
(`X¹³+X¹²+X¹⁰+X⁹+X⁷+X⁶+X⁵+X+1`) starts from a per-SVID initial phase
(`B3I_G2_INITIAL_PHASES`, ICD Table 4-1) and runs its full 8191-chip period.
Both registers reset at the start of each 1 ms period. PRNs 1-63 are supported.

# Example
```julia
b3i = BeiDouB3I()
get_code_length(b3i)            # 10230
get_secondary_code_length(b3i)  # 20
get_band(b3i)                   # B3I()
```
"""
struct BeiDouB3I{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractBeiDouSignal{C}
    codes::C
    secondary_codes::M    # 20 × 63 NH20 overlay: NH20 for D1 (MEO/IGSO) PRNs, all-ones for D2 (GEO)
    lut::SignalLUT        # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

# G1(X) = X¹³ + X⁴ + X³ + X + 1 ; G2(X) = X¹³ + X¹² + X¹⁰ + X⁹ + X⁷ + X⁶ + X⁵ + X + 1
const B3I_G1_FEEDBACK = (1, 3, 4, 13)
const B3I_G2_FEEDBACK = (1, 5, 6, 7, 9, 10, 12, 13)
const B3I_CODE_LENGTH = 10230
const B3I_RESET1_AT = 8190
const B3I_NUM_PRNS = 63

function read_beidou_b3i_codes()
    init1 = fill(Int8(1), 13)   # G1 all-ones start
    codes = Matrix{Int8}(undef, B3I_CODE_LENGTH, B3I_NUM_PRNS)
    for prn = 1:B3I_NUM_PRNS
        init2 = _beidou_state(B3I_G2_INITIAL_PHASES[prn])
        codes[:, prn] = _beidou_gold_code(
            13, B3I_CODE_LENGTH, B3I_G1_FEEDBACK, B3I_G2_FEEDBACK,
            (13,), (13,), init1, init2, B3I_RESET1_AT,
        )
    end
    codes
end

function BeiDouB3I()
    codes = widen_codes_to_storage(read_beidou_b3i_codes())
    secondary = _beidou_nh20_matrix(B3I_NUM_PRNS)
    lut = build_signal_lut(get_modulation(BeiDouB3I), codes, PerPRNSecondaryCode(secondary))
    BeiDouB3I(codes, secondary, lut)
end

get_modulation(::Type{<:BeiDouB3I}) = LOC()
@inline get_modulation(::BeiDouB3I) = LOC()

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

BeiDou B3I is on the 1268.52 MHz B3 frequency, so this returns [`B3I`](@ref).
"""
@inline get_band(::Type{<:BeiDouB3I}) = B3I()

"""
$(SIGNATURES)

Get the human-readable signal name.
"""
get_signal_name(::BeiDouB3I) = "BeiDou B3I"

"""
$(SIGNATURES)

Get the code length for BeiDou B3I (10230 chips).
"""
@inline get_code_length(::Type{<:BeiDouB3I}) = B3I_CODE_LENGTH

"""
$(SIGNATURES)

Get the code chipping rate for BeiDou B3I (10.23 MHz).
"""
@inline get_code_frequency(::Type{<:BeiDouB3I}) = 10_230_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for BeiDou B3I (D1 navigation message, 50 bps).

# Returns
- `Frequency`: 50 Hz
"""
@inline get_data_frequency(::Type{<:BeiDouB3I}) = 50Hz

"""
$(SIGNATURES)

Get the secondary (Neuman-Hoffman NH20) code for BeiDou B3I.

Same per-SVID NH20 overlay as [`BeiDouB1I`](@ref) (BDS-SIS-ICD-B3I-1.0 §5.2.1):
NH20 on the MEO/IGSO PRNs (6-58, D1) for a 20 ms tiered code, and no overlay
(all-ones column) on the GEO PRNs (1-5, 59-63, D2).

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 20 × 63 overlay matrix
"""
@inline get_secondary_code(s::BeiDouB3I) = PerPRNSecondaryCode(s.secondary_codes)
