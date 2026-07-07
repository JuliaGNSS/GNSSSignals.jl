include("b1i_constants.jl")

"""
    BeiDouB1I{C, M} <: AbstractBeiDouSignal{C}

BeiDou B1I signal — the legacy open-service signal on the BeiDou B1 frequency
(1561.098 MHz, reported as the [`B1I`](@ref) band).

BPSK-modulated 2046-chip ranging code at 2.046 Mcps (a balanced Gold code). The
20-bit Neuman-Hoffman secondary code (NH20, at 1 kHz) overlays the ranging code
on the **MEO/IGSO** satellites (PRN 6-58), which broadcast the D1 navigation
message (50 bps), giving a 20 ms tiered code. The **GEO** satellites (PRN 1-5
and 59-63) broadcast the faster D2 message (500 bps) with **no** NH overlay, so
[`get_secondary_code`](@ref) returns a per-SVID overlay that is NH20 for the
MEO/IGSO PRNs and a no-op (all-ones) column for the GEO PRNs (BDS-SIS-ICD-B1I-3.0
Table 4-1 / §5.2.1).

The ranging code is the modulo-2 sum of two 11-stage shift registers G1 and G2
(BDS-SIS-ICD-B1I-3.0 §4.3): both start at `01010101010` and are reset every
2046 chips (the length-2047 m-sequences truncated by one chip). The per-SVID
code is selected by tapping different G2 stages (`B1I_G2_PHASE_SELECT`,
ICD Table 4-1). PRNs 1-63 are supported.

# Example
```julia
b1i = BeiDouB1I()
get_code_length(b1i)            # 2046
get_secondary_code_length(b1i)  # 20
get_band(b1i)                   # B1I()
get_code_frequency(b1i)         # 2046000 Hz
```
"""
struct BeiDouB1I{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractBeiDouSignal{C}
    codes::C
    secondary_codes::M    # 20 × 63 NH20 overlay: NH20 for D1 (MEO/IGSO) PRNs, all-ones for D2 (GEO)
    lut::SignalLUT        # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

# G1(X) = 1 + X + X⁷ + X⁸ + X⁹ + X¹⁰ + X¹¹ ; G2(X) = 1 + X + X² + X³ + X⁴ + X⁵ + X⁸ + X⁹ + X¹¹
const B1I_G1_FEEDBACK = (1, 7, 8, 9, 10, 11)
const B1I_G2_FEEDBACK = (1, 2, 3, 4, 5, 8, 9, 11)
const B1I_INITIAL_PHASE = "01010101010"     # G1 and G2 start state (ICD §4.3)
const B1I_CODE_LENGTH = 2046
const B1I_NUM_PRNS = 63

function read_beidou_b1i_codes()
    init = _beidou_state(B1I_INITIAL_PHASE)
    codes = Matrix{Int8}(undef, B1I_CODE_LENGTH, B1I_NUM_PRNS)
    for prn = 1:B1I_NUM_PRNS
        # G1 output is stage 11; G2 output is the phase-selected modulo-2 sum
        # of the tapped stages. No mid-period reset (reset1_at = 0).
        codes[:, prn] = _beidou_gold_code(
            11, B1I_CODE_LENGTH, B1I_G1_FEEDBACK, B1I_G2_FEEDBACK,
            (11,), B1I_G2_PHASE_SELECT[prn], init, init, 0,
        )
    end
    codes
end

function BeiDouB1I()
    codes = widen_codes_to_storage(read_beidou_b1i_codes())
    secondary = _beidou_nh20_matrix(B1I_NUM_PRNS)
    lut = build_signal_lut(get_modulation(BeiDouB1I), codes, PerPRNSecondaryCode(secondary))
    BeiDouB1I(codes, secondary, lut)
end

get_modulation(::Type{<:BeiDouB1I}) = LOC()
@inline get_modulation(::BeiDouB1I) = LOC()

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

BeiDou B1I is on the 1561.098 MHz B1 frequency, so this returns [`B1I`](@ref).

# Examples
```julia-repl
julia> get_band(BeiDouB1I())
B1I()
```
"""
@inline get_band(::Type{<:BeiDouB1I}) = B1I()

"""
$(SIGNATURES)

Get the human-readable signal name.

# Examples
```julia-repl
julia> get_signal_name(BeiDouB1I())
"BeiDou B1I"
```
"""
get_signal_name(::BeiDouB1I) = "BeiDou B1I"

"""
$(SIGNATURES)

Get the code length for BeiDou B1I (2046 chips).
"""
@inline get_code_length(::Type{<:BeiDouB1I}) = B1I_CODE_LENGTH

"""
$(SIGNATURES)

Get the code chipping rate for BeiDou B1I (2.046 MHz).
"""
@inline get_code_frequency(::Type{<:BeiDouB1I}) = 2_046_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for BeiDou B1I.

The D1 navigation message (MEO/IGSO satellites) is broadcast at 50 bps; the GEO
satellites broadcast the D2 message at 500 bps. This returns the D1 rate.

# Returns
- `Frequency`: 50 Hz
"""
@inline get_data_frequency(::Type{<:BeiDouB1I}) = 50Hz

"""
$(SIGNATURES)

Get the secondary (Neuman-Hoffman NH20) code for BeiDou B1I.

Per BDS-SIS-ICD-B1I-3.0 §5.2.1 the NH20 overlay is applied only on the MEO/IGSO
satellites (PRN 6-58, D1 message), giving a 20 ms tiered code; the GEO
satellites (PRN 1-5 and 59-63, D2 message) carry no secondary code. This is
returned as a per-SVID overlay whose GEO columns are all-ones (a no-op, so a
GEO PRN's tiered code equals its primary code).

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 20 × 63 overlay matrix
"""
@inline get_secondary_code(s::BeiDouB1I) = PerPRNSecondaryCode(s.secondary_codes)
