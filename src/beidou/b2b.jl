include("b2b_constants.jl")

"""
    BeiDouB2bI{C} <: AbstractBeiDouSignal{C}

BeiDou B2b_I signal — the open-service data signal on the BeiDou B2b frequency
(1207.14 MHz, reported as the [`B2b`](@ref) band).

BPSK(10)-modulated 10230-chip ranging code at 10.23 Mcps (a Gold code); the
B-CNAV3 navigation message is broadcast at 1000 symbols/s. B2b_I has no
secondary code (its 1 ms primary period already matches the symbol period).

The ranging code is the modulo-2 sum of two 13-stage shift registers
(BDS-SIS-ICD-B2b-1.0 §5): register 1 (`1+X+X⁹+X¹⁰+X¹³`) starts all-ones and is
reset after 8190 chips; register 2 (`1+X³+X⁴+X⁶+X⁹+X¹²+X¹³`) starts from a
per-SVID initial value (`B2B_REG2_INIT`, ICD Table 5-1) and runs its full
period. The ICD defines 53 codes, for **PRN 6-58**; other PRN indices are not
defined and generate an all-zero code.

# Example
```julia
b2b = BeiDouB2bI()
get_code_length(b2b)     # 10230
get_band(b2b)            # B2b()
gen_code(20460, b2b, 6, 20.46e6 * u"Hz")   # PRN 6 (first defined PRN)
```
"""
struct BeiDouB2bI{C<:AbstractMatrix} <: AbstractBeiDouSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

# g1(x) = 1 + x + x⁹ + x¹⁰ + x¹³ ; g2(x) = 1 + x³ + x⁴ + x⁶ + x⁹ + x¹² + x¹³
const B2B_G1_FEEDBACK = (1, 9, 10, 13)
const B2B_G2_FEEDBACK = (3, 4, 6, 9, 12, 13)
const B2B_CODE_LENGTH = 10230
const B2B_RESET1_AT = 8190
const B2B_NUM_PRNS = 63          # BDS PRN space; only 6..58 are defined by the ICD

function read_beidou_b2b_codes()
    init1 = fill(Int8(1), 13)
    # Undefined PRNs (not in the ICD) stay all-zero — an obviously-invalid code.
    codes = zeros(Int8, B2B_CODE_LENGTH, B2B_NUM_PRNS)
    for (prn, init_str) in B2B_REG2_INIT
        init2 = _beidou_state(init_str)
        codes[:, prn] = _beidou_gold_code(
            13, B2B_CODE_LENGTH, B2B_G1_FEEDBACK, B2B_G2_FEEDBACK,
            (13,), (13,), init1, init2, B2B_RESET1_AT,
        )
    end
    codes
end

function BeiDouB2bI()
    codes = widen_codes_to_storage(read_beidou_b2b_codes())
    lut = build_signal_lut(get_modulation(BeiDouB2bI), codes, NoSecondaryCode())
    BeiDouB2bI(codes, lut)
end

get_modulation(::Type{<:BeiDouB2bI}) = LOC()
@inline get_modulation(::BeiDouB2bI) = LOC()

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

BeiDou B2b is on the 1207.14 MHz B2b frequency, so this returns [`B2b`](@ref).
"""
@inline get_band(::Type{<:BeiDouB2bI}) = B2b()

"""
$(SIGNATURES)

Get the human-readable signal name.
"""
get_signal_name(::BeiDouB2bI) = "BeiDou B2b-I"

"""
$(SIGNATURES)

Get the code length for BeiDou B2b_I (10230 chips).
"""
@inline get_code_length(::Type{<:BeiDouB2bI}) = B2B_CODE_LENGTH

"""
$(SIGNATURES)

Get the code chipping rate for BeiDou B2b_I (10.23 MHz).
"""
@inline get_code_frequency(::Type{<:BeiDouB2bI}) = 10_230_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for BeiDou B2b_I.

The B-CNAV3 message is broadcast at 1000 symbols/s (500 bps with rate-1/2
coding); `get_data_frequency` returns the broadcast symbol rate.

# Returns
- `Frequency`: 1000 Hz
"""
@inline get_data_frequency(::Type{<:BeiDouB2bI}) = 1000Hz

"""
$(SIGNATURES)

Get the secondary code for BeiDou B2b_I.

B2b_I has no secondary code.

# Returns
- [`NoSecondaryCode`](@ref)
"""
@inline get_secondary_code(::BeiDouB2bI) = NoSecondaryCode()
