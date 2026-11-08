"""
    GalileoE5aQP{C} <: AbstractGalileoSignal{C}

Galileo E5a-QP signal (E5a Quasi-Pilot), the short-code acquisition aid on the
E5a carrier.

New in OS SIS ICD Issue 2.2 (§2.3.1.4): a dataless BPSK(5) component at 5.115 Mcps on
1176.45 MHz ([`L5`](@ref)) whose primary code is only **330 chips** long, repeating every
2/31 ms, to enable low-complexity acquisition. The code is repeated 31 times within 2 ms
with no overlay ([`get_secondary_code`](@ref) is [`NoSecondaryCode`](@ref)).

The primary codes are a family of 40 sequences, each the XOR of three short codes of 15,
11 and 10 chips repeated up to 330 chips (ICD §3.4.2, Table 21). SVID `n` uses code `n`,
so PRNs 1-40 are supported.

!!! note
    E5a-QP is transmitted by a *subset* of the constellation (list at
    <https://www.gsc-europa.eu/galileo/services/galileo-open-service/quasipilot>), and the
    ICD states it will eventually be discontinued.

# Example
```julia
e5a_qp = GalileoE5aQP()
get_code_length(e5a_qp)      # 330
get_chip_rate(e5a_qp)        # 5115000 Hz
get_band(e5a_qp)             # L5()
```
"""
struct GalileoE5aQP{C<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

const E5A_QP_CODE_LENGTH = 330
const E5A_QP_NUM_PRNS = 40

#= E5a-QP generative short codes (Galileo OS SIS ICD Issue 2.2 §3.4.2, Table 21).

Each row is one code number's (15-chip, 11-chip, 10-chip) triple in the ICD's hex
notation: MSB-first, four chips per symbol, the last symbol zero-padded at the end (one
pad chip for 15 and 11 chips, two for 10). =#
const E5A_QP_SHORT_CODES = [
    ("335A", "60A", "384"), ("146E", "1C8", "238"), ("3876", "7AE", "4D4"),
    ("499A", "688", "064"), ("3548", "2C0", "588"), ("052E", "504", "070"),
    ("598C", "6E2", "530"), ("2E84", "064", "658"), ("01F2", "298", "6E4"),
    ("3162", "134", "478"), ("3E16", "148", "378"), ("2722", "712", "468"),
    ("4E9E", "5EC", "1BC"), ("20CA", "046", "580"), ("78FA", "74E", "6E0"),
    ("029E", "4D0", "4C0"), ("3D78", "106", "770"), ("4990", "5C4", "00C"),
    ("75D6", "034", "760"), ("03D6", "286", "20C"), ("456A", "114", "49C"),
    ("3AE0", "48E", "658"), ("4EBE", "2BE", "314"), ("7DB8", "67E", "7E4"),
    ("1056", "5DE", "270"), ("0E76", "620", "300"), ("401E", "63E", "1BC"),
    ("1DF2", "73E", "478"), ("691E", "37E", "194"), ("25EE", "328", "668"),
    ("4C08", "0A6", "3F8"), ("442C", "302", "580"), ("0272", "37C", "450"),
    ("42C4", "7A6", "32C"), ("25A2", "11A", "1D8"), ("5598", "4DE", "190"),
    ("2DD6", "508", "658"), ("1D64", "3BA", "5C4"), ("3E4C", "68E", "014"),
    ("7EC0", "60A", "5E8"),
]

"""
$(SIGNATURES)

Decode the first `n_chips` chips of an ICD hexadecimal code string into `0`/`1`
bits, MSB-first, discarding the zeros the ICD pads the last symbol with.
"""
function _e5a_qp_short_code(hex, n_chips)
    bits = Vector{Int8}(undef, n_chips)
    chip = 0
    for c in hex
        nibble = parse(Int, string(c); base = 16)
        for shift = 3:-1:0
            chip += 1
            chip > n_chips && break
            @inbounds bits[chip] = Int8((nibble >> shift) & 1)
        end
    end
    bits
end

"""
$(SIGNATURES)

Build the 330 × 40 Galileo E5a-QP primary code matrix from
`E5A_QP_SHORT_CODES`.

Each code is the XOR of its three short codes, repeated to 330 chips (22, 30 and 33
times for the 15-, 11- and 10-chip codes; Galileo OS SIS ICD Issue 2.2 §3.4.2). Chips are
mapped `0 -> -1`, `1 -> +1`.
"""
function read_galileo_e5a_qp_codes()
    codes = Matrix{Int8}(undef, E5A_QP_CODE_LENGTH, E5A_QP_NUM_PRNS)
    for (prn, (hex15, hex11, hex10)) in enumerate(E5A_QP_SHORT_CODES)
        short = (
            _e5a_qp_short_code(hex15, 15),
            _e5a_qp_short_code(hex11, 11),
            _e5a_qp_short_code(hex10, 10),
        )
        for chip = 1:E5A_QP_CODE_LENGTH
            bit = Int8(0)
            for s in short
                # Short codes repeat: chip `i` is chip `mod1(i, length)` of the short code.
                @inbounds bit ⊻= s[mod1(chip, length(s))]
            end
            @inbounds codes[chip, prn] = Int8(2) * bit - Int8(1)
        end
    end
    codes
end

function GalileoE5aQP()
    codes = read_galileo_e5a_qp_codes()
    lut = build_signal_lut(get_modulation(GalileoE5aQP), codes, NoSecondaryCode())
    GalileoE5aQP(codes, lut)
end

# Shared interface (band, modulation, frequencies).

get_modulation(::Type{<:GalileoE5aQP}) = BPSK()

@inline get_band(::Type{<:GalileoE5aQP}) = L5()

#= No carrier phase offset method, so the 0.0 default applies: the ICD (§2.3.1,
§2.3.1.4) states no phase relation to the E5a-I reference. =#

# A separate signal on the E5a carrier, scaled against the E5a composite (like `GPSL1CA`
# against L1C): −160.75 dBW against −155.25 dBW (Galileo OS SIS ICD Issue 2.2, Table 13).
@inline get_relative_power(::Type{<:GalileoE5aQP}) = 10.0^(-0.55)

@inline get_signal_name(::Type{<:GalileoE5aQP}) = "Galileo E5a-QP"

"""
$(SIGNATURES)

Get the code length for Galileo E5a-QP (330 chips — 2/31 ms at 5.115 Mcps).
"""
@inline get_code_length(::Type{<:GalileoE5aQP}) = E5A_QP_CODE_LENGTH

"""
$(SIGNATURES)

Get the nominal chip rate of Galileo E5a-QP (5.115 MHz, Galileo OS SIS ICD
Issue 2.2, Table 8).
"""
@inline get_chip_rate(::Type{<:GalileoE5aQP}) = 5_115_000Hz

"""
$(SIGNATURES)

Get the bit rate, as tracked, of Galileo E5a-QP.

E5a-QP carries no user data, so this is 0 Hz.
"""
@inline get_bit_rate(::Type{<:GalileoE5aQP}) = 0Hz

"""
$(SIGNATURES)

Get the secondary code for Galileo E5a-QP.

There is none: the 330-chip primary code is repeated 31 times within 2 ms without any
overlay (Galileo OS SIS ICD Issue 2.2, Table 15); returns [`NoSecondaryCode`](@ref).
"""
@inline get_secondary_code(::GalileoE5aQP) = NoSecondaryCode()

@inline get_supported_prns(::Type{<:GalileoE5aQP}) = 1:40
@inline get_secondary_code_length(::Type{<:GalileoE5aQP}) = 1
