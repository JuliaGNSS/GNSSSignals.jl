"""
    GalileoE5aI{C} <: AbstractGNSSSignal{C}

Galileo E5a-I signal (the in-phase, data-carrying component of Galileo E5a).

10230-chip primary code at 10.23 Mcps on the E5 band — which shares the GPS
L5 carrier at 1176.45 MHz, so [`get_band`](@ref) returns [`L5`](@ref). A 20-bit
secondary code (CS20, the same for every SVID) overlays the data channel,
giving a 20 ms tiered code.

Galileo E5a is one half of the wideband E5 AltBOC(15,10) signal, but like
GNSS-SDR and PocketSDR this implementation models the E5a sideband on its own
as BPSK(10) (modulation [`LOC`](@ref)) — the form used when E5a is acquired
and tracked as an independent signal.

The primary code is generated from two 14-stage shift registers per the
Galileo OS SIS ICD (§3.5): a common base register 1 (all-ones start, feedback
`40503₈`) XOR'd with a per-SVID base register 2 (feedback `50661₈`, start
values `E5A_I_X2_INIT`), truncated to 10230 chips. PRNs 1-50 are supported.

# Example
```julia
e5a_i = GalileoE5aI()
get_code_length(e5a_i)            # 10230
get_secondary_code_length(e5a_i)  # 20
get_band(e5a_i)                   # L5()
```
"""
struct GalileoE5aI{C<:AbstractMatrix} <: AbstractGNSSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GalileoE5aQ{C, M} <: AbstractGNSSSignal{C}

Galileo E5a-Q signal (the quadrature, dataless pilot component of Galileo E5a).

10230-chip primary code at 10.23 Mcps on the E5 band (1176.45 MHz, reported as
[`L5`](@ref)), overlaid with a 100-bit per-SVID secondary code (CS100) for a
100 ms tiered code. As the pilot component it carries no navigation data, so
[`get_data_frequency`](@ref) returns 0 Hz.

Like [`GalileoE5aI`](@ref) the E5a sideband is modelled on its own as
BPSK(10) ([`LOC`](@ref)). The primary code uses the same generator as E5a-I
(base register 1 feedback `40503₈`, base register 2 feedback `50661₈`) with
the E5a-Q per-SVID register-2 start values `E5A_Q_X2_INIT`. PRNs 1-50 are
supported.

The struct stores the CS100 overlay matrix, exposed via
[`get_secondary_code`](@ref) as a [`PerPRNSecondaryCode`](@ref).

# Example
```julia
e5a_q = GalileoE5aQ()
get_code_length(e5a_q)            # 10230
get_secondary_code_length(e5a_q)  # 100
get_data_frequency(e5a_q)         # 0 Hz
```
"""
struct GalileoE5aQ{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractGNSSSignal{C}
    codes::C
    secondary_codes::M    # 100 × 50 Int8 ±1 matrix, exposed via PerPRNSecondaryCode
    lut::SignalLUT        # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

#= E5a primary code generation (Galileo OS SIS ICD §3.5).

Each E5 primary code is the modulo-2 sum (chip-wise product of ±1) of two
length-2¹⁴ maximal-length sequences, truncated to 10230 chips. Base register 1
is common to every SVID (all-ones start); base register 2 carries a per-SVID
start value. The taps below are the ICD feedback polynomials in octal; E5a-I
and E5a-Q share the same taps and differ only in the register-2 start values. =#

const E5A_X1_TAP = 0o40503   # base register 1 feedback polynomial
const E5A_X2_TAP = 0o50661   # base register 2 feedback polynomial

# E5a-I base register 2 start values (octal), PRN 1-50, used with E5A_X2_TAP.
const E5A_I_X2_INIT = [
    0o30305, 0o14234, 0o27213, 0o20577, 0o23312, 0o33463, 0o15614, 0o12537, 0o01527, 0o30236,
    0o27344, 0o07272, 0o36377, 0o17046, 0o06434, 0o15405, 0o24252, 0o11631, 0o24776, 0o00630,
    0o11560, 0o17272, 0o27445, 0o31702, 0o13012, 0o14401, 0o34727, 0o22627, 0o30623, 0o27256,
    0o01520, 0o14211, 0o31465, 0o22164, 0o33516, 0o02737, 0o21316, 0o35425, 0o35633, 0o24655,
    0o14054, 0o27027, 0o06604, 0o31455, 0o34465, 0o25273, 0o20763, 0o31721, 0o17312, 0o13277,
]

# E5a-Q base register 2 start values (octal), PRN 1-50 — the E5a-Q counterpart.
const E5A_Q_X2_INIT = [
    0o25652, 0o05142, 0o24723, 0o31751, 0o27366, 0o24660, 0o33655, 0o27450, 0o07626, 0o01705,
    0o12717, 0o32122, 0o16075, 0o16644, 0o37556, 0o02477, 0o02265, 0o06430, 0o25046, 0o12735,
    0o04262, 0o11230, 0o00037, 0o06137, 0o04312, 0o20606, 0o11162, 0o22252, 0o30533, 0o24614,
    0o07767, 0o32705, 0o05052, 0o27553, 0o03711, 0o02041, 0o34775, 0o05274, 0o37356, 0o16205,
    0o36270, 0o06600, 0o26773, 0o17375, 0o35267, 0o36255, 0o12044, 0o26442, 0o21621, 0o25411,
]

"""
$(SIGNATURES)

Reverse the low `n` bits of `register` (bit 0 ↔ bit n-1).

The ICD numbers shift-register stages and feedback taps MSB-first, whereas the
generation loop in `_e5_lfsr_bits` shifts LSB-first. This converts a tap mask
or start value between the two conventions.
"""
function _e5_rev_reg(register, n)
    reversed = 0
    for i = 0:n-1
        reversed = (reversed << 1) | ((register >> i) & 1)
    end
    reversed
end

"""
$(SIGNATURES)

Run a length-`n` Fibonacci LFSR for `count` chips and return its output bit
sequence (`0`/`1`). `register` is the initial state and `tap` the feedback
mask (both LSB-first). Each step outputs the LSB, then shifts right and feeds
the parity of `register & tap` into the top stage.
"""
function _e5_lfsr_bits(count, register, tap, n)
    bits = Vector{Int8}(undef, count)
    mask = (1 << n) - 1
    for i = 1:count
        bits[i] = register & 1
        feedback = count_ones(register & tap) & 1
        register = ((feedback << (n - 1)) | (register >> 1)) & mask
    end
    bits
end

# Base register 1 (common to every SVID, all-ones start) and the per-SVID
# base register 2; both taps are converted to the LSB-first loop convention.
_e5_x1_bits(count) = _e5_lfsr_bits(count, 0b11111111111111, _e5_rev_reg(E5A_X1_TAP >> 1, 14), 14)
_e5_x2_bits(count, x2_init) =
    _e5_lfsr_bits(count, _e5_rev_reg(x2_init, 14), _e5_rev_reg(E5A_X2_TAP >> 1, 14), 14)

"""
$(SIGNATURES)

Build the 10230 × 50 Galileo E5a primary code matrix for the given table of
per-SVID base-register-2 start values (`E5A_I_X2_INIT` or `E5A_Q_X2_INIT`).

The common base register 1 is generated once and reused for every PRN. Chips
are mapped `0 -> -1`, `1 -> +1` to match the package-wide convention (see
[`GPSL5I`](@ref)).
"""
function read_galileo_e5a_codes(x2_init_table)
    code_length = 10230
    x1 = _e5_x1_bits(code_length)
    codes = Matrix{Int8}(undef, code_length, length(x2_init_table))
    for (prn, x2_init) in enumerate(x2_init_table)
        x2 = _e5_x2_bits(code_length, x2_init)
        @inbounds @views codes[:, prn] .= Int8(2) .* (x1 .⊻ x2) .- Int8(1)
    end
    codes
end

function GalileoE5aI()
    codes = widen_codes_to_storage(read_galileo_e5a_codes(E5A_I_X2_INIT))
    lut = build_signal_lut(get_modulation(GalileoE5aI), codes, _galileo_e5a_i_secondary_code())
    GalileoE5aI(codes, lut)
end

#= E5a-I secondary code (CS20, Galileo OS SIS ICD §3.5).

CS20 = 0x842E9 -> bits 1000 0100 0010 1110 1001 (MSB-first), mapped 0 -> -1,
1 -> +1. The same 20-bit overlay is used for the E5a-I data channel of every
SVID, giving a 20 ms tiered code. =#
const E5A_I_SECONDARY_CHIPS = (
    Int8(1), Int8(-1), Int8(-1), Int8(-1), Int8(-1), Int8(1), Int8(-1), Int8(-1), Int8(-1), Int8(-1),
    Int8(1), Int8(-1), Int8(1), Int8(1), Int8(1), Int8(-1), Int8(1), Int8(-1), Int8(-1), Int8(1),
)

#= E5a-Q secondary codes (CS100, Galileo OS SIS ICD §3.5, one per SVID).

CS100_1 .. CS100_50 as 25-hex-character (100-bit) strings, MSB-first. These
match the OS SIS ICD v2.2 secondary-code table; note that some older tables
(e.g. GNSS-SDR) carry stale values for PRNs 37-47. =#
const E5A_Q_SECONDARY_HEX = [
    "83F6F69D8F6E15411FB8C9B1C", "66558BD3CE0C7792E83350525", "59A025A9C1AF0651B779A8381",
    "D3A32640782F7B18E4DF754B7", "B91FCAD7760C218FA59348A93", "BAC77E933A779140F094FBF98",
    "537785DE280927C6B58BA6776", "EFCAB4B65F38531ECA22257E2", "79F8CAE838475EA5584BEFC9B",
    "CA5170FEA3A810EC606B66494", "1FC32410652A2C49BD845E567", "FE0A9A7AFDAC44E42CB95D261",
    "B03062DC2B71995D5AD8B7DBE", "F6C398993F598E2DF4235D3D5", "1BB2FB8B5BF24395C2EF3C5A1",
    "2F920687D238CC7046EF6AFC9", "34163886FC4ED7F2A92EFDBB8", "66A872CE47833FB2DFD5625AD",
    "99D5A70162C920A4BB9DE1CA8", "81D71BD6E069A7ACCBEDC66CA", "A654524074A9E6780DB9D3EC6",
    "C3396A101BEDAF623CFC5BB37", "C3D4AB211DF36F2111F2141CD", "3DFF25EAE761739265AF145C1",
    "994909E0757D70CDE389102B5", "B938535522D119F40C25FDAEC", "C71AB549C0491537026B390B7",
    "0CDB8C9E7B53F55F5B0A0597B", "61C5FA252F1AF81144766494F", "626027778FD3C6BB4BAA7A59D",
    "E745412FF53DEBD03F1C9A633", "3592AC083F3175FA724639098", "52284D941C3DCAF2721DDB1FD",
    "73B3D8F0AD55DF4FE814ED890", "94BF16C83BD7462F6498E0282", "A8C3DE1AC668089B0B45B3579",
    "E23FFC2DD2C14388AD8D6BEC8", "F2AC871CDF89DDC06B5960D2B", "06191EC1F622A77A526868BA1",
    "22D6E2A768E5F35FFC8E01796", "25310A06675EB271F2A09EA1D", "9F7993C621D4BEC81A0535703",
    "D62999EACF1C99083C0B4A417", "F665A7EA441BAA4EA0D01078C", "46F3D3043F24CDEABD6F79543",
    "E2E3E8254616BD96CEFCA651A", "E548231A82F9A01A19DB5E1B2", "265C7F90A16F49EDE2AA706C8",
    "364A3A9EB0F0481DA0199D7EA", "9810A7A898961263A0F749F56",
]

"""
$(SIGNATURES)

Build the 100 × 50 Galileo E5a-Q secondary (CS100) code matrix from
`E5A_Q_SECONDARY_HEX`. Each column is one SVID's 100-chip overlay, decoded
MSB-first and mapped `0 -> -1`, `1 -> +1`.
"""
function _build_galileo_e5a_q_secondary()
    code_length = 100
    codes = Matrix{Int8}(undef, code_length, length(E5A_Q_SECONDARY_HEX))
    for (prn, hex) in enumerate(E5A_Q_SECONDARY_HEX)
        chip = 0
        for c in hex
            nibble = parse(Int, string(c); base = 16)
            for shift = 3:-1:0
                chip += 1
                chip > code_length && break
                @inbounds codes[chip, prn] = Int8(2 * ((nibble >> shift) & 1) - 1)
            end
        end
    end
    codes
end

function GalileoE5aQ()
    codes = widen_codes_to_storage(read_galileo_e5a_codes(E5A_Q_X2_INIT))
    secondary = _build_galileo_e5a_q_secondary()
    # The 100-chip per-SVID CS100 overlay is too long to bake (100·10230·1 > typemax(Int16)),
    # so it stays residual in the SignalLUT and is applied per primary period at gen time.
    lut = build_signal_lut(get_modulation(GalileoE5aQ), codes, PerPRNSecondaryCode(secondary))
    GalileoE5aQ(codes, secondary, lut)
end

# Shared interface (band, modulation, frequencies).

get_modulation(::Type{<:GalileoE5aI}) = LOC()
@inline get_modulation(::GalileoE5aI) = LOC()
get_modulation(::Type{<:GalileoE5aQ}) = LOC()
@inline get_modulation(::GalileoE5aQ) = LOC()

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

Galileo E5a shares the GPS L5 carrier frequency (1176.45 MHz), so this returns
[`L5`](@ref) — band identity is by RF, not by ICD label.

# Examples
```julia-repl
julia> get_band(GalileoE5aI())
L5()
```
"""
@inline get_band(::GalileoE5aI) = L5()
@inline get_band(::GalileoE5aQ) = L5()

"""
$(SIGNATURES)

Get the human-readable signal name.

# Examples
```julia-repl
julia> get_signal_name(GalileoE5aI())
"Galileo E5a-I"
```
"""
get_signal_name(::GalileoE5aI) = "Galileo E5a-I"
get_signal_name(::GalileoE5aQ) = "Galileo E5a-Q"

"""
$(SIGNATURES)

Get the code length for Galileo E5a (10230 chips, both components).
"""
@inline get_code_length(::GalileoE5aI) = 10230
@inline get_code_length(::GalileoE5aQ) = 10230

"""
$(SIGNATURES)

Get the code chipping rate for Galileo E5a (10.23 MHz, both components).
"""
@inline get_code_frequency(::GalileoE5aI) = 10_230_000Hz
@inline get_code_frequency(::GalileoE5aQ) = 10_230_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for Galileo E5a-I.

The E5a-I channel carries the Galileo F/NAV navigation message at a 50
symbols/s rate (25 bps with rate-1/2 convolutional coding).

# Returns
- `Frequency`: 50 Hz
"""
@inline get_data_frequency(::GalileoE5aI) = 50Hz

"""
$(SIGNATURES)

Get the data symbol rate for Galileo E5a-Q.

The E5a-Q component is a dataless pilot.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_data_frequency(::GalileoE5aQ) = 0Hz

"""
$(SIGNATURES)

Get the secondary code for Galileo E5a-I.

CS20 is shared across all SVIDs: every 1 ms primary period is overlaid with
one chip of the 20-bit sequence (Galileo OS SIS ICD §3.5), giving a 20 ms
tiered code.

# Returns
- [`SharedSecondaryCode`](@ref) of length 20
"""
@inline get_secondary_code(::GalileoE5aI) = _galileo_e5a_i_secondary_code()

# CS20 secondary, shared across SVIDs. Factored out so the `GalileoE5aI` constructor can build
# the embedded `SignalLUT` (which needs the secondary) before an instance exists.
@inline _galileo_e5a_i_secondary_code() = SharedSecondaryCode(E5A_I_SECONDARY_CHIPS...)

"""
$(SIGNATURES)

Get the secondary code for Galileo E5a-Q.

Each SVID overlays its primary code with a 100-chip CS100 sequence (Galileo OS
SIS ICD §3.5), giving a 100 ms tiered code.

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 100 × 50 CS100 matrix
"""
@inline get_secondary_code(s::GalileoE5aQ) = PerPRNSecondaryCode(s.secondary_codes)
