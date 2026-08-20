#= Galileo code tables shared by more than one signal.

The CS100 secondary codes live here because three components draw from one
table: the E5a-Q pilot (Galileo OS SIS ICD v2.2 §3.5.2, Table 24) and the E6-C
pilot (E6-B/C Codes Technical Note §2.4, Table 2, which reproduces the ICD's
CS100_1-50 verbatim) take codes 1-50, and the E5b-Q pilot takes 51-100. So do
the E5 primary-code LFSR helpers below, shared by E5a and E5b. One table, one
generator, so the signals cannot drift apart. =#

#= CS100_1 .. CS100_100 as 25-hex-character (100-bit) strings, MSB-first (Galileo
OS SIS ICD v2.2 §3.5.1, Tables 22 and 23). The assignment (§3.5.2, Table 24) is
CS100_n to SVID n for E5a-Q and E6-C, and CS100_(n+50) to SVID n for E5b-Q — so
entries 1-50 serve the first two and 51-100 the third. These match the OS SIS ICD
v2.2 tables; note that some older tables (e.g. GNSS-SDR) carry stale values for
PRNs 37-47. =#
const GALILEO_CS100_HEX = [
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
    "CFF914EE3C6126A49FD5E5C94", "FC317C9A9BF8C6038B5CADAB3", "A2EAD74B6F9866E414393F239",
    "72F2B1180FA6B802CB84DF997", "13E3AE93BC52391D09E84A982", "77C04202B91B22C6D3469768E",
    "FEBC592DD7C69AB103D0BB29C", "0B494077E7C66FB6C51942A77", "DD0E321837A3D52169B7B577C",
    "43DEA90EA6C483E7990C3223F", "0366AB33F0167B6FA979DAE18", "99CCBBFAB1242CBE31E1BD52D",
    "A3466923CEFDF451EC0FCED22", "1A5271F22A6F9A8D76E79B7F0", "3204A6BB91B49D1A2D3857960",
    "32F83ADD43B599CBFB8628E5B", "3871FB0D89DB77553EB613CC1", "6A3CBDFF2D64D17E02773C645",
    "2BCD09889A1D7FC219F2EDE3B", "3E49467F4D4280B9942CD6F8C", "658E336DCFD9809F86D54A501",
    "ED4284F345170CF77268C8584", "29ECCE910D832CAF15E3DF5D1", "456CCF7FE9353D50E87A708FA",
    "FB757CC9E18CBC02BF1B84B9A", "5686229A8D98224BC426BC7FC", "700A2D325EA14C4B7B7AA8338",
    "1210A330B4D3B507D854CBA3F", "438EE410BD2F7DBCDD85565BA", "4B9764CC455AE1F61F7DA432B",
    "BF1F45FDDA3594ACF3C4CC806", "DA425440FE8F6E2C11B8EC1A4", "EE2C8057A7C16999AFA33FED1",
    "2C8BD7D8395C61DFA96243491", "391E4BB6BC43E98150CDDCADA", "399F72A9EADB42C90C3ECF7F0",
    "93031FDEA588F88E83951270C", "BA8061462D873705E95D5CB37", "D24188F88544EB121E963FD34",
    "D5F6A8BB081D8F383825A4DCA", "0FA4A205F0D76088D08EAF267", "272E909FAEBC65215E263E258",
    "3370F35A674922828465FC816", "54EF96116D4A0C8DB0E07101F", "DE347C7B27FADC48EF1826A2B",
    "01B16ECA6FC343AE08C5B8944", "1854DB743500EE94D8FC768ED", "28E40C684C87370CD0597FAB4",
    "5E42C19717093353BCAAF4033", "64310BAD8EB5B36E38646AF01",
]

const GALILEO_CS100_LENGTH = 100

"""
$(SIGNATURES)

Build a 100 × 50 Galileo CS100 secondary-code matrix from the CS100 code numbers
`code_numbers`, decoded MSB-first and mapped `0 -> -1`, `1 -> +1`. Column `n` is
the overlay of SVID `n`, i.e. the `n`-th code of the range.

[`GalileoE5aQ`](@ref) and [`GalileoE6C`](@ref) pass `1:50` (both are assigned
CS100ₙ to SVID `n`), [`GalileoE5bQ`](@ref) passes `51:100` (CS100₍ₙ₊₅₀₎ to SVID
`n`).
"""
function _build_galileo_cs100_secondary(code_numbers)
    codes = Matrix{Int8}(undef, GALILEO_CS100_LENGTH, length(code_numbers))
    for (prn, hex) in enumerate(@view GALILEO_CS100_HEX[code_numbers])
        chip = 0
        for c in hex
            nibble = parse(Int, string(c); base = 16)
            for shift = 3:-1:0
                chip += 1
                chip > GALILEO_CS100_LENGTH && break
                @inbounds codes[chip, prn] = Int8(2 * ((nibble >> shift) & 1) - 1)
            end
        end
    end
    codes
end

#= E5 primary code generation (Galileo OS SIS ICD v2.2 §3.3, §3.4.1).

Every E5 primary code — E5a-I, E5a-Q, E5b-I, E5b-Q — is the modulo-2 sum
(chip-wise product of ±1) of two length-2¹⁴ maximal-length sequences, truncated
to 10230 chips. Base register 1 starts all-ones for every SVID; base register 2
carries a per-SVID start value. Only the two feedback polynomials and the start
values differ between the four components (Table 16 and Tables 17-20), so the
generator is stated once here and each signal supplies its own parameters. =#

const E5_REGISTER_LENGTH = 14
const E5_CODE_LENGTH = 10230

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

"""
$(SIGNATURES)

Build the 10230 × 50 primary-code matrix of one Galileo E5 component from its
two octal feedback polynomials (`x1_tap`, `x2_tap`, ICD Table 16) and its table
of per-SVID base-register-2 start values (`x2_init_table`, ICD Tables 17-20).

Base register 1 does not depend on the SVID, so it is generated once and reused
for every PRN. Chips are mapped `0 -> -1`, `1 -> +1` to match the package-wide
convention (see [`GPSL5I`](@ref)).
"""
function read_galileo_e5_codes(x1_tap, x2_tap, x2_init_table)
    n = E5_REGISTER_LENGTH
    all_ones = (1 << n) - 1
    x1 = _e5_lfsr_bits(E5_CODE_LENGTH, all_ones, _e5_rev_reg(x1_tap >> 1, n), n)
    x2_tap_lsb = _e5_rev_reg(x2_tap >> 1, n)
    codes = Matrix{Int8}(undef, E5_CODE_LENGTH, length(x2_init_table))
    for (prn, x2_init) in enumerate(x2_init_table)
        x2 = _e5_lfsr_bits(E5_CODE_LENGTH, _e5_rev_reg(x2_init, n), x2_tap_lsb, n)
        @inbounds @views codes[:, prn] .= Int8(2) .* (x1 .⊻ x2) .- Int8(1)
    end
    codes
end
