@testset "Galileo E5b-I" begin
    e5b_i = GalileoE5bI()
    @test @inferred(get_band(e5b_i)) == E5b()
    @test @inferred(get_center_frequency(e5b_i)) == 1.20714e9Hz
    @test @inferred(get_code_length(e5b_i)) == 10230
    @test @inferred(get_secondary_code_length(e5b_i)) == 4
    @test @inferred(get_secondary_code(e5b_i)) isa SharedSecondaryCode{4}
    @test @inferred(get_data_frequency(e5b_i)) == 250Hz     # I/NAV symbol rate
    @test @inferred(get_code_frequency(e5b_i)) == 10230e3Hz
    @test get_signal_name(e5b_i) == "Galileo E5b-I"
    @test @inferred(get_modulation(e5b_i)) == GNSSSignals.LOC()
    @test get_code_type(e5b_i) === Int16

    # PRN 1's first chip is 0 -> -1 (ICD Table 19 lists C5BEA1 for code 1, i.e.
    # 1100 0101 ... after the 0 -> -1 mapping the first chip is +1); CS4[0] = +1,
    # so the tiered chip equals the primary chip here.
    @test @inferred(get_code(e5b_i, 0, 1)) == 1
    @test @inferred(get_code(e5b_i, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(e5b_i, 0.0, 1)) == 1

    # BPSK(10) spectrum, and 10.23 MHz against the 1207.14 MHz E5b carrier.
    @test GNSSSignals.get_code_factor(e5b_i) == 1
    @test get_code_spectrum(e5b_i, 0) ≈ 1.0Hz / get_code_frequency(e5b_i)
    @test get_code_center_frequency_ratio(e5b_i) ≈ 1 / 118
end

@testset "Galileo E5b-Q" begin
    e5b_q = GalileoE5bQ()
    @test @inferred(get_band(e5b_q)) == E5b()
    @test @inferred(get_center_frequency(e5b_q)) == 1.20714e9Hz
    @test @inferred(get_code_length(e5b_q)) == 10230
    @test @inferred(get_secondary_code_length(e5b_q)) == 100
    @test @inferred(get_secondary_code(e5b_q)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(e5b_q)) == 0Hz       # dataless pilot
    @test @inferred(get_code_frequency(e5b_q)) == 10230e3Hz
    @test get_signal_name(e5b_q) == "Galileo E5b-Q"
    @test @inferred(get_modulation(e5b_q)) == GNSSSignals.LOC()
    @test get_code_type(e5b_q) === Int16

    @test GNSSSignals.get_code_factor(e5b_q) == 1
    @test get_code_spectrum(e5b_q, 0) ≈ 1.0Hz / get_code_frequency(e5b_q)
    @test get_code_center_frequency_ratio(e5b_q) ≈ 1 / 118
end

# The hexadecimal "Initial Sequence" column of Galileo OS SIS ICD v2.2 Tables 19
# (E5b-I) and 20 (E5b-Q): the first 24 chips of every primary code, MSB-first.
# Pinning all 50 per component checks both feedback polynomials (Table 16) and
# every base-register-2 start value the source carries, independently of the
# full-length fixtures below.
const E5B_I_FIRST24_HEX = [
    "C5BEA1", "4F6248", "FD5488", "86277B", "9E39D5", "EA7EDE",
    "F28321", "4FB0C9", "F0AB64", "79833B", "EC2D91", "409B11",
    "397E16", "1E0FCD", "CB2F5A", "C1079A", "2D9BC6", "BC5146",
    "F848B0", "4F01E8", "B16C9B", "B2827D", "16C809", "7B570F",
    "1969C0", "512FA9", "73F36B", "2D5317", "EC8390", "380374",
    "46B4DE", "6C01D9", "0E0BB6", "2708C7", "265B55", "A68E1C",
    "C3916E", "BDC595", "1327D0", "EA921F", "D45869", "3EB98A",
    "9FDE16", "2A04CA", "3CA56F", "03928A", "4FB5B9", "101EC7",
    "C91D4F", "AFC22B",
]
const E5B_Q_FIRST24_HEX = [
    "E49AF0", "CE701F", "54B709", "641AB1", "FBD0AE", "0D8BC9",
    "805FA5", "D86BA0", "A7E921", "067E55", "3E4B58", "7D82FB",
    "E33BC2", "31372C", "46676F", "D2613E", "EB443C", "3FD0B1",
    "FCB7CF", "B83815", "48224A", "0FEE25", "38D33B", "C135B9",
    "71DE13", "7E8CFB", "B536C3", "B8E68C", "1E7272", "B75B69",
    "533F65", "812B41", "2C4DE1", "759E2C", "7E6434", "B43640",
    "05671B", "8C6FE0", "978D4E", "5319BF", "516499", "6E4292",
    "DC86A3", "8C46BE", "DF0B03", "E9A5B2", "B0E553", "07DBAC",
    "4778E4", "37AF4F",
]

@testset "Galileo E5b primary codes match the ICD first-24-chip tables" begin
    for (signal, table) in
        ((GalileoE5bI(), E5B_I_FIRST24_HEX), (GalileoE5bQ(), E5B_Q_FIRST24_HEX))
        codes = get_codes(signal)
        @test size(codes) == (10230, 50)
        for prn = 1:50
            @test codes[1:24, prn] == GNSSSignals.read_from_documentation(table[prn])
        end
    end
    # The two components are generated from different register-2 polynomials, so
    # they are different code sets even though register 1 is shared.
    @test get_codes(GalileoE5bI()) != get_codes(GalileoE5bQ())
    # ... and neither is the E5a code set (different register-1 polynomial too).
    @test get_codes(GalileoE5bI()) != get_codes(GalileoE5aI())
end

@testset "Galileo E5b primary codes match the ICD electronic annex" begin
    # The fixtures hold the full 10230-chip codes from the OS SIS ICD v2.2
    # electronic annex (C.5 / C.6), which lists the codes in hexadecimal rather
    # than as LFSR parameters. `GalileoE5bI`/`GalileoE5bQ` generate them from the
    # §3.4.1 LFSR definition instead, so this checks the whole code — including
    # the truncation at 10230 chips, which the first-24 table above cannot see.
    fixture_dir = joinpath(@__DIR__, "fixtures")
    for (comp, signal) in (("i", GalileoE5bI()), ("q", GalileoE5bQ()))
        for prn in (1, 25)
            ref = _load_packed_hex_fixture(
                joinpath(fixture_dir, "galileo_e5b_$(comp)_prn$(prn)_primary.hex.gz"),
                10230,
            )
            @test get_codes(signal)[:, prn] == ref
        end
    end
end

@testset "Galileo E5b-I secondary code (CS4)" begin
    e5b_i = GalileoE5bI()
    sec = get_secondary_code(e5b_i)
    # CS4 = 0xE -> 1110 (MSB-first), 0 -> -1, 1 -> +1.
    cs4 = Int8[1, 1, 1, -1]
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:3] == cs4
    # Shared across SVIDs, and wraps with period 4.
    @test [GNSSSignals.secondary_value(sec, 33, k) for k = 0:3] == cs4
    @test GNSSSignals.secondary_value(sec, 1, 4) == cs4[1]

    # One secondary chip per 1 ms primary period -> a 4 ms tiered code.
    prim = get_codes(e5b_i)[:, 1]
    @test get_code(e5b_i, 0, 1) == prim[1] * cs4[1]
    @test get_code(e5b_i, 3 * 10230, 1) == prim[1] * cs4[4]
    @test get_code(e5b_i, 4 * 10230, 1) == prim[1] * cs4[1]
end

@testset "Galileo E5b-Q secondary codes (CS100 51-100)" begin
    e5b_q = GalileoE5bQ()
    sec = get_secondary_code(e5b_q)
    @test size(sec.codes) == (100, 50)
    # E5b-Q is assigned CS100_(n+50), i.e. the half of the table E5a-Q does not
    # use — so no column may coincide with an E5a-Q column.
    e5a_q_codes = get_secondary_code(GalileoE5aQ()).codes
    @test all(sec.codes[:, prn] != e5a_q_codes[:, prn] for prn = 1:50)
    @test isempty(intersect(Set(eachcol(sec.codes)), Set(eachcol(e5a_q_codes))))
    # CS100_51 = 0xCFF914EE3C6126A49FD5E5C94 -> 1100 1111 1111 ... MSB-first.
    cs100_51 = Int8[1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:11] == cs100_51
    @test GNSSSignals.secondary_value(sec, 1, 100) == cs100_51[1]

    prim = get_codes(e5b_q)[:, 1]
    @test get_code(e5b_q, 0, 1) == prim[1] * cs100_51[1]
    @test get_code(e5b_q, 2 * 10230, 1) == prim[1] * cs100_51[3]
end

@testset "Galileo E5b code generation" begin
    # Chip-rate sampling reproduces the tiered code chip-for-chip: E5b-I with its
    # baked CS4 overlay, E5b-Q with the residual per-PRN CS100 one.
    for signal in (GalileoE5bI(), GalileoE5bQ())
        prn = 11
        n = 2 * 10230
        buffer = zeros(Int8, n)
        gen_code!(buffer, signal, prn, 10230e3Hz)
        @test buffer == [get_code(signal, i, prn) for i = 0:n-1]
    end
end
