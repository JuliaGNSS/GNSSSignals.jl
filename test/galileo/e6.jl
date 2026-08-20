@testset "Galileo E6-B" begin
    e6b = GalileoE6B()
    @test @inferred(get_band(e6b)) == E6()
    @test @inferred(get_center_frequency(e6b)) == 1.27875e9Hz
    @test @inferred(get_code_length(e6b)) == 5115
    @test @inferred(get_secondary_code_length(e6b)) == 1
    @test @inferred(get_secondary_code(e6b)) isa NoSecondaryCode
    @test @inferred(get_data_frequency(e6b)) == 1000Hz     # C/NAV symbol rate
    @test @inferred(get_code_frequency(e6b)) == 5115e3Hz
    @test get_signal_name(e6b) == "Galileo E6-B"
    @test @inferred(get_modulation(e6b)) == GNSSSignals.LOC()
    @test get_code_type(e6b) === Int16

    # PRN 1's first hex symbol is E6 -> 1110, so the first chip is +1 (the
    # untiered primary chip: E6-B carries no secondary code).
    @test @inferred(get_code(e6b, 0, 1)) == 1
    @test @inferred(get_code(e6b, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(e6b, 0.0, 1)) == 1
    # 1 ms code period: chip 5115 wraps to chip 0.
    @test get_code(e6b, 5115, 1) == get_code(e6b, 0, 1)

    # BPSK(5) spectrum, and 5.115 MHz against the 1278.75 MHz E6 carrier.
    @test GNSSSignals.get_code_factor(e6b) == 1
    @test get_code_spectrum(e6b, 0) ≈ 1.0Hz / get_code_frequency(e6b)
    @test get_code_center_frequency_ratio(e6b) ≈ 1 / 250
end

@testset "Galileo E6-C" begin
    e6c = GalileoE6C()
    @test @inferred(get_band(e6c)) == E6()
    @test @inferred(get_center_frequency(e6c)) == 1.27875e9Hz
    @test @inferred(get_code_length(e6c)) == 5115
    @test @inferred(get_secondary_code_length(e6c)) == 100
    @test @inferred(get_secondary_code(e6c)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(e6c)) == 0Hz        # dataless pilot
    @test @inferred(get_code_frequency(e6c)) == 5115e3Hz
    @test get_signal_name(e6c) == "Galileo E6-C"
    @test @inferred(get_modulation(e6c)) == GNSSSignals.LOC()
    @test get_code_type(e6c) === Int16

    @test GNSSSignals.get_code_factor(e6c) == 1
    @test get_code_spectrum(e6c, 0) ≈ 1.0Hz / get_code_frequency(e6c)
    @test get_code_center_frequency_ratio(e6c) ≈ 1 / 250
end

@testset "Galileo E6 primary codes match the Technical Note listing" begin
    # The E6-B/C primary codes are memory codes, so there is no generator to
    # check them against: the check that matters is that the shipped binaries
    # decode to the hexadecimal listing published in the Galileo E6-B/C Codes
    # Technical Note (Issue 1, January 2019), Annex A.3 / A.4. The literals
    # below are transcribed from that listing — the first ten and last four hex
    # symbols of three codes per component. The same listing was verified
    # chip-for-chip against GNSS-SDR (`Galileo_E6.h`) and PocketSDR
    # (`sdr_code_gal.py`) for all 50 codes of both components while generating
    # `data/codes_galileo_e6*.bin`.
    #
    # 5115 is not divisible by four, so the ICD pads the last hex symbol with a
    # single zero "at the end in time" (Technical Note A.2): the last four
    # symbols carry chips 5101-5115 plus that pad, which must be 0.
    head_tail = Dict(
        (:B, 1) => ("E6648AA5EF", "338C"),
        (:B, 25) => ("A38DA6E20C", "33E8"),
        (:B, 50) => ("FE62E6E553", "57E4"),
        (:C, 1) => ("F5A3D656F9", "784C"),
        (:C, 25) => ("A1F6FA8D82", "5430"),
        (:C, 50) => ("CE11D9F5C8", "F334"),
    )
    for (component, signal) in ((:B, GalileoE6B()), (:C, GalileoE6C()))
        codes = get_codes(signal)
        @test size(codes) == (5115, 50)
        for prn in (1, 25, 50)
            head_hex, tail_hex = head_tail[(component, prn)]
            head = GNSSSignals.read_from_documentation(head_hex)
            tail = GNSSSignals.read_from_documentation(tail_hex)
            @test codes[1:40, prn] == head
            @test codes[5101:5115, prn] == tail[1:15]
            @test tail[16] == -1     # the padded zero chip, decoded as -1
        end
        # Every code is balanced to one chip (5115 is odd) and distinct from the
        # others, which pins the PRN ordering of the binary as well.
        for prn = 1:50
            @test sum(codes[:, prn]) == -1
        end
        @test length(unique(eachcol(codes))) == 50
    end
    # The two components' codes are different code sets, not one shared table.
    @test get_codes(GalileoE6B()) != get_codes(GalileoE6C())
end

@testset "Galileo E6-C secondary code (CS100)" begin
    e6c = GalileoE6C()
    sec = get_secondary_code(e6c)
    # E6-C is assigned the same CS100_1-50 as E5a-Q (Technical Note §2.4).
    @test sec.codes == get_secondary_code(GalileoE5aQ()).codes
    # CS100_1 = 0x83F6F69D8F6E15411FB8C9B1C -> 1000 0011 1111 ... MSB-first.
    cs100_1 = Int8[1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1]
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:11] == cs100_1
    # Per-SVID, and wraps with period 100.
    @test [GNSSSignals.secondary_value(sec, 2, k) for k = 0:11] != cs100_1
    @test GNSSSignals.secondary_value(sec, 1, 100) == cs100_1[1]

    # The tiered code applies one secondary chip per 1 ms primary period, so the
    # 100-chip overlay gives a 100 ms tiered code.
    prim = get_codes(e6c)[:, 1]
    @test get_code(e6c, 0, 1) == prim[1] * cs100_1[1]
    @test get_code(e6c, 5115, 1) == prim[1] * cs100_1[2]
    @test get_code(e6c, 100 * 5115, 1) == prim[1] * cs100_1[1]
end

@testset "Galileo E6 code generation" begin
    # Chip-rate sampling reproduces the tiered code chip-for-chip for both
    # components (E6-B untiered, E6-C with its CS100 overlay applied).
    for signal in (GalileoE6B(), GalileoE6C())
        prn = 7
        n = 3 * 5115
        buffer = zeros(Int8, n)
        gen_code!(buffer, signal, prn, 5115e3Hz)
        @test buffer == [get_code(signal, i, prn) for i = 0:n-1]
    end
    # 2 samples per chip at the E6 chip rate: each chip appears twice.
    e6b = GalileoE6B()
    buffer = zeros(Int8, 20)
    gen_code!(buffer, e6b, 1, 2 * 5115e3Hz)
    @test buffer == repeat([get_code(e6b, i, 1) for i = 0:9]; inner = 2)
end
