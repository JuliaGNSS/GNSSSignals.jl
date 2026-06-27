@testset "Galileo E5a-I" begin
    e5a_i = GalileoE5aI()
    @test @inferred(get_band(e5a_i)) == L5()
    @test @inferred(get_center_frequency(e5a_i)) == 1.17645e9Hz
    @test @inferred(get_code_length(e5a_i)) == 10230
    @test @inferred(get_secondary_code_length(e5a_i)) == 20
    @test @inferred(get_secondary_code(e5a_i)) isa SharedSecondaryCode{20}
    @test @inferred(get_data_frequency(e5a_i)) == 50Hz
    @test @inferred(get_code_frequency(e5a_i)) == 10230e3Hz
    @test get_signal_name(e5a_i) == "Galileo E5a-I"
    @test @inferred(get_modulation(e5a_i)) == GNSSSignals.LOC()
    @test get_code_type(e5a_i) === Int16

    # PRN 1 starts on a -1 chip; CS20[0] = +1, so the tiered chip equals the
    # primary chip here (and the fixture test below pins the full code).
    @test @inferred(get_code(e5a_i, 0, 1)) == -1
    @test @inferred(get_code(e5a_i, 0.0, 1)) == -1
    @test @inferred(GNSSSignals.get_code_unsafe(e5a_i, 0.0, 1)) == -1

    # BPSK(10) spectrum and code/center-frequency ratio (E5a shares the L5
    # band, so the ratio is the same 1/115 as GPS L5-I).
    @test GNSSSignals.get_code_factor(e5a_i) == 1
    @test get_code_spectrum(e5a_i, 0) ≈ 1.0Hz / get_code_frequency(e5a_i)
    @test get_code_center_frequency_ratio(e5a_i) ≈ 1 / 115
end

@testset "Galileo E5a-Q" begin
    e5a_q = GalileoE5aQ()
    @test @inferred(get_band(e5a_q)) == L5()
    @test @inferred(get_center_frequency(e5a_q)) == 1.17645e9Hz
    @test @inferred(get_code_length(e5a_q)) == 10230
    @test @inferred(get_secondary_code_length(e5a_q)) == 100
    @test @inferred(get_secondary_code(e5a_q)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(e5a_q)) == 0Hz   # dataless pilot
    @test @inferred(get_code_frequency(e5a_q)) == 10230e3Hz
    @test get_signal_name(e5a_q) == "Galileo E5a-Q"
    @test @inferred(get_modulation(e5a_q)) == GNSSSignals.LOC()
    @test get_code_type(e5a_q) === Int16

    @test GNSSSignals.get_code_factor(e5a_q) == 1
    @test get_code_spectrum(e5a_q, 0) ≈ 1.0Hz / get_code_frequency(e5a_q)
    @test get_code_center_frequency_ratio(e5a_q) ≈ 1 / 115
end

@testset "Galileo E5a primary codes match the ICD/GNSS-SDR reference" begin
    # The reference fixtures hold the GNSS-SDR memory codes (which equal the
    # Galileo OS SIS ICD v2.2 E5a primary codes) decoded MSB-first into the
    # package's chip convention. Our `GalileoE5aI`/`GalileoE5aQ` generate the
    # codes independently from the ICD §3.5 LFSR definition, so this is a
    # cross-source check, not a self-comparison.
    fixture_dir = joinpath(@__DIR__, "fixtures")
    e5a_i = GalileoE5aI()
    e5a_q = GalileoE5aQ()

    for (comp, signal) in (("i", e5a_i), ("q", e5a_q))
        for prn in (1, 25)
            ref = _load_packed_hex_fixture(
                joinpath(fixture_dir, "galileo_e5a_$(comp)_prn$(prn)_primary.hex.gz"),
                10230,
            )
            @test get_codes(signal)[:, prn] == ref
        end
    end
end

@testset "Galileo E5a-I secondary code (CS20)" begin
    e5a_i = GalileoE5aI()
    sec = get_secondary_code(e5a_i)
    # CS20 = 0x842E9 -> 1000 0100 0010 1110 1001 (MSB-first), 0 -> -1, 1 -> +1.
    cs20 = Int8[1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1]
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:19] == cs20
    # Shared across SVIDs, and wraps with period 20.
    @test [GNSSSignals.secondary_value(sec, 7, k) for k = 0:19] == cs20
    @test GNSSSignals.secondary_value(sec, 1, 20) == cs20[1]

    # The tiered code applies the secondary chip once per 1 ms primary period:
    # period 0 uses CS20[0] = +1, period 1 uses CS20[1] = -1.
    prim = get_codes(e5a_i)[:, 1]
    @test get_code(e5a_i, 0, 1) == prim[1] * cs20[1]
    @test get_code(e5a_i, 10230, 1) == prim[1] * cs20[2]
end

@testset "Galileo E5a-Q secondary code (CS100)" begin
    e5a_q = GalileoE5aQ()
    sec = get_secondary_code(e5a_q)
    @test sec isa PerPRNSecondaryCode

    # CS100_1 = 0x83F6F69D... -> first nibble 8 = 1000 -> +1,-1,-1,-1.
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:3] == Int8[1, -1, -1, -1]
    # Per-PRN: PRN 2 (CS100_2 = 0x6655...) first nibble 6 = 0110 -> -1,+1,+1,-1.
    @test [GNSSSignals.secondary_value(sec, 2, k) for k = 0:3] == Int8[-1, 1, 1, -1]
    # Wraps with period 100.
    @test GNSSSignals.secondary_value(sec, 1, 100) == GNSSSignals.secondary_value(sec, 1, 0)

    # Tiered code: period 1 multiplies the primary by CS100_1[1] = -1.
    prim = get_codes(e5a_q)[:, 1]
    @test get_code(e5a_q, 10230, 1) == prim[1] * GNSSSignals.secondary_value(sec, 1, 1)
end

@testset "Galileo E5a gen_code! produces a sampled BPSK(10) replica" begin
    # gen_code!'s fixed-point sampler is already covered for LOC signals by the
    # GPS L5-I and Galileo E1B(BOC11) tests; here we confirm the E5a wiring
    # (code table + secondary selection) drives it correctly. At an integer
    # number of samples per chip, every sample equals the corresponding
    # `get_code` chip, so the sampled replica is the tiered code, upsampled.
    for signal in (GalileoE5aI(), GalileoE5aQ())
        samples_per_chip = 2
        code_freq = get_code_frequency(signal)
        # Cover several primary periods so the secondary overlay is exercised.
        n_chips = 3 * get_code_length(signal)
        n_samples = samples_per_chip * n_chips
        buf = zeros(Int16, n_samples)
        gen_code!(buf, signal, 1, samples_per_chip * code_freq, code_freq, 0.0, 0)

        @test all(x -> x == 1 || x == -1, buf)
        expected = [get_code(signal, (k - 1) / samples_per_chip, 1) for k = 1:n_samples]
        @test buf == expected
    end
end
