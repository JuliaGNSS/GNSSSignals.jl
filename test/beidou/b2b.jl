@testset "BeiDou B2b-I" begin
    b2b = BeiDouB2bI()
    @test @inferred(get_band(b2b)) == B2b()
    @test @inferred(get_center_frequency(b2b)) == 1_207_140_000Hz
    @test @inferred(get_code_length(b2b)) == 10230
    @test @inferred(get_secondary_code_length(b2b)) == 1
    @test @inferred(get_secondary_code(b2b)) isa NoSecondaryCode
    @test @inferred(get_data_frequency(b2b)) == 1000Hz
    @test @inferred(get_code_frequency(b2b)) == 10_230_000Hz
    @test get_signal_name(b2b) == "BeiDou B2b-I"
    @test @inferred(get_modulation(b2b)) == GNSSSignals.LOC()
    @test get_code_type(b2b) === Int16
    @test get_code_spectrum(b2b, 0) ≈ 1.0Hz / get_code_frequency(b2b)
end

@testset "BeiDou B2b-I ranging codes match the ICD (PRN 6-58)" begin
    codes = get_codes(BeiDouB2bI())
    @test size(codes) == (10230, 63)
    _bd_verify_octal(codes, ICD_B2B_I)
    # Defined PRNs (6-58) are ±1; PRN indices outside that range are undefined
    # (all-zero), matching the 53-code B2b_I definition.
    for prn = 6:58
        @test all(x -> x == 1 || x == -1, codes[:, prn])
    end
    for prn in (1, 2, 3, 4, 5, 59, 60, 61, 62, 63)
        @test all(iszero, codes[:, prn])
    end
end

@testset "BeiDou B2b-I gen_code! produces a sampled BPSK(10) replica" begin
    b2b = BeiDouB2bI()
    samples_per_chip = 2
    code_freq = get_code_frequency(b2b)
    n_samples = samples_per_chip * 2 * get_code_length(b2b)
    buf = zeros(Int8, n_samples)
    gen_code!(buf, b2b, 6, samples_per_chip * code_freq, code_freq, 0.0, 0)
    @test all(x -> x == 1 || x == -1, buf)
    expected = [get_code(b2b, (k - 1) / samples_per_chip, 6) for k = 1:n_samples]
    @test buf == expected
end
