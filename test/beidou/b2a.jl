@testset "BeiDou B2a data (I)" begin
    b2a_i = BeiDouB2aI()
    @test @inferred(get_band(b2a_i)) == L5()
    @test @inferred(get_center_frequency(b2a_i)) == 1_176_450_000Hz
    @test @inferred(get_code_length(b2a_i)) == 10230
    @test @inferred(get_secondary_code_length(b2a_i)) == 5
    @test @inferred(get_secondary_code(b2a_i)) isa SharedSecondaryCode{5}
    @test @inferred(get_data_frequency(b2a_i)) == 200Hz
    @test @inferred(get_code_frequency(b2a_i)) == 10_230_000Hz
    @test get_signal_name(b2a_i) == "BeiDou B2a data"
    @test @inferred(get_modulation(b2a_i)) == GNSSSignals.LOC()
    @test get_code_type(b2a_i) === Int16
    # Shares the L5 band, so the code/center ratio is the same 1/115 as GPS L5.
    @test get_code_center_frequency_ratio(b2a_i) ≈ 1 / 115
    @test get_code_spectrum(b2a_i, 0) ≈ 1.0Hz / get_code_frequency(b2a_i)
end

@testset "BeiDou B2a pilot (Q)" begin
    b2a_q = BeiDouB2aQ()
    @test @inferred(get_band(b2a_q)) == L5()
    @test @inferred(get_center_frequency(b2a_q)) == 1_176_450_000Hz
    @test @inferred(get_code_length(b2a_q)) == 10230
    @test @inferred(get_secondary_code_length(b2a_q)) == 100
    @test @inferred(get_secondary_code(b2a_q)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(b2a_q)) == 0Hz   # dataless pilot
    @test @inferred(get_code_frequency(b2a_q)) == 10_230_000Hz
    @test get_signal_name(b2a_q) == "BeiDou B2a pilot"
    @test @inferred(get_modulation(b2a_q)) == GNSSSignals.LOC()
end

@testset "BeiDou B2a primary codes match the ICD" begin
    codes_i = get_codes(BeiDouB2aI())
    codes_q = get_codes(BeiDouB2aQ())
    @test size(codes_i) == (10230, 63)
    @test size(codes_q) == (10230, 63)
    _bd_verify_octal(codes_i, ICD_B2A_DATA)
    _bd_verify_octal(codes_q, ICD_B2A_PILOT)
end

@testset "BeiDou B2a data secondary code (00010)" begin
    b2a_i = BeiDouB2aI()
    sec = get_secondary_code(b2a_i)
    # Fixed 5-bit 00010, shared across SVIDs, mapped 0 → -1, 1 → +1.
    s5 = Int8[-1, -1, -1, 1, -1]
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:4] == s5
    @test [GNSSSignals.secondary_value(sec, 30, k) for k = 0:4] == s5
    @test GNSSSignals.secondary_value(sec, 1, 5) == s5[1]
end

@testset "BeiDou B2a pilot secondary codes match the ICD" begin
    b2a_q = BeiDouB2aQ()
    sec_mat = b2a_q.secondary_codes
    @test size(sec_mat) == (100, 63)
    _bd_verify_octal(sec_mat, ICD_B2A_PILOT_SEC)
    # And via the public secondary interface (period-100 per-PRN overlay).
    sec = get_secondary_code(b2a_q)
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:99] == sec_mat[:, 1]
    @test GNSSSignals.secondary_value(sec, 1, 100) == sec_mat[1, 1]
end

@testset "BeiDou B2a gen_code! produces a sampled BPSK(10) replica" begin
    for signal in (BeiDouB2aI(), BeiDouB2aQ())
        samples_per_chip = 2
        code_freq = get_code_frequency(signal)
        n_samples = samples_per_chip * 2 * get_code_length(signal)
        buf = zeros(Int8, n_samples)
        gen_code!(buf, signal, 1, samples_per_chip * code_freq, code_freq, 0.0, 0)
        @test all(x -> x == 1 || x == -1, buf)
        expected = [get_code(signal, (k - 1) / samples_per_chip, 1) for k = 1:n_samples]
        @test buf == expected
    end
end
