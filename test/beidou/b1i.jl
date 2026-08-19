@testset "BeiDou B1I" begin
    b1i = BeiDouB1I()
    @test @inferred(get_band(b1i)) == B1I()
    @test @inferred(get_center_frequency(b1i)) == 1_561_098_000Hz
    @test @inferred(get_code_length(b1i)) == 2046
    @test @inferred(get_secondary_code_length(b1i)) == 20
    @test @inferred(get_secondary_code(b1i)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(b1i)) == 50Hz
    @test @inferred(get_code_frequency(b1i)) == 2_046_000Hz
    @test get_signal_name(b1i) == "BeiDou B1I"
    @test @inferred(get_modulation(b1i)) == GNSSSignals.LOC()
    @test get_code_type(b1i) === Int16
    @test GNSSSignals.get_code_factor(b1i) == 1

    # The B1I chip rate is 2× the C/A rate on the 1561.098 MHz carrier.
    @test get_code_center_frequency_ratio(b1i) ≈ 2_046_000 / 1_561_098_000
    @test get_code_spectrum(b1i, 0) ≈ 1.0Hz / get_code_frequency(b1i)
end

@testset "BeiDou B1I ranging codes are balanced Gold codes" begin
    # BDS-SIS-ICD-B1I-3.0 publishes no chip vectors, so the codes are checked
    # structurally: 63 distinct, ±1, length-2046 balanced Gold codes (the
    # length-2047 m-sequences truncated by one chip give |#(+1) − #(−1)| ≤ 2).
    codes = get_codes(BeiDouB1I())
    @test size(codes) == (2046, 63)
    _bd_check_pm1(codes, 2046)
    for prn = 1:63
        @test abs(sum(codes[:, prn])) <= 2
    end
    @test length(unique(eachcol(codes))) == 63
end

@testset "BeiDou B1I secondary code (NH20, D1 MEO/IGSO only)" begin
    b1i = BeiDouB1I()
    sec = get_secondary_code(b1i)
    # NH20 = 0 0 0 0 0 1 0 0 1 1 0 1 0 1 0 0 1 1 1 0, mapped 0 → -1, 1 → +1.
    nh20 = Int8[-1, -1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1]
    # MEO/IGSO PRNs (6-58, D1) carry NH20 and wrap with period 20.
    @test [GNSSSignals.secondary_value(sec, 7, k) for k = 0:19] == nh20
    @test [GNSSSignals.secondary_value(sec, 58, k) for k = 0:19] == nh20
    @test GNSSSignals.secondary_value(sec, 7, 20) == nh20[1]
    # GEO PRNs (1-5, 59-63, D2) carry no overlay — all-ones column (a no-op).
    for geo in (1, 5, 59, 63)
        @test all(GNSSSignals.secondary_value(sec, geo, k) == 1 for k = 0:19)
    end

    # Tiered code on a D1 PRN: period 0 uses NH20[0] = -1, period 1 uses NH20[1].
    prim7 = get_codes(b1i)[:, 7]
    @test get_code(b1i, 0, 7) == prim7[1] * nh20[1]
    @test get_code(b1i, 2046, 7) == prim7[1] * nh20[2]
    # A GEO PRN's tiered code equals its primary code (no overlay).
    prim1 = get_codes(b1i)[:, 1]
    @test get_code(b1i, 0, 1) == prim1[1]
    @test get_code(b1i, 2046, 1) == prim1[1]
end

@testset "BeiDou B1I gen_code! produces a sampled BPSK(2) replica" begin
    b1i = BeiDouB1I()
    samples_per_chip = 2
    code_freq = get_code_frequency(b1i)
    n_chips = 3 * get_code_length(b1i)           # exercise the NH20 overlay
    n_samples = samples_per_chip * n_chips
    prn = 6                                      # a D1 (MEO/IGSO) PRN, so NH20 applies
    buf = zeros(Int8, n_samples)
    gen_code!(buf, b1i, prn, samples_per_chip * code_freq, code_freq, 0.0, 0)
    @test all(x -> x == 1 || x == -1, buf)
    expected = [get_code(b1i, (k - 1) / samples_per_chip, prn) for k = 1:n_samples]
    @test buf == expected
end
