@testset "BeiDou B3I" begin
    b3i = BeiDouB3I()
    @test @inferred(get_band(b3i)) == B3I()
    @test @inferred(get_center_frequency(b3i)) == 1_268_520_000Hz
    @test @inferred(get_code_length(b3i)) == 10230
    @test @inferred(get_secondary_code_length(b3i)) == 20
    @test @inferred(get_secondary_code(b3i)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(b3i)) == 50Hz
    @test @inferred(get_code_frequency(b3i)) == 10_230_000Hz
    @test get_signal_name(b3i) == "BeiDou B3I"
    @test @inferred(get_modulation(b3i)) == GNSSSignals.LOC()
    @test get_code_type(b3i) === Int16
    @test GNSSSignals.get_code_factor(b3i) == 1
    @test get_code_spectrum(b3i, 0) ≈ 1.0Hz / get_code_frequency(b3i)
end

@testset "BeiDou B3I G2 initial phases match the ICD shift table" begin
    # BDS-SIS-ICD-B3I-1.0 Table 4-1 gives, per PRN, the G2 initial-phase state
    # and the number of shifts from all-ones that produces it. Re-running the G2
    # LFSR (using the package's own feedback taps) validates both the tabulated
    # initial-phase constants and the LFSR convention the generator relies on.
    taps = GNSSSignals.B3I_G2_FEEDBACK
    for (prn, nshift, init) in ICD_B3I_SHIFT
        state = fill(Int8(1), 13)
        for _ = 1:nshift
            fb = Int8(0)
            for t in taps
                fb ⊻= state[t]
            end
            for i = 13:-1:2
                state[i] = state[i-1]
            end
            state[1] = fb
        end
        want = Int8[c == '1' ? 1 : 0 for c in init]
        @test state == want
        @test GNSSSignals.B3I_G2_INITIAL_PHASES[prn] == init
    end
end

@testset "BeiDou B3I ranging codes are ±1 Gold codes" begin
    codes = get_codes(BeiDouB3I())
    @test size(codes) == (10230, 63)
    _bd_check_pm1(codes, 10230)
    @test length(unique(eachcol(codes))) == 63
end

@testset "BeiDou B3I secondary code (NH20, D1 MEO/IGSO only)" begin
    sec = get_secondary_code(BeiDouB3I())
    nh20 = Int8[-1, -1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1]
    @test [GNSSSignals.secondary_value(sec, 6, k) for k = 0:19] == nh20   # D1 (MEO/IGSO)
    for geo in (1, 5, 59, 63)                                             # D2 (GEO): no overlay
        @test all(GNSSSignals.secondary_value(sec, geo, k) == 1 for k = 0:19)
    end
end

@testset "BeiDou B3I gen_code! produces a sampled BPSK(10) replica" begin
    b3i = BeiDouB3I()
    samples_per_chip = 2
    code_freq = get_code_frequency(b3i)
    n_samples = samples_per_chip * 2 * get_code_length(b3i)
    prn = 6                                      # a D1 (MEO/IGSO) PRN, so NH20 applies
    buf = zeros(Int8, n_samples)
    gen_code!(buf, b3i, prn, samples_per_chip * code_freq, code_freq, 0.0, 0)
    @test all(x -> x == 1 || x == -1, buf)
    expected = [get_code(b3i, (k - 1) / samples_per_chip, prn) for k = 1:n_samples]
    @test buf == expected
end
