@testset "GPS L1C-P" begin
    sig = GPSL1C_P()
    @test @inferred(get_band(sig)) == L1()
    @test @inferred(get_center_frequency(sig)) == 1.57542e9Hz
    @test @inferred(get_code_length(sig)) == 10230
    @test @inferred(get_secondary_code_length(sig)) == 1800
    @test @inferred(get_secondary_code(sig)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(sig)) == 0Hz   # pilot is dataless
    @test @inferred(get_code_frequency(sig)) == 1023e3Hz
    @test get_signal_name(sig) == "GPS L1C-P"

    # Modulation: TMBOC over BOC(1,1) and BOC(6,1), pattern length 33.
    mod = @inferred get_modulation(sig)
    @test mod isa GNSSSignals.TMBOC
    @test mod.boc1 == GNSSSignals.BOCsin(1, 1)
    @test mod.boc2 == GNSSSignals.BOCsin(6, 1)
    @test length(mod.pattern) == 33
    # 4 of 33 positions are BOC(6,1) per IS-GPS-800G §3.3.
    @test count(mod.pattern) == 4
    # Specifically positions 0, 4, 6, 29 (zero-based) use BOC(6,1).
    @test mod.pattern[1]  == true   # position 0
    @test mod.pattern[5]  == true   # position 4
    @test mod.pattern[7]  == true   # position 6
    @test mod.pattern[30] == true   # position 29
    @test mod.pattern[2]  == false  # position 1
    @test mod.pattern[33] == false  # position 32

    # Negated codes should be the element-wise negation.
    @test sig.negated_codes == .-sig.codes

    # `_select_codes_for` returns the positive matrix for sec_val > 0 and
    # negated for sec_val < 0; multiplier is `true` in both cases so the
    # inner-loop multiply elides.
    let
        codes_pos, mul_pos = GNSSSignals._select_codes_for(sig, Int8(1))
        codes_neg, mul_neg = GNSSSignals._select_codes_for(sig, Int8(-1))
        @test codes_pos === sig.codes
        @test codes_neg === sig.negated_codes
        @test mul_pos === true
        @test mul_neg === true
    end

    @test GNSSSignals.get_code_factor(sig) == 1
end

@testset "gen_code! integration: GPS L1C-D and L1C-P run without error" begin
    # The TMBOC subcarrier multiply path is exercised only here; make
    # sure it produces a finite, non-zero output of the right type.
    sampling_rate = 25e6Hz
    samples = 4000

    let buf = zeros(get_code_type(GPSL1C_D()), samples)
        gen_code!(buf, GPSL1C_D(), 1, sampling_rate, 1023e3Hz, 0.0, 0)
        @test all(abs.(buf) .> 0)
    end

    let buf = zeros(get_code_type(GPSL1C_P()), samples)
        gen_code!(buf, GPSL1C_P(), 1, sampling_rate, 1023e3Hz, 0.0, 0)
        @test all(abs.(buf) .> 0)
        # The first 33 chips of every primary period mix BOC(1,1) and
        # BOC(6,1); over thousands of samples the average should be
        # close to zero (BOC subcarrier is zero-mean).
        @test abs(sum(buf)) < samples * 0.5
    end
end
