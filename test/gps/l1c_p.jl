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

@testset "TMBOC SIMD fast path matches two-pass fallback" begin
    # The Int16/Float32 fast paths use an explicit SIMD.jl kernel; the
    # generic path falls back to a two-pass scalar/auto-vectorized
    # implementation. They must agree bit-for-bit on Int16 and within
    # floating-point rounding on Float32.
    sig = GPSL1C_P()
    modulation = get_modulation(sig)
    sampling_rate = 15e6Hz
    code_rate = 1023e3Hz

    for (samples, start_phase, start_index) in (
        (2000, 0.0, 0),
        (2000, 1.7, 0),
        (2000, 0.0, 7),
        (2000, 13.5, -3),
        (97,   0.0, 0),  # exercises the scalar tail
        (256,  0.0, 0),  # exact multiple of 16
    )
        buf_simd_i16 = ones(Int16, samples)
        buf_ref_i16  = ones(Int16, samples)
        GNSSSignals.multiply_with_subcarrier!(
            buf_simd_i16, modulation, sampling_rate, code_rate, start_phase, start_index,
        )
        # Force the generic two-pass fallback by passing a typed view.
        v_ref = view(buf_ref_i16, 1:samples)
        GNSSSignals.multiply_with_subcarrier!(
            v_ref, modulation, sampling_rate, code_rate, start_phase, start_index,
        )
        @test buf_simd_i16 == buf_ref_i16

        buf_simd_f32 = ones(Float32, samples)
        buf_ref_f32  = ones(Float32, samples)
        GNSSSignals.multiply_with_subcarrier!(
            buf_simd_f32, modulation, sampling_rate, code_rate, start_phase, start_index,
        )
        v_ref_f32 = view(buf_ref_f32, 1:samples)
        GNSSSignals.multiply_with_subcarrier!(
            v_ref_f32, modulation, sampling_rate, code_rate, start_phase, start_index,
        )
        @test buf_simd_f32 == buf_ref_f32
    end
end
