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

@testset "TMBOC Int16 SIMD fast path matches two-pass fallback" begin
    # The TMBOC `multiply_with_subcarrier!` dispatch picks the 16-lane
    # SIMD.jl kernel (`_tmboc_simd_i16!`) for any contiguous Int16
    # buffer (`Vector`, `FastContiguousSubArray`, or other
    # `DenseVector{Int16}`), and the two-pass fallback
    # (`_tmboc_two_pass_i16!`) only for non-contiguous buffers or
    # patterns that don't fit in `UInt32`. Both kernels must agree
    # bit-for-bit. We exercise each directly so the assertion holds
    # even if the dispatcher is later relaxed further.
    sig = GPSL1C_P()
    modulation = get_modulation(sig)
    sampling_rate = 15e6Hz
    code_rate = 1023e3Hz

    # Replicate the dispatcher's parameter-prep so we can call the
    # private workers directly with matching inputs.
    function _prep(modulation, sampling_rate, code_rate, start_phase, start_index)
        _, sc_phase_boc1, sc_delta_boc1 = GNSSSignals.calc_subcarrier_phase_and_delta(
            modulation.boc1, sampling_rate, code_rate, start_phase, Int32,
        )
        _, sc_phase_boc2, sc_delta_boc2 = GNSSSignals.calc_subcarrier_phase_and_delta(
            modulation.boc2, sampling_rate, code_rate, start_phase, Int32,
        )
        pattern_bits64 = GNSSSignals._pack_tmboc_pattern(modulation.pattern)
        chip_delta_fp = ceil(
            UInt64, Float64(code_rate / sampling_rate) * (UInt64(1) << 32),
        )
        int_chip_offset = mod(floor(Int, start_phase), length(modulation.pattern))
        chip_acc_fp0 = UInt64(floor(UInt64, mod(start_phase, 1.0) * (UInt64(1) << 32)))
        si = reinterpret(UInt32, Int32(start_index))
        return (sc_phase_boc1, sc_delta_boc1, sc_phase_boc2, sc_delta_boc2,
                pattern_bits64, chip_delta_fp, int_chip_offset, chip_acc_fp0, si)
    end

    for (samples, start_phase, start_index) in (
        (2000, 0.0, 0),
        (2000, 1.7, 0),
        (2000, 0.0, 7),
        (2000, 13.5, -3),
        (97,   0.0, 0),  # exercises the scalar tail
        (256,  0.0, 0),  # exact multiple of 16
    )
        (sc1, sd1, sc2, sd2, pat, cdfp, ico, cafp0, si) = _prep(
            modulation, sampling_rate, code_rate, start_phase, start_index,
        )
        NPAT = length(modulation.pattern)

        buf_simd = ones(Int16, samples)
        GNSSSignals._tmboc_simd_i16!(
            buf_simd, sc1, sd1, sc2, sd2, UInt32(pat), cafp0, cdfp, ico, si, Val(NPAT),
        )

        buf_ref = ones(Int16, samples)
        GNSSSignals._tmboc_two_pass_i16!(
            buf_ref, sc1, sd1, sc2, sd2, pat, cafp0, cdfp, ico, si, Val(NPAT),
        )

        @test buf_simd == buf_ref

        # And: the dispatcher must route a contiguous `view` to the
        # SIMD path (regression guard — previously gated by
        # `isa Vector{Int16}`, which excluded views).
        buf_view_backing = ones(Int16, samples)
        v = view(buf_view_backing, 1:samples)
        GNSSSignals.multiply_with_subcarrier!(
            v, modulation, sampling_rate, code_rate, start_phase, start_index,
        )
        @test v == buf_simd
    end
end

@testset "L1C-D / L1C-P sample-stream matches external reference (PRN 1, 12.276 MHz)" begin
    # Independent cross-check against PocketSDR's L1C code generator
    # (primary code, overlay code, and TMBOC sub-carrier modulation).
    # At `sampling_frequency = 12 × code_frequency` the primary chip
    # boundaries and the BOC sub-carrier half-cycle boundaries both
    # land on sample boundaries, so the comparison is bit-exact.
    #
    # The fixtures hold one primary code period
    # (`12 × 10230 = 122760` samples) of ±1 values, packed 1 bit per
    # sample (least-significant-bit-first inside each hex nibble).
    #
    # The fixtures are reproducible from PocketSDR (`sdr_code.py`):
    #
    #   import sdr_code, numpy as np
    #   prn = 1
    #   # L1C-D: PocketSDR outputs 2 samples per chip; upsample by 6×.
    #   l1cd = np.repeat(np.asarray(sdr_code.gen_code_L1CD(prn)), 6)
    #   # L1C-P: PocketSDR is already at 12 samples per chip; apply
    #   # the first overlay bit.
    #   overlay0 = int(sdr_code.sec_code_L1CP(prn)[0])
    #   l1cp = np.asarray(sdr_code.gen_code_L1CP(prn)) * overlay0
    #
    # `l1cd` and `l1cp` reproduce the in-tree fixtures bit-for-bit
    # across all 122760 samples.
    fixture_dir = joinpath(@__DIR__, "fixtures")
    n_samples = 12 * 10230
    ref_d = _load_packed_hex_fixture(joinpath(fixture_dir, "l1c_d_prn1_fs12chip.hex.gz"), n_samples)
    ref_p = _load_packed_hex_fixture(joinpath(fixture_dir, "l1c_p_prn1_fs12chip.hex.gz"), n_samples)

    sampling_rate = 12.276e6Hz
    code_rate = 1023e3Hz

    buf_d = zeros(Int16, n_samples)
    gen_code!(buf_d, GPSL1C_D(), 1, sampling_rate, code_rate, 0.0, 0)
    @test buf_d == ref_d

    buf_p = zeros(Int16, n_samples)
    gen_code!(buf_p, GPSL1C_P(), 1, sampling_rate, code_rate, 0.0, 0)
    @test buf_p == ref_p
end

@testset "TMBOC `multiply_with_subcarrier!` rejects non-Int16 buffers" begin
    # Only `Int16` buffers are supported. Float32 / Float64 / Int32
    # callers get a `MethodError` instead of a silently slow path.
    sig = GPSL1C_P()
    modulation = get_modulation(sig)
    fs = 15e6Hz
    cf = 1023e3Hz

    for buf in (ones(Float32, 100), ones(Float64, 100), ones(Int32, 100))
        @test_throws MethodError GNSSSignals.multiply_with_subcarrier!(
            buf, modulation, fs, cf, 0.0, 0,
        )
    end
end
