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

    @test GNSSSignals.get_code_factor(sig) == 1
end

@testset "gen_code! integration: GPS L1C-D and L1C-P run without error" begin
    # The embedded LUT TMBOC bake/resample path is exercised here; make sure it produces a
    # finite, non-zero ±1 Int8 output.
    sampling_rate = 25e6Hz
    samples = 4000

    let buf = zeros(Int8, samples)
        gen_code!(buf, GPSL1C_D(), 1, sampling_rate, 1023e3Hz, 0.0, 0)
        @test all(abs.(buf) .> 0)
    end

    let buf = zeros(Int8, samples)
        gen_code!(buf, GPSL1C_P(), 1, sampling_rate, 1023e3Hz, 0.0, 0)
        @test all(abs.(buf) .> 0)
        # The first 33 chips of every primary period mix BOC(1,1) and
        # BOC(6,1); over thousands of samples the average should be
        # close to zero (BOC subcarrier is zero-mean).
        @test abs(sum(buf)) < samples * 0.5
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

    # The embedded LUT emits Int8 ±1; at this chip-aligned rate it is byte-exact to the
    # spec-aligned ±1 fixtures (permute regime, integer samples/chip).
    buf_d = zeros(Int8, n_samples)
    gen_code!(buf_d, GPSL1C_D(), 1, sampling_rate, code_rate, 0.0, 0)
    @test buf_d == ref_d

    buf_p = zeros(Int8, n_samples)
    gen_code!(buf_p, GPSL1C_P(), 1, sampling_rate, code_rate, 0.0, 0)
    @test buf_p == ref_p
end
