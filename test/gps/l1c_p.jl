@testset "GPS L1C-P" begin
    sig = GPSL1C_P()
    @test @inferred(get_band(sig)) == L1()
    @test @inferred(get_center_frequency(sig)) == 1.57542e9Hz
    @test @inferred(get_code_length(sig)) == 10230
    @test @inferred(get_secondary_code_length(sig)) == 1800
    @test @inferred(get_secondary_code(sig)) isa PerPRNSecondaryCode
    @test @inferred(get_bit_rate(sig)) == 0Hz   # pilot is dataless
    @test @inferred(get_chip_rate(sig)) == 1023e3Hz
    @test get_signal_name(sig) == "GPS L1C-P"
    @test get_replica_name(sig) == "GPS L1C-P (BOC(1,1) approximation)"
    @test @inferred(get_modulation(sig)) == GNSSSignals.BOCsin(1, 1)

    full = GPSL1C_P_TMBOC()
    @test get_signal_name(full) == "GPS L1C-P"
    @test get_replica_name(full) == "GPS L1C-P (TMBOC)"
    @test full.codes == sig.codes
    @test full.overlay_codes == sig.overlay_codes
    @test get_secondary_code_length(full) == get_secondary_code_length(sig)
    @test get_chip_rate(full) == get_chip_rate(sig)
    @test get_bit_rate(full) == get_bit_rate(sig)
    @test get_band(full) == get_band(sig)
    @test get_relative_power(full) == get_relative_power(sig)

    # Modulation: TMBOC over BOC(1,1) and BOC(6,1), pattern length 33.
    mod = @inferred get_modulation(full)
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

    # Check both subcarrier half-cycles and the overlay boundary independently
    # against the primary and overlay chips.
    for phase in (0.0, 0.25, 0.5, 1.0, 10230.0, 10230.5)
        chip = floor(Int, phase)
        primary = sig.codes[mod1(chip + 1, 10230), 1]
        overlay = sig.overlay_codes[mod1(fld(chip, 10230) + 1, 1800), 1]
        subcarrier = iseven(floor(Int, 2phase)) ? 1 : -1
        @test get_code(sig, phase, 1) == primary * overlay * subcarrier
    end
end

@testset "gen_code! integration: GPS L1C-D and L1C-P run without error" begin
    # Exercise the embedded LUT resample paths for both pilot replicas.
    sampling_rate = 25e6Hz
    samples = 4000

    let buf = zeros(Int8, samples)
        gen_code!(buf, GPSL1C_D(), 1, sampling_rate, 1023e3Hz, 0.0, 0)
        @test all(abs.(buf) .> 0)
    end

    for sig in (GPSL1C_P(), GPSL1C_P_TMBOC())
        buf = zeros(Int8, samples)
        gen_code!(buf, sig, 1, sampling_rate, 1023e3Hz, 0.0, 0)
        @test all(abs.(buf) .> 0)
        @test abs(sum(buf)) < samples * 0.5
    end

    sig = GPSL1C_P()
    buf = zeros(Int8, 1200)
    gen_code!(buf, sig, 1, 12.276e6Hz, 1023e3Hz, 0.0, 0)
    @test buf == [get_code(sig, (k - 1) / 12, 1) for k in eachindex(buf)]
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
    gen_code!(buf_p, GPSL1C_P_TMBOC(), 1, sampling_rate, code_rate, 0.0, 0)
    @test buf_p == ref_p
end
