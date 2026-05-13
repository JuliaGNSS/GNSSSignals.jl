@testset "Modulation" begin
    galileo_e1b = GalileoE1B()
    boc_sin = GNSSSignals.BOCsin(1, 1)
    @test GNSSSignals.get_subcarrier_code(boc_sin, 0.0) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 0.25) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 0.75) == -1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 1.25) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 1.75) == -1.0

    boc_cos = GNSSSignals.BOCcos(1, 1)
    @test GNSSSignals.get_subcarrier_code(boc_cos, 0.0) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 0.5) == -1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 1.0) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 1.5) == -1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 2.0) == 1.0

    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCcos(0, 1)
    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCcos(1, 0)
    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCsin(0, 1)
    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCsin(1, 0)

    @test_throws ErrorException(
        "Power of BOC1 must be between 0 and 1 and n of both BOCs must match",
    ) GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), -1.0)
    @test_throws ErrorException(
        "Power of BOC1 must be between 0 and 1 and n of both BOCs must match",
    ) GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), 2.0)
    @test_throws ErrorException(
        "Power of BOC1 must be between 0 and 1 and n of both BOCs must match",
    ) GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 2), 1 / 11)

    cboc = GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), 1 / 2)
    @test GNSSSignals.get_subcarrier_code(cboc, 0.0) ≈ sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 0.25) ≈ sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 0.75) ≈ -sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 1.25) ≈ sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 1.75) ≈ -sqrt(1 / 2) * 2

    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCcos(2, 1), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCcos(2, 2), 2.3) == 4
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCsin(2, 1), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCsin(2, 2), 2.3) == 4

    @test GNSSSignals.get_code_factor(GNSSSignals.BOCcos(2, 1)) == 1
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCcos(2, 2)) == 2
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCcos(2, 2.5)) == 2.5
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCsin(2, 1)) == 1
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCsin(2, 2)) == 2
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCsin(2, 2.5)) == 2.5
    @test GNSSSignals.get_code_factor(
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), 1 / 2),
    ) == 1
    @test GNSSSignals.get_code_factor(
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 2), GNSSSignals.BOCsin(1, 2), 1 / 2),
    ) == 2
    @test GNSSSignals.get_code_factor(
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 2.5), GNSSSignals.BOCsin(1, 2.5), 1 / 2),
    ) == 2.5

    # Test LOC floored phase
    @test GNSSSignals.get_floored_phase(GNSSSignals.LOC(), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.LOC(), 5.9) == 5

    # Test CBOC floored phase
    cboc_test = GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 2), GNSSSignals.BOCsin(6, 2), 10 / 11)
    @test GNSSSignals.get_floored_phase(cboc_test, 2.3) == 4

    # Test LOC code factor
    @test GNSSSignals.get_code_factor(GNSSSignals.LOC()) == 1
end

@testset "TMBOC accessors" begin
    # The L1C-P TMBOC: BOC(6,1) at chip-positions {0, 4, 6, 29}, BOC(1,1) elsewhere.
    boc1 = GNSSSignals.BOCsin(1, 1)
    boc6 = GNSSSignals.BOCsin(6, 1)
    pattern = ntuple(k -> (k - 1) ∈ (0, 4, 6, 29), Val(33))
    tmboc = GNSSSignals.TMBOC(boc1, boc6, pattern)

    # `get_subcarrier_code`: at each chip position pick the BOC variant
    # selected by the pattern, then evaluate that BOC's subcarrier at
    # the given phase.
    @testset "get_subcarrier_code" begin
        # Chip-position 0 is BOC(6,1): subcarrier flips every 1/12 chip.
        # First half-cycle (phase ∈ [0, 1/12)) is +1.
        @test GNSSSignals.get_subcarrier_code(tmboc, 0.0) == 1
        @test GNSSSignals.get_subcarrier_code(tmboc, 1 / 12 - 0.001) == 1
        # Second half-cycle (phase ∈ [1/12, 2/12)) is −1.
        @test GNSSSignals.get_subcarrier_code(tmboc, 1 / 12 + 0.001) == -1
        @test GNSSSignals.get_subcarrier_code(tmboc, 2 / 12 - 0.001) == -1
        # Third half-cycle returns to +1.
        @test GNSSSignals.get_subcarrier_code(tmboc, 2 / 12 + 0.001) == 1

        # Chip-position 1 is BOC(1,1): subcarrier flips every 1/2 chip.
        @test GNSSSignals.get_subcarrier_code(tmboc, 1.1) == 1   # first half +1
        @test GNSSSignals.get_subcarrier_code(tmboc, 1.6) == -1  # second half −1

        # Chip-position 4 is BOC(6,1) again.
        @test GNSSSignals.get_subcarrier_code(tmboc, 4.0) == 1
        @test GNSSSignals.get_subcarrier_code(tmboc, 4 + 1 / 12 + 0.001) == -1

        # The pattern repeats every 33 chips.
        @test GNSSSignals.get_subcarrier_code(tmboc, 33.0) ==
              GNSSSignals.get_subcarrier_code(tmboc, 0.0)
        @test GNSSSignals.get_subcarrier_code(tmboc, 33 + 1.1) ==
              GNSSSignals.get_subcarrier_code(tmboc, 1.1)
    end

    # `get_floored_phase`: integer chip index from a continuous phase.
    # TMBOC uses the boc1 `n` multiplier (same as CBOC).
    @testset "get_floored_phase" begin
        @test GNSSSignals.get_floored_phase(tmboc, 0.0) == 0
        @test GNSSSignals.get_floored_phase(tmboc, 1.7) == 1
        @test GNSSSignals.get_floored_phase(tmboc, 32.9) == 32
        @test GNSSSignals.get_floored_phase(tmboc, 33.0) == 33

        # With a non-unit `n`, the result scales by `n`.
        boc1_n2 = GNSSSignals.BOCsin(1, 2)
        boc6_n2 = GNSSSignals.BOCsin(6, 2)
        tmboc_n2 = GNSSSignals.TMBOC(boc1_n2, boc6_n2, pattern)
        @test GNSSSignals.get_floored_phase(tmboc_n2, 2.3) == 4
    end

    # `get_code_spectrum`: power-weighted mix of the BOC component
    # spectra, where the weights are the fractions of chips using each.
    @testset "get_code_spectrum" begin
        sig = GPSL1C_P()

        # Sanity: the PSD is real and non-negative at sampled frequencies.
        for f in (-2e6, -100e3, 100e3, 2e6)
            @test get_code_spectrum(sig, f) >= 0
        end

        # PSD value matches the explicit weighted-sum formula at a few
        # representative frequencies. 4 of 33 chips use BOC(6,1) →
        # BOC(6,1) weight is 4/33, BOC(1,1) weight is 29/33.
        n_boc2 = 4
        N = 33
        boc2_frac = n_boc2 / N
        for f in (50e3, 250e3, 1.5e6)
            expected =
                GNSSSignals.get_code_spectrum(boc1, sig, f) * (1 - boc2_frac) +
                GNSSSignals.get_code_spectrum(boc6, sig, f) * boc2_frac
            @test get_code_spectrum(sig, f) ≈ expected
        end

        # Symmetry: BOC subcarriers are zero-mean, so the spectrum
        # vanishes at DC.
        @test get_code_spectrum(sig, 0) == 0.0

        # Spectrum is even.
        @test get_code_spectrum(sig, 500e3) ≈ get_code_spectrum(sig, -500e3)
    end
end

@testset "Spectrum functions with units" begin
    # Test BPSK spectrum with different unit combinations
    fc = 1.023e6
    f = 500e3

    # All numbers
    psd1 = GNSSSignals.get_code_spectrum_BPSK(fc, f)
    @test psd1 > 0

    # fc with units, f number
    psd2 = GNSSSignals.get_code_spectrum_BPSK(fc * 1Hz, f)
    @test psd2 ≈ psd1

    # fc number, f with units
    psd3 = GNSSSignals.get_code_spectrum_BPSK(fc, f * 1Hz)
    @test psd3 ≈ psd1

    # Both with units
    psd4 = GNSSSignals.get_code_spectrum_BPSK(fc * 1Hz, f * 1Hz)
    @test psd4 ≈ psd1
end

@testset "BOC spectrum functions with units" begin
    fc = 1.023e6
    fs = 1.023e6
    f = 500e3

    # BOCsin spectrum
    psd_sin1 = GNSSSignals.get_code_spectrum_BOCsin(fc, fs, f)
    @test psd_sin1 >= 0

    psd_sin2 = GNSSSignals.get_code_spectrum_BOCsin(fc * 1Hz, fs * 1Hz, f)
    @test psd_sin2 ≈ psd_sin1

    psd_sin3 = GNSSSignals.get_code_spectrum_BOCsin(fc, fs, f * 1Hz)
    @test psd_sin3 ≈ psd_sin1

    psd_sin4 = GNSSSignals.get_code_spectrum_BOCsin(fc * 1Hz, fs * 1Hz, f * 1Hz)
    @test psd_sin4 ≈ psd_sin1

    # BOCcos spectrum
    psd_cos1 = GNSSSignals.get_code_spectrum_BOCcos(fc, fs, f)
    @test psd_cos1 >= 0

    psd_cos2 = GNSSSignals.get_code_spectrum_BOCcos(fc * 1Hz, fs * 1Hz, f)
    @test psd_cos2 ≈ psd_cos1

    psd_cos3 = GNSSSignals.get_code_spectrum_BOCcos(fc, fs, f * 1Hz)
    @test psd_cos3 ≈ psd_cos1

    psd_cos4 = GNSSSignals.get_code_spectrum_BOCcos(fc * 1Hz, fs * 1Hz, f * 1Hz)
    @test psd_cos4 ≈ psd_cos1
end

@testset "get_code_spectrum for different modulations" begin
    # Test spectrum through get_code_spectrum interface
    gpsl1ca = GPSL1CA()
    gpsl5i = GPSL5I()
    gal_e1b = GalileoE1B()

    # BPSK signals should have peak at DC
    @test get_code_spectrum(gpsl1ca, 0) ≈ 1.0Hz / get_code_frequency(gpsl1ca)
    @test get_code_spectrum(gpsl5i, 0) ≈ 1.0Hz / get_code_frequency(gpsl5i)

    # CBOC has zero at DC
    @test get_code_spectrum(gal_e1b, 0) == 0.0

    # All spectra should be non-negative at non-zero frequencies
    @test get_code_spectrum(gpsl1ca, 100e3) >= 0
    @test get_code_spectrum(gpsl5i, 100e3) >= 0
    @test get_code_spectrum(gal_e1b, 100e3) >= 0
end
