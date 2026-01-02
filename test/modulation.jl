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
    gpsl1 = GPSL1()
    gpsl5 = GPSL5()
    gal_e1b = GalileoE1B()

    # BPSK systems should have peak at DC
    @test get_code_spectrum(gpsl1, 0) ≈ 1.0Hz / get_code_frequency(gpsl1)
    @test get_code_spectrum(gpsl5, 0) ≈ 1.0Hz / get_code_frequency(gpsl5)

    # CBOC has zero at DC
    @test get_code_spectrum(gal_e1b, 0) == 0.0

    # All spectra should be non-negative at non-zero frequencies
    @test get_code_spectrum(gpsl1, 100e3) >= 0
    @test get_code_spectrum(gpsl5, 100e3) >= 0
    @test get_code_spectrum(gal_e1b, 100e3) >= 0
end
