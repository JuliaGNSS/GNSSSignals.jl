# The PSD of a random ±1 code is |P(f)|² / T_c of its chip pulse P, here a chip of equal
# segments with the given signs. With an odd number of half-periods per chip the subcarrier
# flips the pulse's sign from one chip to the next, which leaves a random code's PSD as is.
function chip_pulse_psd(signs, chip_rate, f)
    d = 1 / (chip_rate * length(signs))
    P = sum(s * d * sinc(f * d) * cispi(-2f * (j - 0.5) * d) for (j, s) in enumerate(signs))
    return abs2(P) * chip_rate
end

@testset "BOC spectrum for k = $k subcarrier half-periods per chip" for k in 1:7
    fc = 1.023e6
    fs = k * fc / 2
    sine_signs = [iseven(j) ? 1 : -1 for j in 0:k-1]
    cosine_signs = [mod(q, 4) in (0, 3) ? 1 : -1 for q in 0:2k-1]  # quarter periods
    fs_grid = range(-20fc, 20fc; length = 801) .+ 123.4  # off the 0/0 points
    for (spectrum, signs) in ((GNSSSignals.get_code_spectrum_BOCsin, sine_signs),
                              (GNSSSignals.get_code_spectrum_BOCcos, cosine_signs))
        expected = [chip_pulse_psd(signs, fc, f) for f in fs_grid]
        got = [spectrum(fc, fs, f) for f in fs_grid]
        @test maximum(abs.(got .- expected)) <= 1e-9 * maximum(expected)
        @test spectrum(fc, fs, 0.0) ≈ chip_pulse_psd(signs, fc, 0.0) atol = 1e-20
    end
end

@testset "BOC spectrum rejects a fractional number of half-periods per chip" begin
    @test_throws ArgumentError GNSSSignals.get_code_spectrum_BOCsin(1.023e6, 0.75 * 1.023e6, 0.0)
    @test_throws ArgumentError GNSSSignals.get_code_spectrum_BOCcos(1.023e6, 0.75 * 1.023e6, 0.0)
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
    gal_e1b = GalileoE1B(CBOC)

    # BPSK signals should have peak at DC
    @test get_code_spectrum(gpsl1ca, 0) ≈ 1.0Hz / get_chip_rate(gpsl1ca)
    @test get_code_spectrum(gpsl5i, 0) ≈ 1.0Hz / get_chip_rate(gpsl5i)

    # CBOC has zero at DC
    @test get_code_spectrum(gal_e1b, 0) == 0.0

    # All spectra should be non-negative at non-zero frequencies
    @test get_code_spectrum(gpsl1ca, 100e3) >= 0
    @test get_code_spectrum(gpsl5i, 100e3) >= 0
    @test get_code_spectrum(gal_e1b, 100e3) >= 0
end

@testset "BOC(m, n): the spectrum takes f_c = n · 1.023 MHz and f_s = m · 1.023 MHz" begin
    for modulation in (BOCsin(10, 5), BOCcos(15, 2.5), BOCcos(3, 1), BOCcos(6, 1))
        signal = TestBOCSignal(modulation)
        fc = modulation.n * 1.023e6
        fs = modulation.m * 1.023e6
        # Away from f = 0 and the cosine form's removable singularities.
        for f in range(-30e6, 30e6; length = 2001) .+ 1234.5
            closed_form = if modulation isa BOCsin
                (sinc(f / fc) * tan(π * f / (2fs)))^2 / fc
            else
                (2 * sinc(f / fc) * sin(π * f / (4fs))^2 / cos(π * f / (2fs)))^2 / fc
            end
            @test get_code_spectrum(signal, f * Hz) ≈ closed_form rtol = 1e-10
        end
    end
end
