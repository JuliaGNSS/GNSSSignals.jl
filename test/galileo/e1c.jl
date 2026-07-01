@testset "Galileo E1C" begin
    galileo_e1c = GalileoE1C()
    a = sqrt(10 / 11)  # BOC(1,1) amplitude
    b = sqrt(1 / 11)   # BOC(6,1) amplitude

    @test @inferred(get_band(galileo_e1c)) == L1()
    @test @inferred(get_center_frequency(galileo_e1c)) == 1.57542e9Hz
    @test @inferred(get_code_length(galileo_e1c)) == 4092
    @test @inferred(get_secondary_code_length(galileo_e1c)) == 25
    @test @inferred(get_secondary_code(galileo_e1c)) isa SharedSecondaryCode{25}
    @test @inferred(get_data_frequency(galileo_e1c)) == 0Hz
    @test @inferred(get_code_frequency(galileo_e1c)) == 1023e3Hz
    @test get_signal_name(galileo_e1c) == "Galileo E1C"

    # E1C uses CBOC(−): the BOC(6,1) component is in anti-phase, so at chip
    # boundaries (subcarrier BOC(1,1) = BOC(6,1) = +1) the subcarrier value
    # is `a − b`, where E1B's would be `a + b`. The primary chip and CS25
    # secondary chip 0 multiply on top.
    prim0 = get_codes(galileo_e1c)[1, 1]
    sec0 = GNSSSignals.secondary_value(get_secondary_code(galileo_e1c), 1, 0)
    @test @inferred(get_code(galileo_e1c, 0, 1)) ≈ prim0 * sec0 * (a - b)
    @test @inferred(get_code(galileo_e1c, 0.0, 1)) ≈ prim0 * sec0 * (a - b)
    # phase 0.5: BOC(1,1) flips to −1, BOC(6,1) stays +1 → subcarrier −a − b.
    @test @inferred(get_code(galileo_e1c, 0.5, 1)) ≈ prim0 * sec0 * (-a - b)
    @test @inferred(get_code(galileo_e1c, 1.0, 1)) ≈
          get_codes(galileo_e1c)[2, 1] * sec0 * (a - b)
    @test @inferred(GNSSSignals.get_code_unsafe(galileo_e1c, 0.0, 1)) ≈
          prim0 * sec0 * (a - b)

    @test GNSSSignals.get_code_factor(galileo_e1c) == 1

    # PSD is independent of the BOC(6,1) sign, so E1C and E1B share a spectrum.
    @test get_code_spectrum(galileo_e1c, 0) == 0.0
    @test get_code_spectrum(galileo_e1c, 100e3) == get_code_spectrum(GalileoE1B(), 100e3)
    @test get_code_spectrum(galileo_e1c, 2 * get_code_frequency(galileo_e1c)) == 0.0
    @test get_code_spectrum(galileo_e1c, -2 * get_code_frequency(galileo_e1c)) == 0.0

    @test get_code_center_frequency_ratio(galileo_e1c) ≈ 1 / 1540

    # CBOC output is floating point (the BOC1/BOC2 amplitudes are irrational).
    @test get_code_type(galileo_e1c) === Float32

    # The embedded per-signal LUT is always populated (the anti-phase CBOC(−)
    # composite is baked into it; see `build_signal_lut` / `gen_code!`).
    @test galileo_e1c.lut isa GNSSSignals.SignalLUT
end

@testset "Galileo E1C secondary code (CS25)" begin
    galileo_e1c = GalileoE1C()
    sec = get_secondary_code(galileo_e1c)

    # CS25 = "0011100000001010110110010" (Galileo OS SIS ICD §2.3.4),
    # mapped to ±1 with bit 1 → +1, bit 0 → −1 (matching the primary-code
    # sign convention).
    cs25_bits = "0011100000001010110110010"
    expected = [c == '1' ? Int8(1) : Int8(-1) for c in cs25_bits]
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:24] == expected

    # Shared across PRNs and wraps modulo 25.
    @test GNSSSignals.secondary_value(sec, 1, 0) == GNSSSignals.secondary_value(sec, 30, 0)
    @test GNSSSignals.secondary_value(sec, 1, 25) == GNSSSignals.secondary_value(sec, 1, 0)

    # Tiered-code overlay: at integer chip phases the CBOC(−) subcarrier is
    # the constant `a − b`, so each 4092-chip primary period is just the
    # primary code scaled by its CS25 chip. Same idea as the GPS L5-I
    # Neuman-Hofman test.
    a = sqrt(10 / 11)
    b = sqrt(1 / 11)
    primary = Float64.(get_codes(galileo_e1c)[:, 1])
    code = get_code.(galileo_e1c, 0:(25 * 4092 - 1), 1)
    for k = 1:25
        period = code[(1 + 4092 * (k - 1)):(4092 * k)]
        @test period ≈ primary .* expected[k] .* (a - b)
    end
    # The tiered code repeats after 25 primary periods (100 ms).
    @test code ≈ get_code.(galileo_e1c, (25 * 4092):(2 * 25 * 4092 - 1), 1)
end

@testset "Galileo E1C primary codes" begin
    e1c = GalileoE1C()
    e1b = GalileoE1B()
    # E1B and E1C use distinct memory codes (Galileo OS SIS ICD Annex C).
    @test get_codes(e1c) != get_codes(e1b)
    @test size(get_codes(e1c)) == (4092, 50)
    # Memory codes are ±1.
    @test all(x -> x == 1 || x == -1, get_codes(e1c))
end

@testset "Galileo E1C BOC(1,1) approximation" begin
    e1c = GalileoE1C_BOC11()
    @test @inferred(get_band(e1c)) == L1()
    @test @inferred(get_center_frequency(e1c)) == 1.57542e9Hz
    @test @inferred(get_code_length(e1c)) == 4092
    @test @inferred(get_secondary_code(e1c)) isa SharedSecondaryCode{25}
    @test @inferred(get_secondary_code_length(e1c)) == 25
    @test @inferred(get_data_frequency(e1c)) == 0Hz
    @test @inferred(get_code_frequency(e1c)) == 1023e3Hz
    @test get_signal_name(e1c) == "Galileo E1C (BOC(1,1) approximation)"

    # Modulation differs from full E1C's CBOC.
    @test @inferred(get_modulation(e1c)) == BOCsin(1, 1)
    @test GNSSSignals.get_code_factor(e1c) == 1

    # Output element type stays integer (the whole point of dropping the
    # BOC(6,1) component).
    @test get_code_type(e1c) === Int16

    # Primary and secondary codes match the full GalileoE1C chip-for-chip.
    @test get_codes(e1c) == get_codes(GalileoE1C())
    @test [GNSSSignals.secondary_value(get_secondary_code(e1c), 1, k) for k = 0:24] ==
          [GNSSSignals.secondary_value(get_secondary_code(GalileoE1C()), 1, k) for k = 0:24]
end

@testset "Galileo E1C gen_code! sign matches get_code" begin
    # The embedded LUT emits an Int8 replica: ±1 for the BOC(1,1) approximation
    # and the multi-level integer approximation of the sqrt-power CBOC amplitudes
    # for full E1C — so it matches `get_code`'s SIGN, not its (irrational) value.
    # Sample at the sub-chip rate `fc · subchip_factor` (12.276 MHz for the CBOC
    # P=12, 2.046 MHz for the BOC(1,1) P=2) so every sample lands on a sub-chip
    # start and the drift-free DDA agrees with `get_code`'s exact floor bit-for-bit.
    for sig in (GalileoE1C(), GalileoE1C_BOC11())
        cf = get_code_frequency(sig)
        fs = cf * sig.lut.subchip_factor
        # Span several primary periods to cross CS25 secondary boundaries,
        # and start past the first period so the secondary index offset is
        # exercised too.
        for (samples, start_phase, prn) in
            ((20000, 0.0, 1), (4000, 3.456, 7), (200, 4092.0 * 3 + 1.5, 12))
            code = zeros(Int8, samples)
            gen_code!(code, sig, prn, fs, cf, start_phase)
            phase = (0:samples-1) .* cf ./ fs .+ start_phase
            @test sign.(Int.(code)) == sign.(get_code.(sig, phase, prn))
        end
    end
end

@testset "GalileoE1C_BOC11 sample-stream matches PocketSDR reference (PRN 1, 12.276 MHz)" begin
    # Independent cross-check against the BOC(1,1) replica produced by
    # PocketSDR's `gen_code_E1C` at 2 samples per chip, upsampled by 6× to
    # 12 samples per chip, with PocketSDR's CS25 secondary code applied per
    # primary period. The fixture spans the full 25-period (100 ms) tiered
    # code, so it exercises the primary memory codes, the CS25 secondary
    # overlay and its transitions, and the BOC(1,1) alignment in one shot.
    # At `sampling_frequency = 12 × code_frequency` the primary-chip and BOC
    # half-cycle boundaries are both sample-aligned, so the comparison is
    # bit-exact.
    #
    # Unlike the E1B fixture, no sign flip is applied at generation time.
    # PocketSDR encodes both the primary code and CS25 with the opposite
    # chip-sign convention from Julia; for E1C those two flips cancel
    # (Julia's primary × secondary equals PocketSDR's primary × secondary),
    # so PocketSDR's stream matches Julia's directly. That cancellation is
    # itself a check that the primary and secondary sign conventions are
    # mutually consistent.
    fixture_dir = joinpath(@__DIR__, "fixtures")
    n_samples = 25 * 12 * 4092
    ref = _load_packed_hex_fixture(
        joinpath(fixture_dir, "galileo_e1c_boc11_prn1_fs12chip.hex.gz"),
        n_samples,
    )

    sampling_rate = 12.276e6Hz
    code_rate = 1023e3Hz

    # The embedded LUT emits Int8 ±1; byte-exact to the ±1 fixture at this chip-aligned rate.
    buf = zeros(Int8, n_samples)
    gen_code!(buf, GalileoE1C_BOC11(), 1, sampling_rate, code_rate, 0.0, 0)
    @test buf == ref
end
