@testset "Modulation" begin
    galileo_e1b = GalileoE1B(CBOC)
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

    # The cosine subcarrier is the sine one shifted by a quarter subcarrier period,
    # 1/(4m) of the 1.023 MHz periods the phase counts: +1 for the first quarter period,
    # -1 for the middle half, +1 for the last quarter, whatever m is.
    for m in (1, 2, 3, 6)
        quarter = 1 / (4m)
        cos_m = GNSSSignals.BOCcos(m, 1)
        @test GNSSSignals.get_subcarrier_code(cos_m, 0.5 * quarter) == 1
        @test GNSSSignals.get_subcarrier_code(cos_m, 1.5 * quarter) == -1
        @test GNSSSignals.get_subcarrier_code(cos_m, 2.5 * quarter) == -1
        @test GNSSSignals.get_subcarrier_code(cos_m, 3.5 * quarter) == 1
        for φ in range(0, 3; length = 1201)
            # Midpoints between the sine's half-period edges, clear of either rounding way.
            φ_mid = (floor(φ * 4m) + 0.5) / 4m
            @test GNSSSignals.get_subcarrier_code(cos_m, φ_mid) ==
                  GNSSSignals.get_subcarrier_code(GNSSSignals.BOCsin(m, 1), φ_mid + quarter)
        end
    end
    @test GNSSSignals.get_subcarrier_code(GNSSSignals.BOCcos(2, 1), 0) == 1   # integer phase

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

    # CBOC defaults to boc2_sign = +1 ("CBOC(+)"); E1C uses boc2_sign = -1
    # ("CBOC(−)"), negating the BOC(6,1) component.
    @test cboc.boc2_sign == 1
    cboc_minus =
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(6, 1), 10 / 11, -1)
    @test cboc_minus.boc2_sign == -1
    a = sqrt(10 / 11)
    b = sqrt(1 / 11)
    cboc_plus =
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(6, 1), 10 / 11)
    # At a chip boundary BOC(1,1) = BOC(6,1) = +1, so the two signs give a ± b.
    @test GNSSSignals.get_subcarrier_code(cboc_plus, 0.0) ≈ a + b
    @test GNSSSignals.get_subcarrier_code(cboc_minus, 0.0) ≈ a - b

    # The phase counts code chips whatever the chip rate, so the chip is its floor.
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCcos(2, 1), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCcos(2, 2), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCsin(2, 1), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCsin(2, 2), 2.3) == 2

    # Test BPSK floored phase
    @test GNSSSignals.get_floored_phase(GNSSSignals.BPSK(), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BPSK(), 5.9) == 5

    # Test CBOC floored phase
    cboc_test = GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 2), GNSSSignals.BOCsin(6, 2), 10 / 11)
    @test GNSSSignals.get_floored_phase(cboc_test, 2.3) == 2
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

        # The phase counts code chips whatever `n` is.
        boc1_n2 = GNSSSignals.BOCsin(1, 2)
        boc6_n2 = GNSSSignals.BOCsin(6, 2)
        tmboc_n2 = GNSSSignals.TMBOC(boc1_n2, boc6_n2, pattern)
        @test GNSSSignals.get_floored_phase(tmboc_n2, 2.3) == 2
    end

    # `get_code_spectrum`: power-weighted mix of the BOC component
    # spectra, where the weights are the fractions of chips using each.
    @testset "get_code_spectrum" begin
        sig = GPSL1C_P(TMBOC)

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

        # BOC(1,1) and BOC(6,1) put whole subcarrier periods in a chip, so the
        # spectrum vanishes at DC.
        @test get_code_spectrum(sig, 0) == 0.0

        # Spectrum is even.
        @test get_code_spectrum(sig, 500e3) ≈ get_code_spectrum(sig, -500e3)
    end
end

@testset "BOC(m, n): the cosine subcarrier is the sine one n/(4m) chips later" begin
    for modulation in (BOCcos(15, 2.5), BOCcos(3, 1), BOCcos(6, 1), BOCcos(10, 5), BOCcos(1, 1))
        (; m, n) = modulation
        sine = BOCsin(m, n)
        quarters_per_chip = 4m / n
        for φ in range(0, 5; length = 2003)
            # The middle of a quarter period, clear of either sign's edges.
            φ_mid = (floor(φ * quarters_per_chip) + 0.5) / quarters_per_chip
            @test GNSSSignals.get_subcarrier_code(modulation, φ_mid) ==
                  GNSSSignals.get_subcarrier_code(sine, φ_mid + n / 4m)
        end
        # +1 for the first quarter subcarrier period, -1 for the next two, +1 for the last.
        quarter = n / 4m
        @test [GNSSSignals.get_subcarrier_code(modulation, (j + 0.5) * quarter) for j in 0:3] ==
              [1, -1, -1, 1]
    end
end

# Counting the phase in chips and carrying n through changes nothing for the signals the package
# defines, all n = 1. Their generated code is pinned by checksums taken before that change: the
# one-shot `gen_code` and a continuing `code_engine` fill, at three sampling frequencies, two PRNs
# and three start phases each. Exact Int8 output, so the checksums hold on any machine.
_output_checksum(v) = foldl(
    (acc, (i, x)) -> acc * 0x100000001b3 ⊻ (UInt64(i) * 0x9e3779b97f4a7c15 + reinterpret(UInt8, Int8(x))),
    enumerate(v); init = UInt64(length(v)))
function _signal_output_checksum(signal)
    chip_rate = get_chip_rate(signal)
    P = signal.lut.subchip_factor
    acc = UInt64(0)
    for kind in (:gen, :engine), (fs_multiple, prn) in ((1.3, 1), (2.7, 2), (5.1, 1)),
        start_phase in (0.0, 0.37, 123.81)

        fs = fs_multiple * P * chip_rate
        output = if kind === :gen
            gen_code(3000, signal, prn, fs, chip_rate + 1234.5Hz, start_phase)
        else
            engine = GNSSSignals.code_engine(signal, prn, fs, chip_rate - 77Hz; start_phase)
            first_fill, second_fill = zeros(Int8, 1700), zeros(Int8, 1300)
            state = gen_code!(first_fill, engine, GNSSSignals.code_state(engine))
            gen_code!(second_fill, engine, state)
            vcat(first_fill, second_fill)
        end
        acc = acc * 0x100000001b3 ⊻ _output_checksum(output)
    end
    acc
end
const _OUTPUT_CHECKSUMS_BEFORE_BOC_MN = Dict(
        "BeiDouB1C_D" => 0x24f7e3101ac0ed5c,
        "BeiDouB1C_P" => 0x3114f26e31854d9e,
        "BeiDouB1I" => 0x3cdfd21947597fd4,
        "BeiDouB2aI" => 0xd36ae4e9b84c791a,
        "BeiDouB2aQ" => 0x16289042cfb6a822,
        "BeiDouB2bI" => 0xef31b0d7a738dca0,
        "BeiDouB3I" => 0x2074e6a33068e70a,
        "GPSL1CA" => 0xdd13225653451bce,
        "GPSL1C_D" => 0xc9320bea1c05c354,
        "GPSL1C_P(TMBOC)" => 0x6fbf2ce8b86b5e5e,
        "GPSL2CL" => 0x1650bbe438e8c388,
        "GPSL2CM" => 0x030686fa055d78d6,
        "GPSL5I" => 0xe1f91276c27c5948,
        "GPSL5Q" => 0x6bb8b1c51563d5da,
        "GalileoE1B(CBOC)" => 0x0d7ca4bcad554266,
        "GalileoE1B" => 0xbe79c0e19b459348,
        "GalileoE1C(CBOC)" => 0xb39a74bc995eccfa,
        "GalileoE1C" => 0x53f2023f5ee5b236,
        "GalileoE5aI" => 0xcef044e1447ba264,
        "GalileoE5aQ" => 0xdd88612a620117c4,
        "GalileoE5aQP" => 0xd698d8875b2bd194,
        "GalileoE5bI" => 0x275c2b3f1b41981c,
        "GalileoE5bQ" => 0xcde355db4b2c3e0a,
        "GalileoE6B" => 0x1e34f18e30df5254,
        "GalileoE6C" => 0xcfea8d536ac3b64a,
)

# The spectrum as it was wired before: the BOC components took `n · f_c` and `m · f_c` from the
# signal's chip rate f_c, which is the same for n = 1.
_spectrum_before(::GNSSSignals.BPSK, signal, f) = GNSSSignals.get_code_spectrum_BPSK(get_chip_rate(signal), f)
_spectrum_before(m::BOCsin, signal, f) =
    GNSSSignals.get_code_spectrum_BOCsin(m.n * get_chip_rate(signal), m.m * get_chip_rate(signal), f)
_spectrum_before(m::BOCcos, signal, f) =
    GNSSSignals.get_code_spectrum_BOCcos(m.n * get_chip_rate(signal), m.m * get_chip_rate(signal), f)
_spectrum_before(m::CBOC, signal, f) =
    _spectrum_before(m.boc1, signal, f) * m.boc1_power +
    _spectrum_before(m.boc2, signal, f) * (1 - m.boc1_power)
function _spectrum_before(m::TMBOC{B1,B2,N}, signal, f) where {B1,B2,N}
    boc2_frac = count(m.pattern) / N
    _spectrum_before(m.boc1, signal, f) * (1 - boc2_frac) +
    _spectrum_before(m.boc2, signal, f) * boc2_frac
end

@testset "BOC(m, n): existing replica outputs and all signal spectra" begin
    # The new L1C-P BOC(1,1) replica has independent primary/overlay checks in gps/l1c_p.jl.
    @test Set(keys(_OUTPUT_CHECKSUMS_BEFORE_BOC_MN)) ==
          Set(signal_label(S) for S in ALL_SIGNALS if S !== GPSL1C_P)
    @testset "$(signal_label(S))" for S in ALL_SIGNALS
        signal = make_signal(S)
        if haskey(_OUTPUT_CHECKSUMS_BEFORE_BOC_MN, signal_label(S))
            @test _signal_output_checksum(signal) == _OUTPUT_CHECKSUMS_BEFORE_BOC_MN[signal_label(S)]
        end
        modulation = get_modulation(signal)
        @test all(f -> get_code_spectrum(signal, f) === _spectrum_before(modulation, signal, f),
            range(-20e6, 20e6; length = 4001))
    end
end
