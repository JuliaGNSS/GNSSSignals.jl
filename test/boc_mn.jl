using Random: MersenneTwister

# BOC(m, n) for n ≠ 1: subcarrier frequency m · 1.023 MHz, chip rate n · 1.023 MHz, the phase
# counted in code chips. No defined signal uses such a modulation, so a test-only signal carries
# it: the modulation is a type parameter (an isbits value), the code a random ±1 sequence.
struct TestBOCSignal{Mod,C<:AbstractMatrix} <: GNSSSignals.AbstractGNSSSignal{C}
    codes::C
    lut::GNSSSignals.SignalLUT
end
# 330 chips: a whole number of subcarrier periods for every modulation below (330 · m/n), and
# of the 33-chip TMBOC pattern, so the code, its subcarrier and the pattern repeat together.
function TestBOCSignal(modulation; code_length = 330)
    primary = rand(MersenneTwister(11), Int8[-1, 1], code_length, 2)
    codes = primary
    lut = GNSSSignals.build_signal_lut(modulation, codes, NoSecondaryCode())
    TestBOCSignal{modulation,typeof(codes)}(codes, lut)
end
GNSSSignals.get_modulation(::Type{<:TestBOCSignal{Mod}}) where {Mod} = Mod
_chip_rate_multiplier(modulation::GNSSSignals.BOC) = modulation.n
_chip_rate_multiplier(modulation::Union{TMBOC,CBOC}) = modulation.boc1.n
GNSSSignals.get_chip_rate(::Type{<:TestBOCSignal{Mod}}) where {Mod} =
    _chip_rate_multiplier(Mod) * 1_023_000Hz
GNSSSignals.get_code_length(signal::TestBOCSignal) = size(signal.codes, 1)
GNSSSignals.get_secondary_code(::TestBOCSignal) = NoSecondaryCode()
GNSSSignals.get_bit_rate(::Type{<:TestBOCSignal}) = 0Hz
GNSSSignals.get_signal_name(::Type{<:TestBOCSignal}) = "test BOC signal"

# Sub-chips per chip the lookup table bakes: 2m/n subcarrier half-periods for sine phase,
# 4m/n quarter periods for cosine phase.
_subchips(m::BOCsin) = 2 * m.m / m.n
_subchips(m::BOCcos) = 4 * m.m / m.n

@testset "BOC(m, n): the code lookup table matches the closed form" begin
    pattern = ntuple(k -> (k - 1) ∈ (0, 4, 6, 29), 33)
    modulations = (
        BOCsin(10, 5), BOCcos(15, 2.5), BOCcos(3, 1), BOCcos(6, 1),
        # 3 half-periods per chip: a chip starts mid-period, so the sign parity alternates.
        BOCsin(1.5, 1),
        TMBOC(BOCsin(1, 2), BOCsin(6, 2), pattern),
    )
    @testset "$modulation" for modulation in modulations
        signal = TestBOCSignal(modulation)
        P = signal.lut.subchip_factor
        modulation isa TMBOC || @test P == _subchips(modulation)
        fc = get_chip_rate(signal)
        L = get_code_length(signal)
        # Sampling frequencies whose step is an exact binary fraction of a sub-chip (3/8, 5/16,
        # 7/64 per sample), started half a step-grid unit off the sub-chip grid: the table's
        # fixed-point step is then exact and no sample falls on a sub-chip edge, so the closed
        # form `get_code` at each sample's phase and the table must agree on every sample.
        for (step_num, log2_den) in ((3, 3), (5, 4), (7, 6))
            step = step_num / 2^log2_den
            fs = fc * P / step
            start_phase = 1 / (2^(log2_den + 1) * P)
            N = ceil(Int, 2 * L * P / step)                   # two code periods
            for prn in (1, 2)
                table = gen_code(N, signal, prn, fs, fc, start_phase)
                closed_form = [get_code(signal, start_phase + i * fc / fs, prn) for i in 0:N-1]
                @test table == closed_form
            end
        end
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

@testset "BOC(m, n): a chip spanning no whole number of sub-chips cannot be baked" begin
    CL = GNSSSignals.CodeLUT
    @test_throws ArgumentError CL.BOC(3, 4)              # 2m/n = 1.5 half-periods per chip
    @test_throws ArgumentError CL.BOCcos(1, 3)           # 4m/n = 4/3 quarter periods per chip
    @test_throws ArgumentError TestBOCSignal(BOCsin(3, 4))
    @test_throws ArgumentError TestBOCSignal(BOCcos(1, 3))
    # The composites bake sine-phased components only.
    pattern = ntuple(k -> isodd(k), 4)
    @test_throws ArgumentError TestBOCSignal(TMBOC(BOCsin(1, 1), BOCcos(2, 1), pattern))
    @test_throws ArgumentError TestBOCSignal(CBOC(BOCcos(1, 1), BOCsin(6, 1), 10 / 11))
    # A code period must hold whole subcarrier periods for a one-period table to repeat: 255
    # chips of BOC(1.5, 1) end half a period in (255 · 1.5 = 382.5 periods).
    @test_throws ArgumentError TestBOCSignal(BOCsin(1.5, 1); code_length = 255)
    @test TestBOCSignal(BOCsin(1.5, 1); code_length = 256).lut.subchip_factor == 3
    # The ones that do divide are accepted, n integer or not.
    @test TestBOCSignal(BOCsin(10, 5)).lut.subchip_factor == 4
    @test TestBOCSignal(BOCcos(15, 2.5)).lut.subchip_factor == 24
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
        :BeiDouB1C_D => 0x24f7e3101ac0ed5c,
        :BeiDouB1C_P => 0x3114f26e31854d9e,
        :BeiDouB1I => 0x3cdfd21947597fd4,
        :BeiDouB2aI => 0xd36ae4e9b84c791a,
        :BeiDouB2aQ => 0x16289042cfb6a822,
        :BeiDouB2bI => 0xef31b0d7a738dca0,
        :BeiDouB3I => 0x2074e6a33068e70a,
        :GPSL1CA => 0xdd13225653451bce,
        :GPSL1C_D => 0xc9320bea1c05c354,
        :GPSL1C_P_TMBOC => 0x6fbf2ce8b86b5e5e,
        :GPSL2CL => 0x1650bbe438e8c388,
        :GPSL2CM => 0x030686fa055d78d6,
        :GPSL5I => 0xe1f91276c27c5948,
        :GPSL5Q => 0x6bb8b1c51563d5da,
        :GalileoE1B_CBOC => 0x0d7ca4bcad554266,
        :GalileoE1B => 0xbe79c0e19b459348,
        :GalileoE1C_CBOC => 0xb39a74bc995eccfa,
        :GalileoE1C => 0x53f2023f5ee5b236,
        :GalileoE5aI => 0xcef044e1447ba264,
        :GalileoE5aQ => 0xdd88612a620117c4,
        :GalileoE5aQP => 0xd698d8875b2bd194,
        :GalileoE5bI => 0x275c2b3f1b41981c,
        :GalileoE5bQ => 0xcde355db4b2c3e0a,
        :GalileoE6B => 0x1e34f18e30df5254,
        :GalileoE6C => 0xcfea8d536ac3b64a,
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
          Set(nameof(S) for S in ALL_SIGNALS if S !== GPSL1C_P)
    @testset "$(nameof(S))" for S in ALL_SIGNALS
        signal = S()
        if haskey(_OUTPUT_CHECKSUMS_BEFORE_BOC_MN, nameof(S))
            @test _signal_output_checksum(signal) == _OUTPUT_CHECKSUMS_BEFORE_BOC_MN[nameof(S)]
        end
        modulation = get_modulation(signal)
        @test all(f -> get_code_spectrum(signal, f) === _spectrum_before(modulation, signal, f),
            range(-20e6, 20e6; length = 4001))
    end
end

