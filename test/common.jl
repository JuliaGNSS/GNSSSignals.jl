@testset "Common functions for $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    if typeof(signal) <: GalileoE1B
        @test get_code_type(signal) == Float32
    else
        @test get_code_type(signal) == Int16
    end
    @test get_codes(signal) == signal.codes
end

@testset "min_bits_for_code_length" begin
    @test min_bits_for_code_length(GPSL1CA()) == 10  # 1023 requires 10 bits
    @test min_bits_for_code_length(GPSL5I()) == 17  # 10230 * 10 = 102300 requires 17 bits
    @test min_bits_for_code_length(GalileoE1B()) == 12  # 4092 requires 12 bits
end

@testset "get_signal_name" begin
    @test get_signal_name(GPSL1CA()) == "GPS L1 C/A"
    @test get_signal_name(GPSL5I()) == "GPS L5-I"
    @test get_signal_name(GalileoE1B()) == "Galileo E1B"
end

@testset "SecondaryCode dispatch" begin
    # L5-I has a SharedSecondaryCode of length 10 (NH10).
    gpsl5i_sec = get_secondary_code(GPSL5I())
    @test gpsl5i_sec isa SharedSecondaryCode{10}
    @test GNSSSignals.secondary_code_length(gpsl5i_sec) == 10
    # NH10 = (1, 1, 1, 1, -1, -1, 1, -1, 1, -1); prn is ignored for SharedSecondaryCode.
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 0) == 1
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 4) == -1
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 9) == -1
    # Index wraps modulo length.
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 10) == 1

    # L1 C/A and E1B have NoSecondaryCode; secondary_value returns `true`
    # (so multiplication is a true no-op preserving eltype).
    for sig in (GPSL1CA(), GalileoE1B())
        sec = get_secondary_code(sig)
        @test sec isa NoSecondaryCode
        @test GNSSSignals.secondary_code_length(sec) == 1
        @test GNSSSignals.secondary_value(sec, 1, 0) === true
        @test GNSSSignals.secondary_value(sec, 7, 42) === true
    end
end

@testset "Base.show for GNSS signals" begin
    io = IOBuffer()
    show(io, GPSL1CA())
    @test occursin("GPSL1CA", String(take!(io)))

    show(io, GPSL5I())
    @test occursin("GPSL5I", String(take!(io)))

    show(io, GalileoE1B())
    @test occursin("GalileoE1B", String(take!(io)))
end

@testset "Broadcasting GNSS signals" begin
    gpsl1ca = GPSL1CA()
    # Test that signals can be broadcast
    phases = [0.0, 1.0, 2.0]
    result = get_code.(gpsl1ca, phases, 1)
    @test length(result) == 3
    @test all(x -> x ∈ [-1, 1], result)
end

@testset "get_modulation type dispatch" begin
    @test get_modulation(GPSL1CA) == GNSSSignals.LOC()
    @test get_modulation(GPSL5I) == GNSSSignals.LOC()
    @test get_modulation(GalileoE1B) isa GNSSSignals.CBOC
end

function conventional_gen_subcarrier(
    code_length,
    modulation::BOCsin,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
)
    return iseven.(
        floor.(
            Int,
            ((0:code_length-1) .+ start_index) .* 2 .* modulation.m .* code_frequency ./
            sampling_frequency .+ start_phase .* 2 .* modulation.m,
        )
    ) .* 2 .- 1
end

function conventional_gen_subcarrier(
    code_length,
    modulation::BOCcos,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
)
    conventional_gen_subcarrier(
        code_length,
        BOCsin(modulation.m, modulation.n),
        sampling_frequency,
        code_frequency,
        start_phase + 0.25,
        start_index,
    )
end

function conventional_gen_subcarrier(
    code_length,
    modulation::CBOC,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
)
    return sqrt(modulation.boc1_power) * conventional_gen_subcarrier(
        code_length,
        modulation.boc1,
        sampling_frequency,
        code_frequency,
        start_phase,
        start_index,
    ) +
           sqrt(1 - modulation.boc1_power) * conventional_gen_subcarrier(
        code_length,
        modulation.boc2,
        sampling_frequency,
        code_frequency,
        start_phase,
        start_index,
    )
end

@testset "Subcarrier generation $modulation" for modulation in [
    BOCsin(2, 1),
    BOCcos(2, 1),
    CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11),
]
    sampling_frequency = 25e6Hz
    num_samples = 4000
    sampled_code = ones(modulation isa CBOC ? Float32 : Int16, num_samples)
    code_frequency = 1023e3Hz

    GNSSSignals.multiply_with_subcarrier!(
        sampled_code,
        modulation,
        sampling_frequency,
        code_frequency,
        3.456,
        -1,
    )

    reference = conventional_gen_subcarrier(
        num_samples,
        modulation,
        BigFloat(sampling_frequency),
        BigFloat(code_frequency),
        BigFloat(3.456),
        -1,
    )

    # Allow a few bounded mismatches: at sampling rates where the BOC
    # half-cycle boundary does not land exactly on a sample boundary,
    # the integer phase accumulator's `ceil`-rounded increment drifts
    # by less than one unit-in-the-last-place per sample relative to
    # the analytical `floor(continuous_phase)` reference. Over the
    # roughly 4 ms covered by 4000 samples this can flip a handful of
    # sign-bit threshold crossings by one sample. The implementation
    # is otherwise spec-aligned and bit-exact at the
    # chip-aligned sampling rate used by the PocketSDR fixture
    # comparison in `test/gps/`.
    mismatches = count(sampled_code .!= reference)
    @test mismatches <= 4
end

@testset "sample_code_tail! applies secondary_start_index for tail-only buffers" begin
    # `sample_code_tail!` had a bug where the secondary-code index it
    # used was computed only from the buffer-local chip offset and
    # ignored `secondary_start_index`. Effect: any buffer small enough
    # to fit entirely in the tail (a few dozen samples at moderate
    # oversampling) silently fell back to secondary chip 0 regardless
    # of `start_phase`. Larger buffers passed through the main worker
    # loop which used the right index, so the bug only surfaced with
    # `start_phase ≥ primary_length` AND small `N`.
    #
    # Test: for each signal that has a secondary code, ask for a
    # small buffer at `start_phase = k * primary_length` and verify
    # that every sample equals
    # `secondary[k] / secondary[0] × gen_code!(start_phase = 0)`.
    # Keep `N` well below one primary period so no secondary-bit
    # transition happens inside the buffer.

    # GPS L1C-P (PerPRNSecondaryCode): cover offsets where the
    # overlay sign flips relative to chip 0.
    let signal = GPSL1C_P(),
        prn = 2,
        primary = get_code_length(signal),
        sr = 12.276e6Hz,
        cf = 1023e3Hz,
        sec = GNSSSignals.get_secondary_code(signal)

        s0 = GNSSSignals.secondary_value(sec, prn, 0)
        for sec_offset in 0:5
            for N in (12, 24, 100)
                buf0 = zeros(Int16, N)
                buf1 = zeros(Int16, N)
                gen_code!(buf0, signal, prn, sr, cf, 0.0, 0)
                gen_code!(buf1, signal, prn, sr, cf, Float64(sec_offset * primary), 0)
                sk = GNSSSignals.secondary_value(sec, prn, sec_offset)
                if Int(sk) * Int(s0) > 0
                    @test buf1 == buf0
                else
                    @test buf1 == .-buf0
                end
            end
        end
    end

    # GPS L5-I (SharedSecondaryCode = Neuman-Hofman NH10).
    let signal = GPSL5I(),
        prn = 1,
        primary = get_code_length(signal),
        sr = 25e6Hz,
        cf = 10230e3Hz,
        sec = GNSSSignals.get_secondary_code(signal)

        s0 = GNSSSignals.secondary_value(sec, prn, 0)
        for sec_offset in 0:9
            for N in (12, 24, 100)
                buf0 = zeros(Int16, N)
                buf1 = zeros(Int16, N)
                gen_code!(buf0, signal, prn, sr, cf, 0.0, 0)
                gen_code!(buf1, signal, prn, sr, cf, Float64(sec_offset * primary), 0)
                sk = GNSSSignals.secondary_value(sec, prn, sec_offset)
                if Int(sk) * Int(s0) > 0
                    @test buf1 == buf0
                else
                    @test buf1 == .-buf0
                end
            end
        end
    end
end

@testset "gen_code! error paths" begin
    gpsl1ca = GPSL1CA()

    # Test sampling frequency too low error
    code = zeros(Int16, 100)
    @test_throws "The sampling frequency must be larger than the code frequency" gen_code!(
        code,
        gpsl1ca,
        1,
        500e3Hz,  # Too low - less than code frequency
        get_code_frequency(gpsl1ca),
    )
end

# Above num_inner_iterations = 64 the dispatcher falls back to a @simd ivdep
# generic worker. Exercise that path so any regression there is caught.
@testset "High oversampling falls back to generic worker" begin
    signal = GPSL1CA()
    # frequency_ratio = 200e6 / 1.023e6 ≈ 195 → num_inner = 196 → generic path
    sampling_rate = 200e6Hz
    samples = 2000
    code = zeros(Int16, samples)
    gen_code!(code, signal, 1, sampling_rate, get_code_frequency(signal), 0.0, 0)
    phase = (0:samples-1) * get_code_frequency(signal) / sampling_rate
    @test code == get_code.(signal, phase, 1)
end

@testset "Code generation $(get_signal_name(signal))" for signal in
                                                          [GalileoE1B(), GPSL1CA(), GPSL5I()]
    sampling_rate = 25e6Hz
    samples = 4000
    code = zeros(get_code_type(signal), samples)
    code = gen_code!(code, signal, 1, sampling_rate, get_code_frequency(signal), 0)
    phase = (0:length(code)-1) * get_code_frequency(signal) / sampling_rate
    @test code ≈ get_code.(signal, phase, 1)
    @test code ≈ gen_code(samples, signal, 1, sampling_rate, get_code_frequency(signal), 0)
end

@testset "Small code generation $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    sampling_rate = 25e6Hz
    samples = 100
    code = zeros(get_code_type(signal), samples)
    code = gen_code!(code, signal, 1, sampling_rate, get_code_frequency(signal), 3.5)
    phase = (0:length(code)-1) * get_code_frequency(signal) / sampling_rate .+ 3.5
    @test code ≈ get_code.(signal, phase, 1)
end

# Regression test for the fixed-point overflow on very short outputs
# (https://github.com/JuliaGNSS/GNSSSignals.jl/issues/63): a handful of
# samples makes `ndigits(length)` tiny, so the fixed-point exponent grows
# until `frequency_ratio * 2^exponent` overflows Int64. High sampling rate
# (large frequency_ratio) is the worst case. Such tiny requests arise in
# multi-signal tracking when two signals' code-block boundaries nearly
# coincide.
@testset "Very short code generation does not overflow $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    sampling_rate = 25e6Hz
    for samples in (1, 2, 3, 5)
        code = zeros(get_code_type(signal), samples)
        code = @test_nowarn gen_code!(
            code, signal, 1, sampling_rate, get_code_frequency(signal), 3.5,
        )
        phase = (0:samples-1) * get_code_frequency(signal) / sampling_rate .+ 3.5
        @test code ≈ get_code.(signal, phase, 1)
    end
end

@testset "Code generation for different units" begin
    sampling_rate = 25MHz
    signal = GPSL1CA()
    samples = 1000
    code = zeros(get_code_type(signal), samples)
    code = gen_code!(code, signal, 1, sampling_rate, get_code_frequency(signal), 0)
    phase = (0:length(code)-1) * get_code_frequency(signal) / sampling_rate
    @test code ≈ get_code.(signal, phase, 1)
    @test code ≈ gen_code(samples, signal, 1, sampling_rate, get_code_frequency(signal), 0)
end

@testset "Code generation with start_phase bigger than code_length" begin
    signal = GPSL1CA()
    sampling_rate = 2.5e6Hz
    num_samples = 4000
    code = zeros(Int16, num_samples)
    code = gen_code!(code, signal, 1, sampling_rate, get_code_frequency(signal), 2065)
    phase = (0:num_samples-1) * get_code_frequency(signal) / sampling_rate .+ 2065
    @test code == get_code.(signal, phase, 1)
    @test code ==
          gen_code(num_samples, signal, 1, sampling_rate, get_code_frequency(signal), 2065)
end

@testset "Code generation $(get_signal_name(signal)) with different index" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    sampling_rate = 25e6Hz
    samples = 4002
    code = zeros(get_code_type(signal), samples)
    code = gen_code!(code, signal, 1, sampling_rate, get_code_frequency(signal), 0.0, -1)
    phase = (-1:4000) * get_code_frequency(signal) / sampling_rate
    @test code ≈ get_code.(signal, phase, 1)
    @test code ≈
          gen_code(samples, signal, 1, sampling_rate, get_code_frequency(signal), 0.0, -1)
end

@testset "Code generation with large start_phase (overflow bug fix)" begin
    # Bug: Large start_phase values cause integer overflow in fixed-point arithmetic,
    # resulting in negative array indices and memory corruption.
    # This was discovered during GPS tracking after ~60 seconds when start_phase
    # accumulates to ~14000 chips.

    signal = GPSL1CA()
    sampling_freq = 5.0e6Hz
    code_frequency = 1.0230022937236385e6Hz  # Actual value from tracking crash
    start_phase = 14113.513288791713  # Large value that caused overflow
    start_index_shift = -2

    sampled_code = Vector{Int16}(undef, 1023)

    # This should not throw a BoundsError
    @test_nowarn gen_code!(
        sampled_code,
        signal,
        1,
        sampling_freq,
        code_frequency,
        start_phase,
        start_index_shift,
    )

    # Verify the result is correct by comparing with BigFloat reference implementation
    phase = (start_index_shift:start_index_shift+length(sampled_code)-1) *
            BigFloat(code_frequency / sampling_freq) .+ BigFloat(start_phase)
    @test sampled_code == get_code.(signal, phase, 1)
end

@testset "Code generation with negative start_phase (overflow bug fix)" begin
    # Bug: Negative start_phase also causes integer overflow in fixed-point arithmetic.

    signal = GPSL1CA()
    sampling_freq = 5.0e6Hz
    code_frequency = get_code_frequency(signal) + 4000Hz * get_code_center_frequency_ratio(signal)
    start_phase = -1000.0  # Negative value that caused overflow
    start_index_shift = 0

    sampled_code = Vector{Int16}(undef, 1023)

    # This should not throw a BoundsError
    @test_nowarn gen_code!(
        sampled_code,
        signal,
        1,
        sampling_freq,
        code_frequency,
        start_phase,
        start_index_shift,
    )

    # Verify the result is correct by comparing with BigFloat reference implementation
    phase = (start_index_shift:start_index_shift+length(sampled_code)-1) *
            BigFloat(code_frequency / sampling_freq) .+ BigFloat(start_phase)
    @test sampled_code == get_code.(signal, phase, 1)
end

# Exactness / no-drift regression for the integer-DDA `sample_code!`
# (https://github.com/JuliaGNSS/GNSSSignals.jl/issues/66). The previous running
# binary fixed-point accumulator reserved `ndigits(length)` headroom bits, which
# coarsened the realised code frequency on long buffers and let the phase drift
# against the true frequency. The current DDA reduces the accumulator every chip,
# so `fixed_point` is maximal and length-independent: drift-free for any buffer.
#
# Reference: treat the Float64 frequency *values* as exact and compute the chip
# with exact rational arithmetic (no rounding of the ratio, exact floor).
@testset "sample_code! is drift-free and exact (issue #66)" begin
    # Exact chip for sample `j` (0-based): floor of the absolute phase in
    # `Rational{BigInt}` from the Float64 frequency values.
    function exact_codes(signal, prn, N, sampling, codefreq, start_phase, shift)
        s = Rational{BigInt}(Float64(sampling / 1Hz))
        c = Rational{BigInt}(Float64(codefreq / 1Hz))
        sp = Rational{BigInt}(Float64(start_phase))
        step = c // s
        out = Vector{Int}(undef, N)
        for j = 0:N-1
            chip = floor(BigInt, sp + (shift + j) * step)
            out[j+1] = get_code(signal, Float64(chip), prn)
        end
        out
    end

    cf_l1 = get_code_frequency(GPSL1CA())
    cf_l5 = get_code_frequency(GPSL5I())

    # Doppler-shifted (non-integer) code rates: bit-exact at every oversampling,
    # length, start phase and index shift. These never land a sample on an exact
    # integer chip boundary within a realistic buffer (the rationalised period is
    # ~2^53 samples), so there are no boundary ties — only drift could cause a
    # mismatch, and the maximal `fixed_point` keeps that far below one sample.
    # This is the realistic GNSS case (the code rate is always Doppler-shifted).
    cf_l1_dop = 1.0230022937236385e6Hz
    cf_l5_dop = 1.0229734476348808e7Hz
    doppler_cases = [
        # (signal, sampling, code_frequency, N, start_phase, shift)
        (GPSL1CA(), 5e6Hz, cf_l1_dop, 2_000_000, 0.0, 0),        # long buffer (drift)
        (GPSL1CA(), 1.2 * cf_l1_dop, cf_l1_dop, 100_000, 0.0, 0),# ~1.2x oversampling
        (GPSL1CA(), 5e6Hz, cf_l1_dop, 1_000, 3.456, -1),         # fractional phase + shift
        (GPSL1CA(), 5e6Hz, cf_l1_dop, 1_023, -1000.0, 3),        # negative phase
        (GPSL5I(), 30e6Hz, cf_l5_dop, 500_000, 12.7, 1),
        (GPSL5I(), 2.4 * cf_l5_dop, cf_l5_dop, 100_000, 0.0, -2),
        (GPSL1CA(), 200e6Hz, cf_l1_dop, 4_000, 0.0, 0),          # generic worker
        (GPSL5I(), 25e6Hz, cf_l5_dop, 3, 7.5, -2),               # tail-only tiny buffer
    ]
    for (signal, sampling, cf, N, sp, shift) in doppler_cases
        buf = zeros(get_code_type(signal), N)
        gen_code!(buf, signal, 1, sampling, cf, sp, shift)
        @test buf == exact_codes(signal, 1, N, sampling, cf, sp, shift)
    end

    # Perfectly rational (zero-Doppler) ratios *do* place samples exactly on chip
    # boundaries, where the realised frequency's last-bit rounding may put a
    # sample one index either side. Those ties are isolated (their count tracks
    # the number of coincidences, ∝ length, NOT growing per sample) and bounded
    # well under 0.1% — crucially they do not accumulate, so this is not the drift
    # of the old accumulator. (The exact count depends on `fixed_point`, hence on
    # the architecture-tuned inner-loop padding, so only a loose bound is checked.)
    clean_cases = [
        (GPSL1CA(), 25e6Hz, cf_l1, 2_000_000, 0.0, 0),
        (GPSL5I(), 25e6Hz, cf_l5, 2_000_000, 0.0, 0),
        (GPSL1CA(), 2.5e6Hz, cf_l1, 4_000, 2065.0, 0),
        (GPSL5I(), 9.207e7Hz, cf_l5, 5_000, 2053.0, -3),  # ratio exactly 9
    ]
    for (signal, sampling, cf, N, sp, shift) in clean_cases
        buf = zeros(get_code_type(signal), N)
        gen_code!(buf, signal, 1, sampling, cf, sp, shift)
        mismatches = count(buf .!= exact_codes(signal, 1, N, sampling, cf, sp, shift))
        @test mismatches <= cld(N, 1000)   # < 0.1%, isolated boundary ties
    end
end
