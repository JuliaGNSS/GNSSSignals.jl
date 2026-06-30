@testset "Common functions for $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GalileoE1C(),
    GPSL1CA(),
    GPSL5I(),
]
    if get_modulation(signal) isa GNSSSignals.CBOC
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
    @test min_bits_for_code_length(GalileoE1C()) == 17  # 4092 * 25 = 102300 requires 17 bits
end

@testset "get_signal_name" begin
    @test get_signal_name(GPSL1CA()) == "GPS L1 C/A"
    @test get_signal_name(GPSL5I()) == "GPS L5-I"
    @test get_signal_name(GalileoE1B()) == "Galileo E1B"
    @test get_signal_name(GalileoE1C()) == "Galileo E1C"
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

    show(io, GalileoE1C())
    @test occursin("GalileoE1C", String(take!(io)))
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
    @test get_modulation(GalileoE1C) isa GNSSSignals.CBOC
    # E1C is CBOC(−): the BOC(6,1) component is in anti-phase.
    @test get_modulation(GalileoE1C).boc2_sign == -1
    @test get_modulation(GalileoE1B).boc2_sign == 1
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
           modulation.boc2_sign *
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
    CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11, -1),  # E1C CBOC(−)
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
                                                          [GalileoE1B(), GalileoE1C(), GPSL1CA(), GPSL5I()]
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
    GalileoE1C(),
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
    GalileoE1C(),
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
    GalileoE1C(),
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
