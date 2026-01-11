@testset "Common functions for $(get_system_string(system))" for system in [
    GalileoE1B(),
    GPSL1(),
    GPSL5(),
]
    if typeof(system) <: GalileoE1B
        @test get_code_type(system) == Float32
    else
        @test get_code_type(system) == Int16
    end
    @test get_codes(system) == system.codes
end

@testset "min_bits_for_code_length" begin
    @test min_bits_for_code_length(GPSL1()) == 10  # 1023 requires 10 bits
    @test min_bits_for_code_length(GPSL5()) == 17  # 10230 * 10 = 102300 requires 17 bits
    @test min_bits_for_code_length(GalileoE1B()) == 12  # 4092 requires 12 bits
end

@testset "get_system_string" begin
    @test get_system_string(GPSL1()) == "GPSL1"
    @test get_system_string(GPSL5()) == "GPSL5"
    @test get_system_string(GalileoE1B()) == "GalileoE1B"
end

@testset "get_secondary_code with phase" begin
    gpsl5 = GPSL5()
    # GPS L5 has a 10-element secondary code: (1, 1, 1, 1, -1, -1, 1, -1, 1, -1)
    @test get_secondary_code(gpsl5, 0.0) == 1      # First period
    @test get_secondary_code(gpsl5, 10230.0) == 1  # Second period
    @test get_secondary_code(gpsl5, 40920.0) == -1 # Fifth period (index 4)
    @test get_secondary_code(gpsl5, 51150.0) == -1 # Sixth period (index 5)

    # GPS L1 has no secondary code (returns 1)
    gpsl1 = GPSL1()
    @test get_secondary_code(gpsl1, 0.0) == 1
    @test get_secondary_code(gpsl1, 1000.0) == 1

    # Galileo E1B has no secondary code (returns 1)
    gal_e1b = GalileoE1B()
    @test get_secondary_code(gal_e1b, 0.0) == 1
    @test get_secondary_code(gal_e1b, 5000.0) == 1
end

@testset "Base.show for GNSS systems" begin
    io = IOBuffer()
    show(io, GPSL1())
    @test occursin("GPSL1", String(take!(io)))

    show(io, GPSL5())
    @test occursin("GPSL5", String(take!(io)))

    show(io, GalileoE1B())
    @test occursin("GalileoE1B", String(take!(io)))
end

@testset "Broadcasting GNSS systems" begin
    gpsl1 = GPSL1()
    # Test that systems can be broadcast
    phases = [0.0, 1.0, 2.0]
    result = get_code.(gpsl1, phases, 1)
    @test length(result) == 3
    @test all(x -> x ∈ [-1, 1], result)
end

@testset "get_modulation type dispatch" begin
    @test get_modulation(GPSL1) == GNSSSignals.LOC()
    @test get_modulation(GPSL5) == GNSSSignals.LOC()
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

    @test sampled_code == conventional_gen_subcarrier(
        num_samples,
        modulation,
        BigFloat(sampling_frequency),
        BigFloat(code_frequency),
        BigFloat(3.456),
        -1,
    )
end

@testset "Failed in Tracking.jl" begin
    code = zeros(Int16, 2502)
    gpsl1 = GPSL1()
    @test_throws "The code frequency 3.069e6 Hz is larger than expected (1031000 Hz)). Please increase the expected maximum Doppler frequency 8000 Hz" gen_code!(
        code,
        gpsl1,
        1,
        7.5e6Hz,
        1023e3Hz * 3,
        2.0,
    )
end

@testset "gen_code! error paths" begin
    gpsl1 = GPSL1()

    # Test sampling frequency too low error
    code = zeros(Int16, 100)
    @test_throws "The sampling frequency must be larger than the code frequency" gen_code!(
        code,
        gpsl1,
        1,
        500e3Hz,  # Too low - less than code frequency
        get_code_frequency(gpsl1),
    )
end

@testset "Code generation $(get_system_string(system))" for system in
                                                            [GalileoE1B(), GPSL1(), GPSL5()]
    sampling_rate = 25e6Hz
    samples = 4000
    code = zeros(get_code_type(system), samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 0)
    phase = (0:length(code)-1) * get_code_frequency(system) / sampling_rate
    @test code ≈ get_code.(system, phase, 1)
    @test code ≈ gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 0)
end

@testset "Small code generation $(get_system_string(system))" for system in [
    GalileoE1B(),
    GPSL1(),
    GPSL5(),
]
    sampling_rate = 25e6Hz
    samples = 100
    code = zeros(get_code_type(system), samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 3.5)
    phase = (0:length(code)-1) * get_code_frequency(system) / sampling_rate .+ 3.5
    @test code ≈ get_code.(system, phase, 1)
end

@testset "Code generation for different units" begin
    sampling_rate = 25MHz
    system = GPSL1()
    samples = 1000
    code = zeros(get_code_type(system), samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 0)
    phase = (0:length(code)-1) * get_code_frequency(system) / sampling_rate
    @test code ≈ get_code.(system, phase, 1)
    @test code ≈ gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 0)
end

@testset "Code generation with start_phase bigger than code_length" begin
    system = GPSL1()
    sampling_rate = 2.5e6Hz
    num_samples = 4000
    code = zeros(Int16, num_samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 2065)
    phase = (0:num_samples-1) * get_code_frequency(system) / sampling_rate .+ 2065
    @test code == get_code.(system, phase, 1)
    @test code ==
          gen_code(num_samples, system, 1, sampling_rate, get_code_frequency(system), 2065)
end

@testset "Code generation $(get_system_string(system)) with different index" for system in [
    GalileoE1B(),
    GPSL1(),
    GPSL5(),
]
    sampling_rate = 25e6Hz
    samples = 4002
    code = zeros(get_code_type(system), samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 0.0, -1)
    phase = (-1:4000) * get_code_frequency(system) / sampling_rate
    @test code ≈ get_code.(system, phase, 1)
    @test code ≈
          gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 0.0, -1)
end

@testset "Code generation with large start_phase (overflow bug fix)" begin
    # Bug: Large start_phase values cause integer overflow in fixed-point arithmetic,
    # resulting in negative array indices and memory corruption.
    # This was discovered during GPS tracking after ~60 seconds when start_phase
    # accumulates to ~14000 chips.

    system = GPSL1()
    sampling_freq = 5.0e6Hz
    code_frequency = 1.0230022937236385e6Hz  # Actual value from tracking crash
    start_phase = 14113.513288791713  # Large value that caused overflow
    start_index_shift = -2

    sampled_code = Vector{Int16}(undef, 1023)

    # This should not throw a BoundsError
    @test_nowarn gen_code!(
        sampled_code,
        system,
        1,
        sampling_freq,
        code_frequency,
        start_phase,
        start_index_shift,
        Val(sampling_freq),
        Val(8000Hz),
    )

    # Verify the result is correct by comparing with BigFloat reference implementation
    phase = (start_index_shift:start_index_shift+length(sampled_code)-1) *
            BigFloat(code_frequency / sampling_freq) .+ BigFloat(start_phase)
    @test sampled_code == get_code.(system, phase, 1)
end

@testset "Code generation with negative start_phase (overflow bug fix)" begin
    # Bug: Negative start_phase also causes integer overflow in fixed-point arithmetic.

    system = GPSL1()
    sampling_freq = 5.0e6Hz
    code_frequency = get_code_frequency(system) + 4000Hz * get_code_center_frequency_ratio(system)
    start_phase = -1000.0  # Negative value that caused overflow
    start_index_shift = 0

    sampled_code = Vector{Int16}(undef, 1023)

    # This should not throw a BoundsError
    @test_nowarn gen_code!(
        sampled_code,
        system,
        1,
        sampling_freq,
        code_frequency,
        start_phase,
        start_index_shift,
        Val(sampling_freq),
        Val(8000Hz),
    )

    # Verify the result is correct by comparing with BigFloat reference implementation
    phase = (start_index_shift:start_index_shift+length(sampled_code)-1) *
            BigFloat(code_frequency / sampling_freq) .+ BigFloat(start_phase)
    @test sampled_code == get_code.(system, phase, 1)
end
