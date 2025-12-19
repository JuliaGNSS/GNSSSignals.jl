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
