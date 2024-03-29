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

@testset "Code generation $(get_system_string(system))" for system in
                                                            [GalileoE1B(), GPSL1(), GPSL5()]
    sampling_rate = 25e6Hz
    samples = 1000
    code = zeros(get_code_type(system), samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 0)
    phase = (0:length(code)-1) * get_code_frequency(system) / sampling_rate
    @test code ≈ get_code.(system, phase, 1)
    @test code ≈ gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 0)
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
    samples = 2500
    code = zeros(Int16, samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 2065)
    phase = (0:samples-1) * get_code_frequency(system) / sampling_rate .+ 2065
    @test code == get_code.(system, phase, 1)
    @test code ==
          gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 2065)
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
