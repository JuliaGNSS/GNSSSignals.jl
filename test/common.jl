@testset "Common" begin
    @test get_code_center_frequency_ratio(GPSL1()) ≈ 1/1540
end

@testset "Code spectra $(get_system_string(system))" for system = [GPSL1(), GPSL5(), BOCcos(GPSL1(),0,15)]
    @test get_code_spectrum(system, 0) ≈ 1.0Hz/get_code_frequency(system)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(system, m * get_code_frequency(system)) == 0
        @test get_code_spectrum(system, -m * get_code_frequency(system)) == 0
    end
    @test sum(get_code_spectrum.(system, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5
end

@testset "Code spectra $(get_system_string(system))" for system = [GalileoE1B(), BOCcos(GPSL1(),15,2.5)]
    @test get_code_spectrum(system, 0) == 0.0
    @test sum(get_code_spectrum.(system, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5
end

@testset "Code generation $(get_system_string(system))" for system = [GalileoE1B(), GPSL1(), GPSL5()]
    sampling_rate = 25e6Hz
    samples = 1000
    code = zeros(Int16, samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 0)
    phase = (0:length(code) - 1) * get_code_frequency(system) / sampling_rate
    @test code == get_code.(system, phase, 1)
    @test code == gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 0)
end

@testset "Code generation with start_phase bigger than code_length" begin
    system = GPSL1()
    sampling_rate = 2.5e6Hz
    samples = 2500
    code = zeros(Int16, samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 2065)
    phase = (0:samples - 1) * get_code_frequency(system) / sampling_rate .+ 2065
    @test code == get_code.(system, phase, 1)
    @test code == gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 2065)
end

@testset "Code generation $(get_system_string(system)) with different index" for system = [GalileoE1B(), GPSL1(), GPSL5()]
    sampling_rate = 25e6Hz
    samples = 4002
    code = zeros(Int16, samples)
    code = gen_code!(code, system, 1, sampling_rate, get_code_frequency(system), 0.0, -1)
    phase = (-1:4000) * get_code_frequency(system) / sampling_rate
    @test code == get_code.(system, phase, 1)
    @test code == gen_code(samples, system, 1, sampling_rate, get_code_frequency(system), 0.0, -1)
end