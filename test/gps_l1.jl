@testset "GPS L1" begin

    gpsl1 = GPSL1()
    @test @inferred(get_center_frequency(gpsl1)) == 1.57542e9Hz
    @test @inferred(get_code_length(gpsl1)) == 1023
    @test @inferred(get_secondary_code_length(gpsl1)) == 1
    @test @inferred(get_code(gpsl1, 0, 1)) == 1
    @test @inferred(get_code(gpsl1, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(gpsl1, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl1)) == 50Hz
    @test @inferred(get_code_frequency(gpsl1)) == 1023e3Hz
    @test get_code.(gpsl1, 0:1022, 1) == L1_SAT1_CODE
    @test @inferred(get_modulation(gpsl1)) == GNSSSignals.LOC()

    @test GNSSSignals.get_code_factor(gpsl1) == 1

    @test get_code_spectrum(gpsl1, 0) ≈ 1.0Hz / get_code_frequency(gpsl1)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(gpsl1, m * get_code_frequency(gpsl1)) == 0
        @test get_code_spectrum(gpsl1, -m * get_code_frequency(gpsl1)) == 0
    end
    @test sum(get_code_spectrum.(gpsl1, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5

    @test get_code_center_frequency_ratio(gpsl1) ≈ 1/1540
end
