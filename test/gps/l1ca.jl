@testset "GPS L1 C/A" begin
    gpsl1ca = GPSL1CA()
    @test @inferred(get_band(gpsl1ca)) == L1()
    @test @inferred(get_center_frequency(gpsl1ca)) == 1.57542e9Hz
    @test @inferred(get_code_length(gpsl1ca)) == 1023
    @test @inferred(get_secondary_code_length(gpsl1ca)) == 1
    @test @inferred(get_secondary_code(gpsl1ca)) isa NoSecondaryCode
    @test @inferred(get_code(gpsl1ca, 0, 1)) == 1
    @test @inferred(get_code(gpsl1ca, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(gpsl1ca, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl1ca)) == 50Hz
    @test @inferred(get_code_frequency(gpsl1ca)) == 1023e3Hz
    @test get_code.(gpsl1ca, 0:1022, 1) == L1_SAT1_CODE
    @test @inferred(get_modulation(gpsl1ca)) == GNSSSignals.LOC()
    @test get_signal_name(gpsl1ca) == "GPS L1 C/A"

    @test GNSSSignals.get_code_factor(gpsl1ca) == 1

    @test get_code_spectrum(gpsl1ca, 0) ≈ 1.0Hz / get_code_frequency(gpsl1ca)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(gpsl1ca, m * get_code_frequency(gpsl1ca)) == 0
        @test get_code_spectrum(gpsl1ca, -m * get_code_frequency(gpsl1ca)) == 0
    end
    @test sum(get_code_spectrum.(gpsl1ca, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5

    @test get_code_center_frequency_ratio(gpsl1ca) ≈ 1 / 1540
end
