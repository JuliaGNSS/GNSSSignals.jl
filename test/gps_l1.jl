@testset "GPS L1" begin

    @test @inferred(get_center_frequency(GPSL1)) == 1.57542e9Hz
    @test @inferred(get_code_length(GPSL1)) == 1023
    @test @inferred(get_secondary_code_length(GPSL1)) == 1
    @test @inferred(get_code(GPSL1, 0, 1)) == 1
    @test @inferred(get_code(GPSL1, 0.0, 1)) == 1
    @test @inferred(get_code_unsafe(GPSL1, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(GPSL1)) == 50Hz
    @test @inferred(get_code_frequency(GPSL1)) == 1023e3Hz
    @test get_code.(GPSL1, 0:1022, 1) == L1_SAT1_CODE

end
