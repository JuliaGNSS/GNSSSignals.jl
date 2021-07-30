@testset "GPS L1" begin

    gpsl1 = GPSL1(use_gpu = Val(true))
    @test @inferred(get_center_frequency(gpsl1)) == 1.57542e9Hz
    @test @inferred(get_code_length(gpsl1)) == 1023
    @test @inferred(get_secondary_code_length(gpsl1)) == 1
    CUDA.@allowscalar @test @inferred(get_code(gpsl1, 0, 1)) == 1
    CUDA.@allowscalar @test @inferred(get_code(gpsl1, 0.0, 1)) == 1
    CUDA.@allowscalar @test @inferred(get_code_unsafe(gpsl1, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl1)) == 50Hz
    @test @inferred(get_code_frequency(gpsl1)) == 1023e3Hz
    CUDA.@allowscalar @test gpsl1.codes[:,1] == L1_SAT1_CODE_GPU

end
