@testset "Common" begin
    @test get_code_center_frequency_ratio(GPSL1(use_gpu = Val(true))) ≈ 1/1540
end
