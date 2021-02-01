@testset "Common" begin
    @test get_code_center_frequency_ratio(GPSL1()) ≈ 1/1540
end
