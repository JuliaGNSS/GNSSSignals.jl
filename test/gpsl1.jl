@testset "GPS L1" begin
    gps_l1 = @inferred GPSL1()

    code = @inferred gen_code(gps_l1, 0:1022, 1023, 0, 1023, 1)
    power = code.' * code / 1023
    @test power ≈ 1
    @test code == SAT1_CODE

    early = @inferred gen_code(gps_l1, 1:4092, 1023e3, 3.5, 4 * 1023e3, 1)
    prompt = @inferred gen_code(gps_l1, 1:4092, 1023e3, 4, 4 * 1023e3, 1)
    late = @inferred gen_code(gps_l1, 1:4092, 1023e3, 4.5, 4 * 1023e3, 1)
    @test early' * prompt == late' * prompt
end