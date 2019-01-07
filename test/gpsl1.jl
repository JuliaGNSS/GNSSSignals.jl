@testset "GPS L1" begin
    gps_l1 = @inferred GPSL1()

    @test gps_l1.code_length == 1023
    @test gps_l1.code_period == 1ms
    @test gps_l1.code_freq == 1023e3Hz
    @test gps_l1.center_freq == 1.57542e9Hz
    @test gps_l1.num_prns_per_bit == 20

    @inferred gen_code(gps_l1, 0, 1023, 0, 1023, 1)
    code = gen_code.(Ref(gps_l1), 0:1022, 1023, 0, 1023, 1)
    power = dot(Float64.(code), Float64.(code)) / 1023
    @test power ≈ 1
    @test code == SAT1_CODE

    early = gen_code.(Ref(gps_l1), 1:4092, 1023e3, 3.5, 4 * 1023e3, 1)
    prompt = gen_code.(Ref(gps_l1), 1:4092, 1023e3, 4, 4 * 1023e3, 1)
    late = gen_code.(Ref(gps_l1), 1:4092, 1023e3, 4.5, 4 * 1023e3, 1)
    @test early' * prompt == late' * prompt
end
