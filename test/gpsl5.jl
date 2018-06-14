@testset "shift_registers" begin
reg_xb = GNSSSignals.init_shift_register([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                                                    XB -> (XB[1] + XB[3] + XB[4] + XB[6] + XB[7] + XB[8] + XB[12] + XB[13]) % 2
                                                    )
    registers = []
    for i = 1:8191
        reg_xb, output, registers = reg_xb()
        if (i in [266, 804, 1559, 3471, 5343])
            @test registers in [ [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0] , [0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1], [1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1] ]
        end
    end
    @test registers == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

end

@testset "GPS L5" begin
gen_sampled_code, get_code_phase = @inferred GNSSSignals.init_gpsl5_code()

code = @inferred gen_sampled_code(0:10229, 10230, 0, 10230, 1)
power = code' * code / 10230
@test power ≈ 1 atol = 1e-4
@test code == L5_SAT1_I5_Code

early = @inferred gen_sampled_code(1:40920, 1023e4, 3.5, 4 * 1023e4, 2)
prompt = @inferred gen_sampled_code(1:40920, 1023e4, 4, 4 * 1023e4, 2)
late = @inferred gen_sampled_code(1:40920, 1023e4, 4.5, 4 * 1023e4, 2)
@test early' * prompt == late' * prompt
end