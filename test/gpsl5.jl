@testset "GPS L5" begin
gen_sampled_code, get_code_phase = @inferred GNSSSignals.init_gpsl5_codes(1)

code = @inferred gen_sampled_code(0:10229, 10230, 0, 10230)
power = code' * code / 10230
@test power ≈ 1
@test code == L5_SAT1_I5_Code

early = @inferred gen_sampled_code(1:40920, 1023e4, 3.5, 4 * 1023e4)
prompt = @inferred gen_sampled_code(1:40920, 1023e4, 4, 4 * 1023e4)
late = @inferred gen_sampled_code(1:40920, 1023e4, 4.5, 4 * 1023e4)
@test early' * prompt == late' * prompt
end