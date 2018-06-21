@testset "GPS L1" begin
gen_sampled_code, get_code_phase = @inferred GNSSSignals.init_gpsl1_codes()

code = @inferred gen_sampled_code(0:1022, 1023, 0, 1023, 1)
f_code = float(code)
power = f_code' * f_code / 1023
@test power ≈ 1
@test code == SAT1_CODE

early = @inferred gen_sampled_code(1:4092, 1023e3, 3.5, 4 * 1023e3, 1)
prompt = @inferred gen_sampled_code(1:4092, 1023e3, 4, 4 * 1023e3, 1)
late = @inferred gen_sampled_code(1:4092, 1023e3, 4.5, 4 * 1023e3, 1)
@test early' * prompt == late' * prompt
end