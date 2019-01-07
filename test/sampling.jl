@testset "Carrier" begin
    @inferred gen_carrier(1, 1e3, 30π / 180, 4e6)
    carrier = gen_carrier.(1:1000, 1e3, 30π / 180, 4e6)
    power = 1e-3 * carrier' * carrier
    @test power ≈ 1

    @test carrier ≈ cis.((2 * π * 1e3 / 4e6) .* (1:1000) .+ (30π / 180))

    @test @inferred(calc_carrier_phase(1, 100, 0.2, 4e6)) == mod2pi((2 * π * 100 / 4e6) * 1 + 0.2)

    @inferred GNSSSignals.cis_fast(3)

    @test GNSSSignals.cis_fast.(0:0.01:7) ≈ cis.(0:0.01:7) atol = 1.5

    @test gen_carrier_fast.(1:1000, 1e3, 30π / 180, 4e6) ≈ carrier atol = 1.8
end

@testset "Subcarrier" begin
    code = [1, 1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, -1]
    codes = hcat(code, circshift(code, 2))
    @inferred gen_code(1, 1e6, 2, 4e6, codes, 1)
    code_sampled = [gen_code(x, 1e6, 2, 4e6, codes, 1) for x = 1:1000]
    power = 1e-3 * code_sampled' * code_sampled
    @test power ≈ 1

    @test code_sampled == code[1 .+ mod.(floor.(Int, 1e6 ./ 4e6 .* (1:1000) .+ 2), length(code)), 1]

    @test @inferred(calc_code_phase(5000, 1023e3, 100, 4e6, 1023)) == mod(1023e3 / 4e6 * 5000 + 100, 1023)

    @test @inferred(calc_code_phase_unsafe(5000, 1023e3, 100, 4e6)) == 1023e3 / 4e6 * 5000 + 100
end
