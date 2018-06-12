@testset "Carrier" begin
carrier = @inferred GNSSSignals.gen_carrier(1:1000, 1e3, 30 * pi / 180, 4e6)
power = 1e-3 * carrier' * carrier
@test power ≈ 1

@test carrier ≈ cis.((2 * π * 1e3 / 4e6) .* (1:1000) .+ (30 * pi / 180))
end

@testset "Subcarrier" begin
prn = [1 1 -1 1 -1 -1 1 1 1 -1 1 -1 -1 -1]
code = @inferred GNSSSignals.gen_sat_code(1:1000, 1e6, 4e-7, 4e6, prn)
power = 1e-3 * code' * code
@test power ≈ 1
end
