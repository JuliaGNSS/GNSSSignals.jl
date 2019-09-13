@testset "Fast cis" begin
    @inferred GNSSSignals.cis_fast(3)

    @test GNSSSignals.cis_fast.(0:0.01:7) ≈ cis.(0:0.01:7) atol = 1.5
end
