@testset "Bands" begin
    @test @inferred(get_center_frequency(L1())) == 1_575_420_000Hz
    @test @inferred(get_center_frequency(L5())) == 1_176_450_000Hz

    @test @inferred(get_band(GPSL1CA())) isa L1
    @test @inferred(get_band(GPSL5I())) isa L5
    @test @inferred(get_band(GalileoE1B())) isa L1

    # Center frequency is inherited from the band — same Hz for everything on L1.
    @test get_center_frequency(GPSL1CA()) ==
          get_center_frequency(GalileoE1B()) ==
          get_center_frequency(L1())

    # Inference smoke test: the dispatch chain collapses to a literal.
    @test @inferred(get_center_frequency(GPSL1CA())) == 1_575_420_000Hz
end
