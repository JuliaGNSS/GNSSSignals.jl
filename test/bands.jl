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

@testset "get_band_id" begin
    # On the band itself (instance and type).
    @test @inferred(get_band_id(L1())) === :L1
    @test @inferred(get_band_id(L5())) === :L5
    @test @inferred(get_band_id(L1)) === :L1

    # On a signal: id follows the signal's band, so everything on one carrier
    # shares an id regardless of constellation.
    @test @inferred(get_band_id(GPSL1CA())) === :L1
    @test @inferred(get_band_id(GalileoE1B())) === :L1
    @test @inferred(get_band_id(GPSL5I())) === :L5
    @test get_band_id(GPSL1CA()) === get_band_id(GalileoE1B()) === get_band_id(L1())

    # Type-level signal dispatch works without constructing the signal.
    @test @inferred(get_band_id(GPSL1CA)) === :L1
end
