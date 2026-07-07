@testset "Bands" begin
    @test @inferred(get_center_frequency(L1())) == 1_575_420_000Hz
    @test @inferred(get_center_frequency(L5())) == 1_176_450_000Hz
    @test @inferred(get_center_frequency(B1I())) == 1_561_098_000Hz
    @test @inferred(get_center_frequency(B3I())) == 1_268_520_000Hz
    @test @inferred(get_center_frequency(B2b())) == 1_207_140_000Hz

    @test @inferred(get_band(GPSL1CA())) isa L1
    @test @inferred(get_band(GPSL5I())) isa L5
    @test @inferred(get_band(GalileoE1B())) isa L1
    @test @inferred(get_band(BeiDouB1I())) isa B1I
    @test @inferred(get_band(BeiDouB3I())) isa B3I
    @test @inferred(get_band(BeiDouB2bI())) isa B2b
    # BeiDou B1C / B2a share the GPS L1 / L5 carriers.
    @test @inferred(get_band(BeiDouB1C_D())) isa L1
    @test @inferred(get_band(BeiDouB2aI())) isa L5

    # Center frequency is inherited from the band — same Hz for everything on L1.
    @test get_center_frequency(GPSL1CA()) ==
          get_center_frequency(GalileoE1B()) ==
          get_center_frequency(BeiDouB1C_D()) ==
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

    # BeiDou: signals sharing a carrier share the band id (B1C on L1, B2a on L5),
    # while the BeiDou-only carriers get their own ids.
    @test @inferred(get_band_id(B1I())) === :B1I
    @test @inferred(get_band_id(B3I())) === :B3I
    @test @inferred(get_band_id(B2b())) === :B2b
    @test @inferred(get_band_id(BeiDouB1I())) === :B1I
    @test @inferred(get_band_id(BeiDouB2bI())) === :B2b
    @test get_band_id(BeiDouB1C_D()) === get_band_id(GPSL1CA()) === :L1
    @test get_band_id(BeiDouB2aI()) === get_band_id(GPSL5I()) === :L5

    # Type-level signal dispatch works without constructing the signal.
    @test @inferred(get_band_id(GPSL1CA)) === :L1
end
