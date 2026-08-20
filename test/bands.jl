@testset "Bands" begin
    @test @inferred(get_center_frequency(L1())) == 1_575_420_000Hz
    @test @inferred(get_center_frequency(L5())) == 1_176_450_000Hz
    @test @inferred(get_center_frequency(B1I())) == 1_561_098_000Hz
    @test @inferred(get_center_frequency(B3I())) == 1_268_520_000Hz
    @test @inferred(get_center_frequency(B2b())) == 1_207_140_000Hz
    @test @inferred(get_center_frequency(E6())) == 1_278_750_000Hz

    # Band type entry, as get_band_id / get_band_name accept.
    @test @inferred(get_center_frequency(L1)) == 1_575_420_000Hz
    @test @inferred(get_center_frequency(L2)) == 1_227_600_000Hz
    @test @inferred(get_center_frequency(L5)) == 1_176_450_000Hz
    @test @inferred(get_center_frequency(B1I)) == 1_561_098_000Hz
    @test @inferred(get_center_frequency(B3I)) == 1_268_520_000Hz
    @test @inferred(get_center_frequency(B2b)) == 1_207_140_000Hz
    @test @inferred(get_center_frequency(E6)) == 1_278_750_000Hz
    @test get_center_frequency(L2()) == get_center_frequency(L2)

    @test @inferred(get_band(GPSL1CA())) isa L1
    @test @inferred(get_band(GPSL5I())) isa L5
    @test @inferred(get_band(GalileoE1B())) isa L1
    @test @inferred(get_band(BeiDouB1I())) isa B1I
    @test @inferred(get_band(BeiDouB3I())) isa B3I
    @test @inferred(get_band(BeiDouB2bI())) isa B2b
    @test @inferred(get_band(GalileoE6B())) isa E6
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
    # Galileo E6 gets its own id: 1278.75 MHz is shared with QZSS L6, not with any
    # band above.
    @test @inferred(get_band_id(E6())) === :E6
    @test @inferred(get_band_id(GalileoE6C())) === :E6
    @test get_band_id(BeiDouB1C_D()) === get_band_id(GPSL1CA()) === :L1
    @test get_band_id(BeiDouB2aI()) === get_band_id(GPSL5I()) === :L5

    # Type-level signal dispatch works without constructing the signal.
    @test @inferred(get_band_id(GPSL1CA)) === :L1
end

@testset "get_band_name" begin
    # On the band itself (instance and type).
    @test @inferred(get_band_name(L1())) == "L1"
    @test @inferred(get_band_name(L2())) == "L2"
    @test @inferred(get_band_name(L5())) == "L5"
    @test @inferred(get_band_name(L1)) == "L1"

    # On a signal: the name follows the signal's band, so it is the RF band's
    # label, not the constellation's — Galileo E1B reports "L1", not "E1".
    @test @inferred(get_band_name(GPSL1CA())) == "L1"
    @test @inferred(get_band_name(GalileoE1B())) == "L1"
    @test @inferred(get_band_name(GPSL2CM())) == "L2"
    @test @inferred(get_band_name(GPSL5I())) == "L5"
    @test get_band_name(GalileoE1B()) == get_band_name(GPSL1CA()) == get_band_name(L1())

    # Type-level signal dispatch works without constructing the signal.
    @test @inferred(get_band_name(GalileoE5aI)) == "L5"

    # Display counterpart of the id — same granularity, String instead of Symbol.
    for band in (L1(), L2(), L5(), B1I(), B3I(), B2b(), E6())
        @test get_band_name(band) == String(get_band_id(band))
    end

    # Every band states its name as a literal rather than taking the derived
    # fallback: that is what makes the call allocation-free and constant-foldable,
    # and value equality alone cannot tell the two apart ("L1" == String(:L1)). So
    # check the method dispatch actually reaches, not just what it returns.
    derived = which(get_band_name, Tuple{Type{Band}})
    for B in ALL_BANDS
        @test which(get_band_name, Tuple{Type{B}}) !== derived
        # Stated, but still exactly String() of the id — the two cannot drift apart.
        @test get_band_name(B) == String(get_band_id(B))
    end
end

# A band declared outside the package: it defines no accessor of its own, so both
# the id and the name come from the fallbacks.
struct TestOnlyBand <: Band end

@testset "get_band_name fallback for a new Band" begin
    @test get_band_id(TestOnlyBand) === :TestOnlyBand
    @test get_band_name(TestOnlyBand) == "TestOnlyBand"
    @test get_band_name(TestOnlyBand()) == "TestOnlyBand"
    @test get_band_name(TestOnlyBand) == String(get_band_id(TestOnlyBand))
    # It really is the derived fallback doing the work here.
    @test which(get_band_name, Tuple{Type{TestOnlyBand}}) ===
          which(get_band_name, Tuple{Type{Band}})
end

@testset "every band answers every band accessor" begin
    # `TestOnlyBand` above is deliberately not in ALL_BANDS: the invariant is the
    # package's to keep.
    for B in ALL_BANDS
        @test get_band_id(B) isa Symbol
        @test get_band_name(B) isa String
        # The carrier frequency is the one fact a band has to state for itself —
        # on the band type, the entry get_band_id / get_band_name accept.
        @test get_center_frequency(B) isa Frequency
        @test get_center_frequency(B()) == get_center_frequency(B)
    end
end
