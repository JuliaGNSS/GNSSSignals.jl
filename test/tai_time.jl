import Unitful
using Dates: DateTime
using AstroTime: TAIEpoch, value

@testset "TAITime" begin
    # The J2000 origin and a fraction normalised into [0, 1).
    @test TAITime(DateTime(2000, 1, 1, 12)) == TAITime(0, 0.0)
    @test TAITime(10, 1.25) == TAITime(11, 0.25)
    @test TAITime(10, -0.25) == TAITime(9, 0.75)
    # A fraction just below zero rounds to exactly `1.0` after the floor is taken out;
    # the invariant `fraction ∈ [0, 1)` still holds.
    @test TAITime(10, -1e-17) == TAITime(10, 0.0)
    @test TAITime(10, prevfloat(1.0)).fraction < 1.0

    @test TAITime(0, big(2)^55 + big"0.25") == TAITime(Int64(2)^55, 0.25)

    # A calendar label round-trips at millisecond resolution, before J2000 as well.
    for datetime in (DateTime(2021, 5, 31, 12, 53, 14, 118), DateTime(1980, 1, 6, 0, 0, 19))
        @test DateTime(TAITime(datetime)) == datetime
    end

    t = TAITime(675737594, 0.75)
    @test t + 0.5 == TAITime(675737595, 0.25)
    @test t - 1.0 == TAITime(675737593, 0.75)
    @test (t + 1234.5) - t == 1234.5
    @test 0.5 + t == t + 0.5
    @test t + 500Unitful.ms == t + 0.5
    @test 500Unitful.ms + t == t + 0.5
    @test t - 500Unitful.ms == t - 0.5
    @test t < t + 1e-6
    # `≈` is an absolute tolerance on the difference, not a relative one on ~6.8e8 s.
    @test t ≈ t + 1e-10
    @test !(t ≈ t + 1e-6)
    @test !(t ≈ t + 10.0)
    @test t ≈ t + 1e-6 atol = 1e-5
    # A negative zero fraction is the same epoch as a positive one.
    @test TAITime(10, -0.0) == TAITime(10, 0.0)
    @test TAITime(10, -0.0).fraction === 0.0
    @test sprint(show, TAITime(DateTime(2021, 5, 31, 12, 53, 14, 118))) ==
          "2021-05-31T12:53:14.118 TAI"

    @testset "AstroTime extension" begin
        epoch = TAIEpoch(2021, 5, 31, 12, 53, 14.1183385390904732)
        @test abs(value(TAIEpoch(TAITime(epoch)) - epoch)) <= eps(epoch.fraction) / 2
        # A compensation smaller than our resolution may be lost, but not ignored
        # when it is representable. AstroTime's == alone cannot check this.
        subresolution = TAIEpoch(Int64(42), 0.75, 1e-17)
        @test abs(value(TAIEpoch(TAITime(subresolution)) - subresolution)) <=
              eps(subresolution.fraction) / 2
        compensated = TAIEpoch(Int64(42), 0.25, 2.0^-54)
        @test TAITime(compensated) - TAITime(42, 0.25) == 2.0^-54
        @test TAITime(compensated).fraction == compensated.fraction + compensated.error
        # The same calendar label on both sides, so the two share their origin.
        @test TAITime(DateTime(2021, 5, 31, 12, 53, 14)) ==
              TAITime(TAIEpoch(2021, 5, 31, 12, 53, 14.0))
    end
end
