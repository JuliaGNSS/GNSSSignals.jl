using Unitful

# Express a linear power (W) back in dBW for readable comparison with the ICD.
_to_dbw(p) = 10 * log10(ustrip(u"W", p))

@testset "get_min_received_power" begin
    # Per-component values on the Galileo 0 dBi RHCP reference. GPS values equal
    # their IS-GPS ICD figures (net-zero re-referencing).
    @test _to_dbw(get_min_received_power(GPSL1CA())) ≈ -158.5 atol = 1e-6
    @test _to_dbw(get_min_received_power(GPSL5I())) ≈ -157.9 atol = 1e-6
    @test _to_dbw(get_min_received_power(GPSL5Q())) ≈ -157.9 atol = 1e-6
    @test _to_dbw(get_min_received_power(GPSL1C_D())) ≈ -163.0 atol = 1e-6
    @test _to_dbw(get_min_received_power(GPSL1C_P())) ≈ -158.25 atol = 1e-6
    @test _to_dbw(get_min_received_power(GPSL2CM())) ≈ -160.0 + 10log10(0.5) atol = 1e-6
    @test _to_dbw(get_min_received_power(GPSL2CL())) ≈ -160.0 + 10log10(0.5) atol = 1e-6
    @test _to_dbw(get_min_received_power(GalileoE1B())) ≈ -157.25 + 10log10(0.5) atol = 1e-6
    @test _to_dbw(get_min_received_power(GalileoE1C())) ≈ -157.25 + 10log10(0.5) atol = 1e-6
    @test _to_dbw(get_min_received_power(GalileoE5aI())) ≈ -155.25 + 10log10(0.5) atol = 1e-6
    @test _to_dbw(get_min_received_power(GalileoE5aQ())) ≈ -155.25 + 10log10(0.5) atol = 1e-6

    # Linear watts, type-stable, and folds identically on instance and type.
    @test get_min_received_power(GPSL1CA()) isa typeof(1.0u"W")
    @test @inferred(get_min_received_power(GPSL1CA())) == get_min_received_power(GPSL1CA)
    @test @inferred(get_min_received_power(GPSL1CA)) isa typeof(1.0u"W")

    # Component-split relationships.
    @test get_min_received_power(GalileoE1B()) == get_min_received_power(GalileoE1C())
    @test get_min_received_power(GPSL5I()) == get_min_received_power(GPSL5Q())
    # L1C pilot is 4.75 dB above the data component (IS-GPS-800J explicit
    # −158.25 vs −163.0; ≈ the 75/25 split, which is 10·log10(3) = 4.77 dB).
    @test _to_dbw(get_min_received_power(GPSL1C_P())) -
          _to_dbw(get_min_received_power(GPSL1C_D())) ≈ 4.75 atol = 1e-6

    # BOC(1,1) approximations inherit their parent's power.
    @test get_min_received_power(GalileoE1B_BOC11()) == get_min_received_power(GalileoE1B())
    @test get_min_received_power(GalileoE1C_BOC11()) == get_min_received_power(GalileoE1C())
end
