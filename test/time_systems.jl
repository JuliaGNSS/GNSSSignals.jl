using Dates: DateTime
import Unitful: s

@testset "Time systems" begin
    gps_signals = (GPSL1CA, GPSL1C_D, GPSL1C_P, GPSL2CM, GPSL2CL, GPSL5I, GPSL5Q)
    galileo_signals =
        (GalileoE1B, GalileoE1B_BOC11, GalileoE1C, GalileoE1C_BOC11, GalileoE5aI, GalileoE5aQ)

    @test GPST <: TimeSystem
    @test GST <: TimeSystem

    # Every signal belongs to its constellation supertype; both supertypes are
    # signal subtypes. This is what lets `get_time_system` be stated once per
    # constellation.
    @test AbstractGPSSignal <: AbstractGNSSSignal
    @test AbstractGalileoSignal <: AbstractGNSSSignal
    for S in gps_signals
        @test S <: AbstractGPSSignal
    end
    for S in galileo_signals
        @test S <: AbstractGalileoSignal
    end

    # Signal → time system (a per-constellation fact, dispatched via the
    # constellation supertype).
    for S in gps_signals
        @test @inferred(get_time_system(S)) === GPST()
    end
    for S in galileo_signals
        @test @inferred(get_time_system(S)) === GST()
    end
    # Instance path forwards to the type.
    @test get_time_system(GPSL1CA()) === get_time_system(GPSL1CA) === GPST()
    @test get_time_system(GalileoE1B()) === get_time_system(GalileoE1B) === GST()

    # Epoch / offset defined on the time system (IS-GPS-200 / Galileo OS SIS
    # ICD 2.2 §5.1.2). The Galileo epoch is not a UTC minute boundary.
    @test @inferred(get_system_start_time(GPST())) == DateTime(1980, 1, 6, 0, 0, 0)
    @test @inferred(get_system_start_time(GST())) == DateTime(1999, 8, 21, 23, 59, 47)
    @test @inferred(get_tai_offset(GPST())) == 19s
    @test @inferred(get_tai_offset(GST())) == 19s

    # Signal-level access forwards through the time system, on type or instance.
    for S in gps_signals
        @test get_system_start_time(S) == DateTime(1980, 1, 6, 0, 0, 0)
        @test get_system_start_time(S()) == get_system_start_time(S)
        @test get_tai_offset(S) == 19s
    end
    for S in galileo_signals
        @test get_system_start_time(S) == DateTime(1999, 8, 21, 23, 59, 47)
        @test get_system_start_time(S()) == get_system_start_time(S)
        @test get_tai_offset(S) == 19s
    end

    # The TAI offset is identical across every implemented signal — both GPST
    # and GST are defined as TAI − 19 s.
    all_signals = (gps_signals..., galileo_signals...)
    @test length(unique(get_tai_offset(S) for S in all_signals)) == 1

    # GPS and Galileo share the offset but differ in epoch / time system.
    @test get_tai_offset(GPSL1CA) == get_tai_offset(GalileoE1B)
    @test get_system_start_time(GPSL1CA) != get_system_start_time(GalileoE1B)
    @test get_time_system(GPSL1CA) !== get_time_system(GalileoE1B)
end
