using Dates: Dates, DateTime
import Unitful
import Unitful: s, Time

@testset "Time systems" begin
    gps_signals = (GPSL1CA, GPSL1C_D, GPSL1C_P, GPSL2CM, GPSL2CL, GPSL5I, GPSL5Q)
    galileo_signals =
        (GalileoE1B, GalileoE1B_BOC11, GalileoE1C, GalileoE1C_BOC11, GalileoE5aI, GalileoE5aQ)
    beidou_signals =
        (BeiDouB1I, BeiDouB3I, BeiDouB2bI, BeiDouB2aI, BeiDouB2aQ, BeiDouB1C_D, BeiDouB1C_P)

    @test GPST <: TimeSystem
    @test GST <: TimeSystem
    @test BDT <: TimeSystem

    # Every signal belongs to its constellation supertype; both supertypes are
    # signal subtypes. This is what lets `get_time_system` be stated once per
    # constellation.
    @test AbstractGPSSignal <: AbstractGNSSSignal
    @test AbstractGalileoSignal <: AbstractGNSSSignal
    @test AbstractBeiDouSignal <: AbstractGNSSSignal
    for S in gps_signals
        @test S <: AbstractGPSSignal
    end
    for S in galileo_signals
        @test S <: AbstractGalileoSignal
    end
    for S in beidou_signals
        @test S <: AbstractBeiDouSignal
    end

    # Signal → time system (a per-constellation fact, dispatched via the
    # constellation supertype).
    for S in gps_signals
        @test @inferred(get_time_system(S)) === GPST()
    end
    for S in galileo_signals
        @test @inferred(get_time_system(S)) === GST()
    end
    for S in beidou_signals
        @test @inferred(get_time_system(S)) === BDT()
    end
    # Instance path forwards to the type.
    @test get_time_system(GPSL1CA()) === get_time_system(GPSL1CA) === GPST()
    @test get_time_system(GalileoE1B()) === get_time_system(GalileoE1B) === GST()
    @test get_time_system(BeiDouB1I()) === get_time_system(BeiDouB1I) === BDT()

    # Epoch / offset defined on the time system (IS-GPS-200 / Galileo OS SIS
    # ICD 2.2 §5.1.2). The Galileo epoch is not a UTC minute boundary.
    @test @inferred(get_system_start_time(GPST())) == DateTime(1980, 1, 6, 0, 0, 0)
    @test @inferred(get_system_start_time(GST())) == DateTime(1999, 8, 21, 23, 59, 47)
    @test @inferred(get_system_start_time(BDT())) == DateTime(2006, 1, 1, 0, 0, 0)
    @test @inferred(get_tai_offset(GPST())) == 19s
    @test @inferred(get_tai_offset(GST())) == 19s
    @test @inferred(get_tai_offset(BDT())) == 33s

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
    for S in beidou_signals
        @test get_system_start_time(S) == DateTime(2006, 1, 1, 0, 0, 0)
        @test get_system_start_time(S()) == get_system_start_time(S)
        @test get_tai_offset(S) == 33s
    end

    # GPS and Galileo share the TAI − 19 s offset; BeiDou uses TAI − 33 s.
    gps_gal = (gps_signals..., galileo_signals...)
    @test length(unique(get_tai_offset(S) for S in gps_gal)) == 1

    # GPS and Galileo share the offset but differ in epoch / time system.
    @test get_tai_offset(GPSL1CA) == get_tai_offset(GalileoE1B)
    @test get_system_start_time(GPSL1CA) != get_system_start_time(GalileoE1B)
    @test get_time_system(GPSL1CA) !== get_time_system(GalileoE1B)

    # BeiDou differs from both GPS and Galileo in offset, epoch and time system.
    @test get_tai_offset(BeiDouB1I) != get_tai_offset(GPSL1CA)
    @test get_time_system(BeiDouB1I) !== get_time_system(GalileoE1B)
end

@testset "get_time_system_id / get_time_system_name" begin
    # On the time system itself (instance and type).
    @test @inferred(get_time_system_id(GPST())) === :GPST
    @test @inferred(get_time_system_id(GST())) === :GST
    @test @inferred(get_time_system_id(GST)) === :GST
    @test @inferred(get_time_system_id(BDT())) === :BDT
    @test @inferred(get_time_system_name(GPST())) == "GPS Time"
    @test @inferred(get_time_system_name(GST())) == "Galileo System Time"
    @test @inferred(get_time_system_name(BDT())) == "BeiDou Time"
    @test @inferred(get_time_system_name(GPST)) == "GPS Time"

    # On a signal: forwards through get_time_system, on type or instance.
    for S in (GPSL1CA, GPSL2CM, GPSL5Q)
        @test @inferred(get_time_system_id(S)) === :GPST
        @test @inferred(get_time_system_name(S)) == "GPS Time"
        @test get_time_system_id(S()) === get_time_system_id(S)
        @test get_time_system_name(S()) == get_time_system_name(S)
    end
    for S in (GalileoE1B, GalileoE1C, GalileoE5aI)
        @test @inferred(get_time_system_id(S)) === :GST
        @test @inferred(get_time_system_name(S)) == "Galileo System Time"
        @test get_time_system_id(S()) === get_time_system_id(S)
        @test get_time_system_name(S()) == get_time_system_name(S)
    end
    for S in (BeiDouB1I, BeiDouB3I, BeiDouB2aQ, BeiDouB1C_D)
        @test @inferred(get_time_system_id(S)) === :BDT
        @test @inferred(get_time_system_name(S)) == "BeiDou Time"
        @test get_time_system_id(S()) === get_time_system_id(S)
        @test get_time_system_name(S()) == get_time_system_name(S)
    end

    # Every signal referenced to one scale shares id and name, across bands.
    @test get_time_system_id(GPSL1CA) === get_time_system_id(GPSL5I)
    @test get_time_system_id(GPSL1CA) !== get_time_system_id(GalileoE1B)
    # The name is the ICD's spelled-out one, not the acronym.
    @test get_time_system_name(GalileoE1B) != String(get_time_system_id(GalileoE1B))
end

@testset "every time system answers every time-system accessor" begin
    for T in ALL_TIME_SYSTEMS
        @test get_time_system_id(T) isa Symbol
        # Deliberately the one name with no fallback to its id: the ICDs' spelled-out
        # names bear no relation to the acronyms, so each system states its own and a
        # new one has to say it here rather than report "BDT" and call it a name.
        @test get_time_system_name(T) isa String
        @test get_time_system_name(T) != String(get_time_system_id(T))
        # The two constants a time scale is fixed by.
        @test get_system_start_time(T()) isa DateTime
        @test get_tai_offset(T()) isa Time
    end
end

@testset "TAI epochs make cross-scale arithmetic exact" begin
    # The values themselves: UTC label + (TAI − UTC at the epoch). GPST and BDT
    # were set equal to UTC at their epochs, so their TAI labels are the UTC
    # label plus their own TAI offset; GST's is not — its epoch is anchored to
    # GPS Time at the week-1024 rollover (OS SIS ICD §5.1.2), 13 s after the
    # naive derivation, which is the whole reason this accessor exists.
    @test get_tai_system_start_time(GPST()) == DateTime(1980, 1, 6, 0, 0, 19)
    @test get_tai_system_start_time(GST()) == DateTime(1999, 8, 22, 0, 0, 19)
    @test get_tai_system_start_time(BDT()) == DateTime(2006, 1, 1, 0, 0, 33)
    for T in ALL_TIME_SYSTEMS
        naive =
            get_system_start_time(T()) +
            Dates.Second(Int(Unitful.ustrip(s, get_tai_offset(T()))))
        if T === GST
            @test get_tai_system_start_time(T()) == naive + Dates.Second(13)
        else
            @test get_tai_system_start_time(T()) == naive
        end
        @test @inferred(get_tai_system_start_time(T())) isa DateTime
    end
    # Signal forwarding, type and instance, like the other epoch accessors.
    @test @inferred(get_tai_system_start_time(GalileoE1B)) ==
          get_tai_system_start_time(GST())
    @test get_tai_system_start_time(GPSL1CA()) == get_tai_system_start_time(GPST())

    # The one number cross-scale (WN, TOW) conversion needs. GST week 0 IS GPS
    # week 1024 (the rollover instant); BDT week 0 begins 14 s into GPS week
    # 1356 — the same 14 s that separates the two scales' counts.
    @test get_epoch_offset(GST(), GPST()) == 1024 * 604800 * s
    @test get_epoch_offset(BDT(), GPST()) == (1356 * 604800 + 14) * s
    @test get_epoch_offset(BDT(), GST()) == (332 * 604800 + 14) * s
    for A in ALL_TIME_SYSTEMS, B in ALL_TIME_SYSTEMS
        @test @inferred(get_epoch_offset(A(), B())) == -get_epoch_offset(B(), A())
        A === B && @test get_epoch_offset(A(), B()) == 0s
    end
end
