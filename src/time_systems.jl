"""
    TimeSystem

Abstract supertype for a GNSS time scale.

A *time system* is a continuous atomic time scale (no leap seconds) that a
constellation references its broadcast time to: [`GPST`](@ref) for the GPS
signals, [`GST`](@ref) for the Galileo signals, [`BDT`](@ref) for the BeiDou
signals. Signals that report the same time system share a receiver clock bias —
that is the architectural reason this abstraction exists, mirroring
[`Band`](@ref) for the RF carrier.

A time scale is fixed by its epoch ([`get_week_start_epoch`](@ref)`(time_system, 0)`, on
the continuous TAI scale) and its offset from TAI ([`get_tai_offset`](@ref)). Between two scales, the difference of
their week starts is the elapsed seconds between them, and
[`get_time_system_offset`](@ref) the whole seconds between their counts within a week.

Those two offsets are Unitful seconds. Instants use [`TAITime`](@ref), whose arithmetic
accepts plain seconds or Unitful durations; their differences and times of week
([`get_time_of_week`](@ref)) are plain seconds, for full `Float64` resolution.
"""
abstract type TimeSystem end

"""
    GPST <: TimeSystem

GPS Time. Epoch `1980-01-06T00:00:00` UTC, `GPST = TAI − 19 s`, no leap
seconds (IS-GPS-200). Used by every GPS signal.
"""
struct GPST <: TimeSystem end

"""
    GST <: TimeSystem

Galileo System Time. Epoch is 13 s before midnight 21/22 August 1999, i.e.
`1999-08-21T23:59:47` UTC; `GST = TAI − 19 s`, no leap seconds (Galileo OS SIS
ICD, Issue 2.2, §5.1.2). Used by every Galileo signal.
"""
struct GST <: TimeSystem end

"""
    BDT <: TimeSystem

BeiDou Time. Epoch `2006-01-01T00:00:00` UTC, `BDT = TAI − 33 s`, no leap
seconds (BDS-SIS-ICD-B1I-3.0 §5.2 / BDS-SIS-ICD, "BDT is synchronised to UTC
within 100 ns"). Used by every BeiDou signal. Because its epoch is 26 years
later than GPS's while both are `TAI − 19 s`/`TAI − 33 s`, `BDT = GPST − 14 s`
in elapsed time (14 leap seconds accrued between the two epochs).
"""
struct BDT <: TimeSystem end

# Prints what Base's default does, but without its reflective fallback: a message that
# interpolates a time system (as `Dictionary`'s missing-key error does with its key)
# would otherwise reach a dynamic `show` that `juliac --trim` rejects.
Base.show(io::IO, x::TimeSystem) = print(io, typeof(x), "()")

"""
$(SIGNATURES)

Get the constant offset between a GNSS time scale and TAI, as a Unitful time
quantity.

The value is the leap-second-free amount by which TAI leads the time scale:

    TAI = system_time + get_tai_offset(time_system)

equivalently `system_time = TAI − get_tai_offset(time_system)`. Works on a
[`TimeSystem`](@ref), or on a signal / signal type (which forwards through
[`get_time_system`](@ref)).

The two GPS/Galileo time scales — [`GPST`](@ref) and [`GST`](@ref) — are
defined as `TAI − 19 s`, so this returns `19s` for either; [`BDT`](@ref)
(BeiDou Time) is `TAI − 33 s`, so it returns `33s`. It is defined per time
system (not as a shared fallback) so each states its own value.

# Examples

```julia-repl
julia> get_tai_offset(GST())
19 s

julia> get_tai_offset(BDT())
33 s

julia> get_tai_offset(GPSL1CA) == get_tai_offset(GalileoE1B)
true
```
"""
@inline get_tai_offset(::GPST) = 19s
@inline get_tai_offset(::GST) = 19s
@inline get_tai_offset(::BDT) = 33s

"""
$(SIGNATURES)

Get the [`TimeSystem`](@ref) a signal's measurements are referenced to.

The time system is a per-constellation fact, so it is defined once per
constellation on the signal supertype, e.g.
`get_time_system(::Type{<:AbstractGPSSignal}) = GPST()`; every concrete GPS
signal inherits it through subtype dispatch. Available on the type (without
constructing a signal); the instance method forwards to the type.

# Examples

```julia-repl
julia> get_time_system(GPSL1CA)
GPST()

julia> get_time_system(GalileoE1B())
GST()
```
"""
function get_time_system end

@inline get_time_system(s::AbstractGNSSSignal) = get_time_system(typeof(s))

# Signal → time system. This is a per-constellation fact, so it is stated once
# per constellation through the signal supertype (unlike the genuinely
# per-signal `get_band`); the epoch and TAI offset are properties of the time
# system, defined once above.
@inline get_time_system(::Type{<:AbstractGPSSignal}) = GPST()         # GPS Time
@inline get_time_system(::Type{<:AbstractGalileoSignal}) = GST()      # Galileo System Time
@inline get_time_system(::Type{<:AbstractBeiDouSignal}) = BDT()       # BeiDou Time

# Signal-level access forwards through the time system — same dual entry (type
# or instance) as `get_center_frequency` through `get_band`.
@inline get_tai_offset(::Type{S}) where {S<:AbstractGNSSSignal} =
    get_tai_offset(get_time_system(S))
@inline get_tai_offset(s::AbstractGNSSSignal) = get_tai_offset(get_time_system(s))

"""
    SECONDS_PER_WEEK

The length of a GNSS week in seconds, `604800`. Every time system here counts its
broadcast time as a week number plus a time of week of this length.

Not exported: bind it with `using GNSSSignals: SECONDS_PER_WEEK`.
"""
const SECONDS_PER_WEEK = 604_800

"""
    fold_week_crossover(Δt) -> seconds

Fold a difference `Δt` of two seconds-of-week counts into `[-302400, 302400]`: across a
week crossover (IS-GPS-200's term, §20.3.3.3.3.1) the raw difference is off by ±604800 s.
The result has the type of `Δt`, so integer seconds stay integer.
A difference of two counts within half a week of each other passes through this; one
meant to span whole weeks must not (a broadcast time offset adds the weeks back itself).

Not exported, as [`SECONDS_PER_WEEK`](@ref GNSSSignals.SECONDS_PER_WEEK) is not: bind it
with `using GNSSSignals: fold_week_crossover`.

```julia-repl
julia> GNSSSignals.fold_week_crossover(604799.5 - 0.5)
-1.0
```
"""
function fold_week_crossover(Δt)
    # Integer week arithmetic, so an integer difference stays an integer and a
    # floating-point one keeps its type.
    half_week = SECONDS_PER_WEEK ÷ 2
    # The week in `Δt`'s own type, so the result has that type whichever branch it takes.
    Δt > half_week ? Δt - oftype(Δt, SECONDS_PER_WEEK) :
    Δt < -half_week ? Δt + oftype(Δt, SECONDS_PER_WEEK) : Δt
end

"""
    get_time_system_offset(from::TimeSystem, to::TimeSystem) -> Unitful seconds

The seconds to add to a time of week counted in `from` to read the same instant in
the count of `to`: the fixed difference between the definitions of the two scales —
not the leap seconds between their epochs, of which 13 lie between the GPST and GST
epochs while the offset is `0 s`,

    get_time_system_offset(from, to) = get_tai_offset(from) - get_tai_offset(to)

`0 s` between [`GPST`](@ref) and [`GST`](@ref), both `TAI − 19 s`; a [`BDT`](@ref)
reading is 14 s behind, so `get_time_system_offset(BDT(), GPST()) == 14 s`. It is the
difference of the two scales' week starts ([`get_week_start_epoch`](@ref)) without its
whole weeks, so it moves a time of week, not a whole count, and a time of week moved by it
can cross a week boundary ([`GNSSSignals.fold_week_crossover`](@ref) folds differences of
such times). To convert a whole `(WN, TOW)`, go through TAI instead:
`get_week_and_time_of_week(to, get_week_start_epoch(from, WN) + TOW)` already includes
this offset, so do not add it again. Antisymmetric; zero from a scale to itself.

It follows from the scales' definitions and is never broadcast: a broadcast time
offset (GGTO, BGTO) carries only the nanosecond steering residual on top of it
(IS-GPS-800J §3.5.4.2.1.1: users "must also apply any integer seconds difference
between the systems using definitions of each system time scale").

# Examples

```julia-repl
julia> get_time_system_offset(BDT(), GPST())
14 s

julia> get_time_system_offset(GST(), GPST())
0 s
```
"""
@inline get_time_system_offset(from::TimeSystem, to::TimeSystem) =
    get_tai_offset(from) - get_tai_offset(to)

"""
    get_week_start_epoch(time_system, week) -> TAITime

The TAI epoch at which `week` of a GNSS time system begins: the system's origin (week 0)
plus `week` whole weeks, exactly. A time of week `t` counted in that system dates as
`get_week_start_epoch(time_system, week) + t`, with full `Float64` resolution in the
fraction. TAI is continuous, so — unlike a difference of UTC labels, which misses the leap
seconds between them — the difference of two of these is an elapsed duration, across
scales too: it is what converts a week and time of week from one scale to another. No
leap-second table is consulted, so it is safe at precompile time. Works on a [`TimeSystem`](@ref), or on a
signal / signal type (which forwards through [`get_time_system`](@ref)).

The ICDs state the origins as UTC labels: `1980-01-06T00:00:00` for GPST (IS-GPS-200),
`1999-08-21T23:59:47` for GST, 13 s before midnight and so not on a UTC minute (Galileo OS
SIS ICD, Issue 2.2, §5.1.2), and `2006-01-01T00:00:00` for BDT. They are stated here on
TAI rather than derived from those labels, because deriving them would need a
leap-second table, and because GST's is not derivable from this package's other
constants at all: the Galileo ICD anchors the GST epoch to GPS Time at the week-1024
rollover — `GST − UTC` was already 13 s at the epoch — so its origin is
`1999-08-21T23:59:47 UTC + 32 s (TAI − UTC in 1999)`, while
`get_tai_offset(GST()) == 19 s` reflects none of that. GPST and BDT were each set equal to
UTC at their epochs, so their origins are simply the UTC label plus their own
[`get_tai_offset`](@ref) — a coincidence GST breaks, which is why the naive
`UTC difference + TAI-offset difference` arithmetic is 13 s wrong exactly and only for
Galileo. The origins' differences are whole weeks, or whole weeks and the 14 s of
[`get_time_system_offset`](@ref): GST week 0 is GPS week 1024, BDT week 0 is GPS week 1356
and 14 s.

```julia-repl
julia> get_week_start_epoch(GPST(), 0)
1980-01-06T00:00:19 TAI

julia> get_week_start_epoch(GST(), 0)
1999-08-22T00:00:19 TAI

julia> get_week_start_epoch(BDT(), 0)
2006-01-01T00:00:33 TAI

julia> get_week_start_epoch(GPST(), 2364)
2025-04-27T00:00:19 TAI

julia> get_week_start_epoch(BDT(), 1008) - get_week_start_epoch(GPST(), 2364)
14.0
```
"""
@inline get_week_start_epoch(time_system::TimeSystem, week::Integer) =
    _tai_origin(time_system) + week * SECONDS_PER_WEEK

# The scales' origins on TAI; why they are stated rather than derived is in
# `get_week_start_epoch`'s docstring.
@inline _tai_origin(::GPST) = TAITime(DateTime(1980, 1, 6, 0, 0, 19))
@inline _tai_origin(::GST) = TAITime(DateTime(1999, 8, 22, 0, 0, 19))
@inline _tai_origin(::BDT) = TAITime(DateTime(2006, 1, 1, 0, 0, 33))
@inline get_week_start_epoch(::Type{S}, week::Integer) where {S<:AbstractGNSSSignal} =
    get_week_start_epoch(get_time_system(S), week)
@inline get_week_start_epoch(s::AbstractGNSSSignal, week::Integer) =
    get_week_start_epoch(get_time_system(s), week)

"""
    get_week_and_time_of_week(time_system, epoch::TAITime) -> (week, seconds)

The week and time of week (s) of the TAI `epoch` on the count of a GNSS time system,
together: [`get_week`](@ref) and [`get_time_of_week`](@ref) in one call, the inverse of
[`get_week_start_epoch`](@ref). The time of week is in `[0, 604800]`, as
[`get_time_of_week`](@ref) explains.

Works on a [`TimeSystem`](@ref), or on a signal / signal type (which forwards through
[`get_time_system`](@ref)).

```julia-repl
julia> get_week_and_time_of_week(BDT(), get_week_start_epoch(GPST(), 2364) + 3.5)
(1007, 604789.5)
```
"""
function get_week_and_time_of_week(time_system::TimeSystem, epoch::TAITime)
    # The whole weeks come off in integer seconds before any floating point, so the time
    # of week is good to half a `Float64` step at up to 604800 s, about 6e-11 s (an
    # absolute count of about 1.2e9 s would leave only 0.2 µs). The origins fall on whole
    # TAI seconds, so the fraction is the epoch's own.
    origin = get_week_start_epoch(time_system, 0)
    elapsed_seconds = epoch.second - origin.second
    week = fld(elapsed_seconds, SECONDS_PER_WEEK)
    # An epoch in the last ~6e-11 s of a week rounds up to exactly 604800.0 here, the same
    # instant as week + 1 at 0.0. Left as is: every consumer reduces a time of week modulo
    # a day or the week, or folds differences of it (`fold_week_crossover`), all of which
    # treat 604800.0 as 0.0, so carrying it would only serve the stated range.
    time_of_week = Float64(
        (elapsed_seconds - week * SECONDS_PER_WEEK) + (epoch.fraction - origin.fraction),
    )
    Int(week), time_of_week
end

@inline get_week_and_time_of_week(::Type{S}, epoch::TAITime) where {S<:AbstractGNSSSignal} =
    get_week_and_time_of_week(get_time_system(S), epoch)
@inline get_week_and_time_of_week(s::AbstractGNSSSignal, epoch::TAITime) =
    get_week_and_time_of_week(get_time_system(s), epoch)

"""
    get_week(time_system, epoch::TAITime) -> Int

The week, counted from the system's origin, that the TAI `epoch` falls in on the count
of a GNSS time system; with [`get_time_of_week`](@ref) the inverse of
[`get_week_start_epoch`](@ref):

    get_week_start_epoch(time_system, get_week(time_system, epoch)) +
        get_time_of_week(time_system, epoch) ≈ epoch

to within about 6e-11 s, the resolution of a `Float64` time of week near the end of a
week, a sub-millimetre range; it is exact where that resolution holds the fraction.

Works on a [`TimeSystem`](@ref), or on a signal / signal type (which forwards through
[`get_time_system`](@ref)).

```julia-repl
julia> get_week(GPST(), get_week_start_epoch(GPST(), 2364) + 3.5)
2364

julia> get_week(BDT(), get_week_start_epoch(GPST(), 2364) + 3.5)
1007
```
"""
@inline get_week(time_system::TimeSystem, epoch::TAITime) =
    first(get_week_and_time_of_week(time_system, epoch))
@inline get_week(::Type{S}, epoch::TAITime) where {S<:AbstractGNSSSignal} =
    get_week(get_time_system(S), epoch)
@inline get_week(s::AbstractGNSSSignal, epoch::TAITime) = get_week(get_time_system(s), epoch)

"""
    get_time_of_week(time_system, epoch::TAITime) -> Float64

The time of week (s), in `[0, 604800]`, of the TAI `epoch` on the count of a GNSS time
system; with [`get_week`](@ref) the inverse of [`get_week_start_epoch`](@ref). The whole
weeks come off in integer seconds before any floating point, so the time of week is good
to half a `Float64` step at its magnitude, at most about 6e-11 s, where the absolute count
since the origin would leave only 0.2 µs. An epoch in the last ~6e-11 s of a week rounds
up to 604800.0, the same instant as the next week's 0.0.

Works on a [`TimeSystem`](@ref), or on a signal / signal type (which forwards through
[`get_time_system`](@ref)).

```julia-repl
julia> get_time_of_week(GPST(), get_week_start_epoch(GPST(), 2364) + 3.5)
3.5

julia> get_time_of_week(BDT(), get_week_start_epoch(GPST(), 2364) + 3.5)
604789.5
```
"""
@inline get_time_of_week(time_system::TimeSystem, epoch::TAITime) =
    last(get_week_and_time_of_week(time_system, epoch))
@inline get_time_of_week(::Type{S}, epoch::TAITime) where {S<:AbstractGNSSSignal} =
    get_time_of_week(get_time_system(S), epoch)
@inline get_time_of_week(s::AbstractGNSSSignal, epoch::TAITime) =
    get_time_of_week(get_time_system(s), epoch)

"""
$(SIGNATURES)

Get the `Symbol` identifier of a [`TimeSystem`](@ref) — `:GPST`, `:GST` or `:BDT`.

The machine-readable key for a time scale: the level at which signals share a
receiver clock bias, so every signal referenced to the same scale maps to the
same id. Its counterpart along the identity axis is [`get_constellation_id`](@ref)
— currently one-to-one with it, but a distinct fact (a constellation can
broadcast a scale it does not own).

Defaults to `nameof` of the time system type, so a new [`TimeSystem`](@ref) gets
a sensible id for free; override `get_time_system_id(::Type{MySystem})` if you
need a different symbol. Works on a time system or signal, instance or type, and
folds to a compile-time constant.

```julia-repl
julia> get_time_system_id(GPST())
:GPST

julia> get_time_system_id(GalileoE1B)
:GST
```
"""
@inline get_time_system_id(::Type{T}) where {T<:TimeSystem} = nameof(T)
@inline get_time_system_id(t::TimeSystem) = get_time_system_id(typeof(t))
@inline get_time_system_id(::Type{S}) where {S<:AbstractGNSSSignal} =
    get_time_system_id(get_time_system(S))
@inline get_time_system_id(s::AbstractGNSSSignal) = get_time_system_id(get_time_system(s))

"""
$(SIGNATURES)

Get the human-readable name of a [`TimeSystem`](@ref) — `"GPS Time"`,
`"Galileo System Time"` or `"BeiDou Time"`.

The display counterpart to [`get_time_system_id`](@ref) — same granularity, but
a `String` meant for log lines and user-facing output rather than a key to branch
or dictionary on (prefer the `Symbol` id for that).

Stated per time system rather than derived from the id, because the spelled-out
names are the ICDs' own (`"Galileo System Time"`, Galileo OS SIS ICD §5.1.2) and
bear no relation to the acronyms. Works on a time system or signal, instance or
type, and folds to a compile-time constant.

```julia-repl
julia> get_time_system_name(GST())
"Galileo System Time"

julia> get_time_system_name(GPSL1CA)
"GPS Time"
```
"""
@inline get_time_system_name(::Type{GPST}) = "GPS Time"
@inline get_time_system_name(::Type{GST}) = "Galileo System Time"
@inline get_time_system_name(::Type{BDT}) = "BeiDou Time"
@inline get_time_system_name(t::TimeSystem) = get_time_system_name(typeof(t))
@inline get_time_system_name(::Type{S}) where {S<:AbstractGNSSSignal} =
    get_time_system_name(get_time_system(S))
@inline get_time_system_name(s::AbstractGNSSSignal) =
    get_time_system_name(get_time_system(s))
