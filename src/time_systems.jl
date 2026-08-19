"""
    TimeSystem

Abstract supertype for a GNSS time scale.

A *time system* is a continuous atomic time scale (no leap seconds) that a
constellation references its broadcast time to: [`GPST`](@ref) for the GPS
signals, [`GST`](@ref) for the Galileo signals. Signals that report the same
time system share a receiver clock bias — that is the architectural reason
this abstraction exists, mirroring [`Band`](@ref) for the RF carrier.

A time scale is fixed by two constants, both queried through it:
[`get_system_start_time`](@ref) (its epoch) and [`get_tai_offset`](@ref) (its
offset from TAI).
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

"""
$(SIGNATURES)

Get the start epoch of a GNSS time scale, as a UTC `DateTime`.

This is the instant at which the broadcast time count (week number 0, time of
week 0) is zero. Works on a [`TimeSystem`](@ref), or on a signal / signal type
(which forwards through [`get_time_system`](@ref)).

The Galileo epoch is **not** on a UTC minute boundary: the Galileo OS SIS ICD
(Issue 2.2, §5.1.2) defines it as 13 s before midnight between 21 and 22
August 1999, i.e. `1999-08-21T23:59:47` UTC. The GPS epoch is
`1980-01-06T00:00:00` UTC (IS-GPS-200).

See also [`get_tai_offset`](@ref) for the time scale's constant offset from TAI.

# Examples
```julia-repl
julia> get_system_start_time(GPST())
1980-01-06T00:00:00

julia> get_system_start_time(GalileoE1B())
1999-08-21T23:59:47
```
"""
@inline get_system_start_time(::GPST) = DateTime(1980, 1, 6, 0, 0, 0)
@inline get_system_start_time(::GST) = DateTime(1999, 8, 21, 23, 59, 47)
@inline get_system_start_time(::BDT) = DateTime(2006, 1, 1, 0, 0, 0)

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
@inline get_system_start_time(::Type{S}) where {S<:AbstractGNSSSignal} =
    get_system_start_time(get_time_system(S))
@inline get_system_start_time(s::AbstractGNSSSignal) =
    get_system_start_time(get_time_system(s))
@inline get_tai_offset(::Type{S}) where {S<:AbstractGNSSSignal} =
    get_tai_offset(get_time_system(S))
@inline get_tai_offset(s::AbstractGNSSSignal) = get_tai_offset(get_time_system(s))

"""
$(SIGNATURES)

Get the `Symbol` identifier of a [`TimeSystem`](@ref) — `:GPST` or `:GST`.

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

Get the human-readable name of a [`TimeSystem`](@ref) — `"GPS Time"` or
`"Galileo System Time"`.

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
