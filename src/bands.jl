"""
    Band

Abstract supertype for GNSS RF bands.

A *band* represents a shared RF carrier frequency. Signals that report the
same band can share a carrier NCO in a receiver — that is the architectural
reason this abstraction exists.

Band identity here is by RF frequency, not by ICD label: GPS L1 and Galileo E1
both report [`L1`](@ref) because they share 1575.42 MHz. If you need the
ICD-specific name, look at the concrete signal instead (e.g. via
[`get_signal_name`](@ref)).
"""
abstract type Band end

"""
    L1 <: Band

The 1575.42 MHz GNSS band. Shared by GPS L1 C/A, GPS L1C (data and pilot),
Galileo E1 (B and C), BeiDou B1C, QZSS L1C, and others.
"""
struct L1 <: Band end

"""
    L2 <: Band

The 1227.6 MHz GNSS band. Carries the GPS L2 civil signals (L2C: the L2 CM
and L2 CL codes) and the legacy L2 P(Y) signal.
"""
struct L2 <: Band end

"""
    L5 <: Band

The 1176.45 MHz GNSS band. Shared by GPS L5 (I and Q), Galileo E5a, and
BeiDou B2a.
"""
struct L5 <: Band end

"""
    B1I <: Band

The 1561.098 MHz GNSS band (the legacy BeiDou B1 frequency). Carries the
BeiDou B1I signal. Distinct from the [`L1`](@ref) band (1575.42 MHz) that
carries BeiDou B1C.
"""
struct B1I <: Band end

"""
    B3I <: Band

The 1268.52 MHz GNSS band (the BeiDou B3 frequency). Carries the BeiDou B3I
signal.
"""
struct B3I <: Band end

"""
    E5b <: Band

The 1207.14 MHz GNSS band. Carries the Galileo E5b signals (E5b-I / E5b-Q) and
the BeiDou B2b signal — the two constellations share the carrier exactly: 118 ×
10.23 MHz, named E5b by the Galileo OS SIS ICD (Table 2) and B2b (formerly B2I)
by the BeiDou ICDs. One of the two labels has to name the band, and this package
uses the Galileo one, as it uses the GPS ones for [`L1`](@ref) and [`L5`](@ref).
"""
struct E5b <: Band end

"""
    E6 <: Band

The 1278.75 MHz GNSS band (the Galileo E6 frequency). Carries the Galileo E6-B
and E6-C signals — and with them the C/NAV message and the High Accuracy
Service — as well as the QZSS L6 signals, which share the carrier.
"""
struct E6 <: Band end

"""
$(SIGNATURES)

Get the center (carrier) frequency of a band.

One method per concrete [`Band`](@ref) — `L1` 1575.42 MHz, `L2` 1227.6 MHz,
`L5` 1176.45 MHz, `B1I` 1561.098 MHz, `B3I` 1268.52 MHz, `E5b` 1207.14 MHz,
`E6` 1278.75 MHz.
Works on a band instance or its type, as [`get_band_id`](@ref) /
[`get_band_name`](@ref) do.

```julia-repl
julia> get_center_frequency(L1())
1575420000 Hz

julia> get_center_frequency(L5)
1176450000 Hz
```
"""
@inline get_center_frequency(::Type{L1}) = 1_575_420_000Hz
@inline get_center_frequency(::Type{L2}) = 1_227_600_000Hz
@inline get_center_frequency(::Type{L5}) = 1_176_450_000Hz
@inline get_center_frequency(::Type{B1I}) = 1_561_098_000Hz
@inline get_center_frequency(::Type{B3I}) = 1_268_520_000Hz
@inline get_center_frequency(::Type{E5b}) = 1_207_140_000Hz
@inline get_center_frequency(::Type{E6}) = 1_278_750_000Hz
@inline get_center_frequency(b::Band) = get_center_frequency(typeof(b))

"""
$(SIGNATURES)

Get the center (carrier) frequency of a signal.

Dispatches through the signal's band ([`get_band`](@ref)), so all signals on
the same band return the same value by construction. Works on either a signal
instance or its type — `get_center_frequency(GPSL1CA)` avoids constructing a
signal just to read the carrier.

```julia-repl
julia> get_center_frequency(GPSL1CA())
1575420000 Hz

julia> get_center_frequency(GPSL1CA)
1575420000 Hz
```
"""
@inline get_center_frequency(::Type{S}) where {S<:AbstractGNSSSignal} =
    get_center_frequency(get_band(S))
@inline get_center_frequency(s::AbstractGNSSSignal) = get_center_frequency(get_band(s))

"""
$(SIGNATURES)

Get the [`Band`](@ref) a signal is transmitted on.

Concrete signal types define one method each on the type, e.g.
`get_band(::Type{<:GPSL1CA}) = L1()`, so the band is available without
constructing a signal. The instance method forwards to the type.

Band identity is by RF frequency, so signals of different constellations sharing
a carrier report the same band — `get_band(GalileoE1B())` is `L1()`. See
[`get_band_id`](@ref) / [`get_band_name`](@ref) for the band's key and label, and
[`get_center_frequency`](@ref) for its carrier frequency.

# Examples
```julia-repl
julia> get_band(GPSL1CA())
L1()

julia> get_band(GalileoE5aI)
L5()
```
"""
function get_band end

@inline get_band(s::AbstractGNSSSignal) = get_band(typeof(s))

"""
$(SIGNATURES)

Get the `Symbol` identifier of a [`Band`](@ref).

This is the machine-readable key for a band (e.g. `:L1`, `:L5`) — the level at
which signals share a carrier NCO and a receiver shares an inter-frequency bias.
Because band identity is by RF frequency, every signal on the same carrier maps
to the same id regardless of constellation: GPS L1 C/A, GPS L1C and Galileo E1
are all `:L1`.

Defaults to `nameof` of the band type, so a new [`Band`](@ref) gets a sensible
id for free; override `get_band_id(::Type{MyBand})` if you need a different
symbol. Works on a band or signal, instance or type, and dispatches to a
compile-time constant.

Distinct from [`get_signal_id`](@ref) (per-signal, e.g. `:GPSL1CA`) and from
[`get_center_frequency`](@ref) (the band's numeric carrier frequency).

```julia-repl
julia> get_band_id(L1())
:L1

julia> get_band_id(GPSL1CA())
:L1
```
"""
@inline get_band_id(::Type{B}) where {B<:Band} = nameof(B)
@inline get_band_id(b::Band) = get_band_id(typeof(b))
@inline get_band_id(::Type{S}) where {S<:AbstractGNSSSignal} = get_band_id(get_band(S))
@inline get_band_id(s::AbstractGNSSSignal) = get_band_id(get_band(s))

"""
$(SIGNATURES)

Get the human-readable name of a [`Band`](@ref), e.g. `"L1"`, `"L5"`.

The display counterpart to [`get_band_id`](@ref) — same granularity, but a
`String` meant for log lines and user-facing output rather than a key to
branch or dictionary on (prefer the `Symbol` id for that; comparing symbols is
cheaper and the name is free to be respelled).

Band identity is by RF frequency, so the name is the *band's* label, not the
constellation's label for it: `get_band_name(GalileoE1B())` is `"L1"`, not
`"E1"`. Use [`get_signal_name`](@ref) (`"Galileo E1B"`) when you need the
ICD-specific naming.

Stated as a literal per band — `get_band_name(::Type{L1}) = "L1"` — so naming a
band is free: the call allocates nothing and folds to a compile-time constant.
Each literal still spells exactly what `String` of [`get_band_id`](@ref) yields,
and a test pins that, so the two identity layers cannot drift apart.

A [`Band`](@ref) declared elsewhere needs no method of its own: the
`String(get_band_id(B))` fallback names it for free. That path allocates a fresh
`String` on every call (32 bytes, measured) even where the band type is
statically known, because a heap-allocated `String` — unlike an interned `Symbol`
or a literal already in the method's source — cannot be baked into the IR as a
constant. State the name as a literal, as the bands here do, if you are naming
such a band in a hot loop rather than a log line.

Works on a band or signal, instance or type.

```julia-repl
julia> get_band_name(L1())
"L1"

julia> get_band_name(GPSL1CA())
"L1"
```
"""
@inline get_band_name(::Type{L1}) = "L1"
@inline get_band_name(::Type{L2}) = "L2"
@inline get_band_name(::Type{L5}) = "L5"
@inline get_band_name(::Type{B1I}) = "B1I"
@inline get_band_name(::Type{B3I}) = "B3I"
@inline get_band_name(::Type{E5b}) = "E5b"
@inline get_band_name(::Type{E6}) = "E6"
@inline get_band_name(::Type{B}) where {B<:Band} = String(get_band_id(B))
@inline get_band_name(b::Band) = get_band_name(typeof(b))
@inline get_band_name(::Type{S}) where {S<:AbstractGNSSSignal} = get_band_name(get_band(S))
@inline get_band_name(s::AbstractGNSSSignal) = get_band_name(get_band(s))
