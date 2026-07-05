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
$(SIGNATURES)

Get the center (carrier) frequency of a band.

```julia-repl
julia> get_center_frequency(L1())
1575420000 Hz
```
"""
@inline get_center_frequency(::L1) = 1_575_420_000Hz

@inline get_center_frequency(::L2) = 1_227_600_000Hz

@inline get_center_frequency(::L5) = 1_176_450_000Hz

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
