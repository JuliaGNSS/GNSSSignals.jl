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

@inline get_center_frequency(::L5) = 1_176_450_000Hz

"""
$(SIGNATURES)

Get the center (carrier) frequency of a signal.

Dispatches through the signal's band ([`get_band`](@ref)), so all signals on
the same band return the same value by construction.

```julia-repl
julia> get_center_frequency(GPSL1CA())
1575420000 Hz
```
"""
@inline get_center_frequency(s::AbstractGNSSSignal) = get_center_frequency(get_band(s))

"""
$(SIGNATURES)

Get the [`Band`](@ref) a signal is transmitted on.

Concrete signal types define one method each, e.g. `get_band(::GPSL1CA) = L1()`.
"""
function get_band end
