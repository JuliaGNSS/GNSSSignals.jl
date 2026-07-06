# Received signal power.
#
# The per-signal `get_min_received_power` values live in each signal's file
# (next to `get_code_frequency`, `get_band`, …). This file holds the generic
# instance→type forwarding and the dBW→W conversion used to write those values.

# Convert a value in dBW to a linear power in watts. Folds to a compile-time
# constant when `dbw` is a literal, so the accessors below stay as cheap as
# the other type-level constants.
@inline _dbw_to_watts(dbw) = (10.0^(dbw / 10)) * u"W"

"""
$(SIGNATURES)

Get the **minimum received power** of a GNSS signal on the ground, as a linear
power (`u"W"`).

This is the ICD-guaranteed *minimum* (worst-case satellite orientation, low
elevation), **per signal component**.

# Reference conditions

All values are stated on the **Galileo convention**: received at an ideal
`0 dBi` RHCP antenna, satellite elevation `> 5°` (Galileo OS SIS ICD v2.2,
Table 13). GPS values equal their IS-GPS ICD figures unchanged: IS-GPS also
specifies `≥ 5°`, and its reference antenna (`3 dBi`, linearly polarized) is
power-equivalent to `0 dBi` RHCP for an ideal RHCP wave — the `+3 dB` of antenna
gain exactly offsets the `−3 dB` linear-vs-circular polarization mismatch (equal
effective aperture), so no re-referencing offset is applied. Cross-constellation
comparisons are therefore on one footing; the residual from non-ideal axial
ratio has no clean single-dB correction and is treated as zero.

# Per component

Signals modelled as separate data/pilot or I/Q components carry the *component*
power. Where the ICD tabulates the component directly (GPS L1C-D/L1C-P, GPS L5
I5/Q5) that value is used; otherwise it is the total signal power times the ICD
power split.

Defined per concrete signal type; there is no generic default (a signal without
a value raises a `MethodError`). Works on an instance or a type and folds to a
compile-time constant.

```julia-repl
julia> using Unitful

julia> get_min_received_power(GPSL1CA())
1.4125375446227554e-16 W

julia> uconvert(u"dBm", get_min_received_power(GPSL1CA()))   # −128.5 dBm
-128.5 dBm
```
"""
@inline get_min_received_power(s::AbstractGNSSSignal) = get_min_received_power(typeof(s))
