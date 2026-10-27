# Bands and Time Systems

## Bands

```@docs
GNSSSignals.Band
GNSSSignals.L1
GNSSSignals.L2
GNSSSignals.L5
GNSSSignals.B1I
GNSSSignals.B3I
GNSSSignals.E5b
GNSSSignals.E6
GNSSSignals.get_band
GNSSSignals.get_band_id
GNSSSignals.get_band_name
```

## Time Systems

Two kinds of time quantity appear here. The fixed offsets between and of the time
scales are Unitful seconds: [`get_tai_offset`](@ref), [`get_epoch_offset`](@ref) and
[`get_leap_second_offset`](@ref), and the dates of their epochs are `DateTime` labels.
Instants and times of week are plain seconds, so that they keep full `Float64`
resolution and add without unit conversions: a [`TAITime`](@ref) is shifted by and
differenced into plain seconds, [`get_time_of_week`](@ref) returns them, and
[`GNSSSignals.SECONDS_PER_WEEK`](@ref) and
[`GNSSSignals.fold_week_crossover`](@ref) work in them. `ustrip(s, offset)` turns an
offset into plain seconds where the two meet.

```@docs
GNSSSignals.TimeSystem
GNSSSignals.GPST
GNSSSignals.GST
GNSSSignals.BDT
GNSSSignals.get_time_system
GNSSSignals.get_time_system_id
GNSSSignals.get_time_system_name
GNSSSignals.get_system_start_time
GNSSSignals.get_tai_system_start_time
GNSSSignals.get_epoch_offset
GNSSSignals.get_leap_second_offset
GNSSSignals.get_week_start_epoch
GNSSSignals.get_week
GNSSSignals.get_time_of_week
GNSSSignals.TAITime
GNSSSignals.get_tai_offset
GNSSSignals.SECONDS_PER_WEEK
GNSSSignals.fold_week_crossover
```
