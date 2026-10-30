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
[`get_time_scale_offset`](@ref), and the dates of their epochs are `DateTime` labels.
Instants use [`TAITime`](@ref), which accepts either plain seconds or Unitful durations
for shifts and returns differences in plain seconds. [`get_time_of_week`](@ref) returns
plain seconds too, and
[`GNSSSignals.SECONDS_PER_WEEK`](@ref) and
[`GNSSSignals.fold_week_crossover`](@ref) work in them. `ustrip(s, offset)` turns an
offset into plain seconds where needed. [`get_week_and_time_of_week`](@ref) returns a
normalized pair: a TOW rounded to the week endpoint carries into the next week.

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
GNSSSignals.get_time_scale_offset
GNSSSignals.get_week_start_epoch
GNSSSignals.get_week_and_time_of_week
GNSSSignals.get_week
GNSSSignals.get_time_of_week
GNSSSignals.TAITime
GNSSSignals.get_tai_offset
GNSSSignals.SECONDS_PER_WEEK
GNSSSignals.fold_week_crossover
```
