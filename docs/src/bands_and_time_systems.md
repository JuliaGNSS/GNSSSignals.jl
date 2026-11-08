# Bands and Time Systems

## Bands

```@docs
GNSSSignals.Band
GNSSSignals.L1
GNSSSignals.L2
GNSSSignals.L5
GNSSSignals.B1I
GNSSSignals.B3
GNSSSignals.E5b
GNSSSignals.E6
GNSSSignals.get_band
GNSSSignals.get_band_id
GNSSSignals.get_band_name
```

## Time Systems

Every time quantity here is plain seconds, so that it keeps full `Float64` resolution: a
[`TAITime`](@ref), such as [`get_week_start_epoch`](@ref) returns, is shifted by plain
seconds or Unitful durations and differenced into plain seconds,
[`get_time_of_week`](@ref), [`get_time_system_offset`](@ref) and
[`GNSSSignals.get_tai_offset`](@ref) return them, and
[`GNSSSignals.SECONDS_PER_WEEK`](@ref) and [`fold_week_crossover`](@ref) work in them.
The UTC dates of the scales' epochs are `DateTime` labels.

```@docs
GNSSSignals.TimeSystem
GNSSSignals.GPST
GNSSSignals.GST
GNSSSignals.BDT
GNSSSignals.get_time_system
GNSSSignals.get_time_system_id
GNSSSignals.get_time_system_name
GNSSSignals.get_time_system_offset
GNSSSignals.get_week_start_epoch
GNSSSignals.get_week
GNSSSignals.get_time_of_week
GNSSSignals.TAITime
GNSSSignals.get_tai_offset
GNSSSignals.SECONDS_PER_WEEK
GNSSSignals.fold_week_crossover
```
