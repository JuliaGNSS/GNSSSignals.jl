# API Reference

## GNSS Signal Types

```@docs
GNSSSignals.AbstractGNSSSignal
GNSSSignals.AbstractGPSSignal
GNSSSignals.AbstractGalileoSignal
GNSSSignals.GPSL1CA
GNSSSignals.GPSL1C_D
GNSSSignals.GPSL1C_P
GNSSSignals.GPSL2CM
GNSSSignals.GPSL2CL
GNSSSignals.GPSL5I
GNSSSignals.GPSL5Q
GNSSSignals.GalileoE1B
GNSSSignals.GalileoE1B_BOC11
GNSSSignals.GalileoE1C
GNSSSignals.GalileoE1C_BOC11
GNSSSignals.GalileoE5aI
GNSSSignals.GalileoE5aQ
```

## Bands

```@docs
GNSSSignals.Band
GNSSSignals.L1
GNSSSignals.L2
GNSSSignals.L5
GNSSSignals.get_band
GNSSSignals.get_band_id
```

## Time Systems

```@docs
GNSSSignals.TimeSystem
GNSSSignals.GPST
GNSSSignals.GST
GNSSSignals.get_time_system
GNSSSignals.get_system_start_time
GNSSSignals.get_tai_offset
```

## Secondary Codes

```@docs
GNSSSignals.SecondaryCode
GNSSSignals.NoSecondaryCode
GNSSSignals.SharedSecondaryCode
GNSSSignals.PerPRNSecondaryCode
```

## Modulation Types

```@docs
GNSSSignals.LOC
GNSSSignals.BOC
GNSSSignals.BOCsin
GNSSSignals.BOCcos
GNSSSignals.CBOC
GNSSSignals.TMBOC
```

## Code Generation

```@docs
GNSSSignals.gen_code
GNSSSignals.gen_code!
```

## LUT Code Generation

`gen_code!`/`gen_code` (above) resample PRN `prn`'s fully-modulated replica from the signal's
embedded `Int8` LUT. For continuing, block-to-block generation (tracking) build a
[`code_engine`](@ref)`(signal, prn, fs, fc)` once and thread the immutable state: seed it with
[`code_state`](@ref)`(eng)`, then call `gen_code!(out, eng, st)` per integration, threading the
returned `CodeFillState`. For an allocation-free, register-resident fused loop, use the
value-based [`code_engine`](@ref)`(signal, prn, fs, fc, Val(K))` and its state stepping
([`code_state`](@ref) / [`code_lookup`](@ref) / [`code_advance`](@ref) / [`code_width`](@ref)).
Both continuing paths use immutable, explicit state — nothing mutates behind your back.

```@docs
GNSSSignals.CodeFillEngine
GNSSSignals.CodeFillState
GNSSSignals.code_engine
GNSSSignals.code_state
GNSSSignals.code_lookup
GNSSSignals.code_advance
GNSSSignals.code_width
```

## Code Access

```@docs
GNSSSignals.get_code
GNSSSignals.get_code_unsafe
GNSSSignals.get_codes
```

## Signal Properties

```@docs
GNSSSignals.get_code_length
GNSSSignals.get_secondary_code_length
GNSSSignals.get_secondary_code
GNSSSignals.get_center_frequency
GNSSSignals.get_code_frequency
GNSSSignals.get_data_frequency
GNSSSignals.get_code_center_frequency_ratio
GNSSSignals.get_modulation
GNSSSignals.get_signal_name
GNSSSignals.get_signal_id
GNSSSignals.min_bits_for_code_length
GNSSSignals.get_code_type
```

## Spectrum Functions

```@docs
GNSSSignals.get_code_spectrum
GNSSSignals.get_code_spectrum_BPSK
GNSSSignals.get_code_spectrum_BOCsin
GNSSSignals.get_code_spectrum_BOCcos
```

## Utility Functions

```@docs
GNSSSignals.read_in_codes
```
