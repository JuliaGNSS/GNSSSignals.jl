# API Reference

## GNSS Signal Types

```@docs
GNSSSignals.AbstractGNSSSignal
GNSSSignals.GPSL1CA
GNSSSignals.GPSL1C_D
GNSSSignals.GPSL1C_P
GNSSSignals.GPSL5I
GNSSSignals.GalileoE1B
GNSSSignals.GalileoE1B_BOC11
```

## Bands

```@docs
GNSSSignals.Band
GNSSSignals.L1
GNSSSignals.L5
GNSSSignals.get_band
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

Fast register-resident SIMD code resampling: build a [`CodeReplicaLUT`](@ref) plan once per
`(signal, prn)`, then resample it with `gen_code!`, the continuing [`CodeGeneratorLUT`](@ref),
or — for an allocation-free, fused loop — the value-based [`code_engine`](@ref) and its
state stepping ([`code_state`](@ref) / [`code_lookup`](@ref) / [`code_advance`](@ref) /
[`code_width`](@ref)).

```@docs
GNSSSignals.CodeReplicaLUT
GNSSSignals.CodeGeneratorLUT
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
