# Codes

## Secondary Codes

```@docs
GNSSSignals.SecondaryCode
GNSSSignals.NoSecondaryCode
GNSSSignals.SharedSecondaryCode
GNSSSignals.PerPRNSecondaryCode
```

## Modulation Types

```@docs
GNSSSignals.Modulation
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
GNSSSignals.get_code_amplitude
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
