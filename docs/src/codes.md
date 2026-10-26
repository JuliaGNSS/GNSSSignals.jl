# Codes

## Secondary Codes

```@docs
GNSSSignals.SecondaryCode
GNSSSignals.NoSecondaryCode
GNSSSignals.SharedSecondaryCode
GNSSSignals.PerPRNSecondaryCode
GNSSSignals.secondary_value
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
[`GNSSSignals.code_engine`](@ref)`(signal, prn, fs, fc)` once and thread the immutable state: seed it with
[`GNSSSignals.code_state`](@ref)`(eng)`, then call `gen_code!(out, eng, st)` per integration, threading the
returned `CodeFillState`. For an allocation-free, register-resident fused loop, use the
value-based [`GNSSSignals.code_engine`](@ref)`(signal, prn, fs, fc, Val(K))` and its state stepping
([`GNSSSignals.code_state`](@ref) / [`GNSSSignals.code_lookup`](@ref) / [`GNSSSignals.code_advance`](@ref) / [`GNSSSignals.code_width`](@ref)).
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
