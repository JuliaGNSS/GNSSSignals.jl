# API Reference

## GNSS System Types

```@docs
GNSSSignals.AbstractGNSS
GNSSSignals.GPSL1
GNSSSignals.GPSL5
GNSSSignals.GalileoE1B
```

## Modulation Types

```@docs
GNSSSignals.LOC
GNSSSignals.BOCsin
GNSSSignals.BOCcos
GNSSSignals.CBOC
```

## Code Generation

```@docs
GNSSSignals.gen_code
GNSSSignals.gen_code!
```

## Code Access

```@docs
GNSSSignals.get_code
GNSSSignals.get_code_unsafe
GNSSSignals.get_codes
```

## System Properties

```@docs
GNSSSignals.get_code_length
GNSSSignals.get_secondary_code_length
GNSSSignals.get_secondary_code
GNSSSignals.get_center_frequency
GNSSSignals.get_code_frequency
GNSSSignals.get_data_frequency
GNSSSignals.get_code_center_frequency_ratio
GNSSSignals.get_modulation
GNSSSignals.get_system_string
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
