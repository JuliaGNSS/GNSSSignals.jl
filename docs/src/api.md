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
