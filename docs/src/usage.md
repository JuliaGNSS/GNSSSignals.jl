# Usage

## Basic Code Access

The simplest way to get spreading code values is with [`get_code`](@ref):

```julia
using GNSSSignals

# Create a GNSS system instance
gpsl1 = GPSL1()

# Get the code value at phase 0 for PRN 1
code_value = get_code(gpsl1, 0.0, 1)  # Returns 1 or -1

# Get a full code period using broadcasting
code_phases = 0:1022
full_code = get_code.(gpsl1, code_phases, 1)
```

The phase is specified in chips and automatically wraps around the code length.

## Generating Sampled Code Signals

For signal processing applications, you typically need the code sampled at a specific frequency. Use [`gen_code`](@ref) or [`gen_code!`](@ref):

```julia
using GNSSSignals
using Unitful: Hz, kHz, MHz

gpsl1 = GPSL1()
prn = 1
sampling_frequency = 4MHz
num_samples = 4000  # 1 ms at 4 MHz

# Allocating version
sampled_code = gen_code(num_samples, gpsl1, prn, sampling_frequency)

# In-place version (more efficient for repeated calls)
buffer = zeros(Int16, num_samples)
gen_code!(buffer, gpsl1, prn, sampling_frequency)
```

### With Doppler Shift

To generate code with a Doppler-shifted code frequency:

```julia
using GNSSSignals
using Unitful: Hz, MHz

gpsl1 = GPSL1()
carrier_doppler = 1000Hz

# Calculate code Doppler from carrier Doppler
code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1)
code_frequency = get_code_frequency(gpsl1) + code_doppler

# Generate Doppler-shifted code
sampled_code = gen_code(4000, gpsl1, 1, 4MHz, code_frequency)
```

### With Phase Offset

You can specify a starting code phase:

```julia
start_phase = 100.5  # Start at chip 100.5
sampled_code = gen_code(4000, gpsl1, 1, 4MHz, get_code_frequency(gpsl1), start_phase)
```

## Working with Different GNSS Systems

### GPS L1

GPS L1 uses BPSK modulation with a 1023-chip C/A code:

```julia
gpsl1 = GPSL1()
get_code_length(gpsl1)           # 1023
get_center_frequency(gpsl1)      # 1575420000 Hz
get_code_frequency(gpsl1)        # 1023000 Hz
get_secondary_code_length(gpsl1) # 1 (no secondary code)
get_modulation(gpsl1)            # LOC()
```

### GPS L5

GPS L5 uses BPSK modulation with a 10230-chip code and 10-bit Neuman-Hofman secondary code:

```julia
gpsl5 = GPSL5()
get_code_length(gpsl5)           # 10230
get_center_frequency(gpsl5)      # 1176450000 Hz
get_code_frequency(gpsl5)        # 10230000 Hz
get_secondary_code_length(gpsl5) # 10
get_secondary_code(gpsl5)        # (1, 1, 1, 1, -1, -1, 1, -1, 1, -1)
```

### Galileo E1B

Galileo E1B uses CBOC(6,1,1/11) modulation:

```julia
gal_e1b = GalileoE1B()
get_code_length(gal_e1b)         # 4092
get_center_frequency(gal_e1b)    # 1575420000 Hz
get_code_frequency(gal_e1b)      # 1023000 Hz
get_modulation(gal_e1b)          # CBOC(BOCsin(1,1), BOCsin(6,1), 10/11)
```

Note that due to CBOC modulation, Galileo E1B code values are floating-point rather than integer.

## Accessing Raw Codes

To get the full code matrix directly:

```julia
gpsl1 = GPSL1()
codes = get_codes(gpsl1)  # Matrix of size (code_length, num_prns)
size(codes)               # (1023, 37)
```

Each column represents a different PRN.

## Signal Spectrum

To compute the power spectral density at a given frequency:

```julia
using GNSSSignals
using Unitful: kHz

gpsl1 = GPSL1()
psd = get_code_spectrum(gpsl1, 0kHz)  # PSD at DC
```

For custom spectrum calculations:

```julia
using GNSSSignals
using Unitful: MHz, kHz

# BPSK spectrum
psd_bpsk = get_code_spectrum_BPSK(1.023MHz, 500kHz)

# BOC spectrum
psd_boc = get_code_spectrum_BOCsin(1.023MHz, 1.023MHz, 500kHz)
```

## Performance Tips

1. **Use `gen_code!`** instead of `gen_code` when generating codes repeatedly to avoid allocations.

2. **Pre-allocate buffers** for signal generation:
   ```julia
   buffer = zeros(Int16, num_samples)
   for i in 1:num_iterations
       gen_code!(buffer, gpsl1, prn, sampling_frequency)
       # process buffer...
   end
   ```
