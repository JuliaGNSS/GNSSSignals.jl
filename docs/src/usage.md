# Usage

## Generating Sampled Code Signals (Recommended)

For most applications, use [`gen_code`](@ref) or [`gen_code!`](@ref) to generate sampled codes. These functions are highly optimized for real-time GNSS signal processing and are significantly faster than calling [`get_code`](@ref) in a loop.

Each signal bakes its fully-modulated ±1 replica into an embedded `Int8` lookup table at construction; `gen_code!` resamples that table with a drift-free fixed-point DDA and a single SIMD sliding-window permute (AVX-512 / AVX2 / NEON, with a scalar fallback). The output buffer is `Int8`, and the sampling frequency must satisfy `fs ≥ code_frequency · subchip_factor`.

```julia
using GNSSSignals
using Unitful: Hz, kHz, MHz

gpsl1ca = GPSL1CA()
prn = 1
sampling_frequency = 4MHz
num_samples = 4000  # 1 ms at 4 MHz

# Allocating version
sampled_code = gen_code(num_samples, gpsl1ca, prn, sampling_frequency)

# In-place version (more efficient for repeated calls). The buffer is Int8.
buffer = zeros(Int8, num_samples)
gen_code!(buffer, gpsl1ca, prn, sampling_frequency)
```

### With Doppler Shift

To generate code with a Doppler-shifted code frequency:

```julia
using GNSSSignals
using Unitful: Hz, MHz

gpsl1ca = GPSL1CA()
prn = 1
sampling_frequency = 4MHz
carrier_doppler = 1000Hz

# Calculate code Doppler from carrier Doppler
code_doppler = carrier_doppler * get_code_center_frequency_ratio(gpsl1ca)
code_frequency = get_code_frequency(gpsl1ca) + code_doppler

# Generate Doppler-shifted code
sampled_code = gen_code(4000, gpsl1ca, prn, sampling_frequency, code_frequency)
```

### With Phase Offset

You can specify a starting code phase:

```julia
start_phase = 100.5  # Start at chip 100.5
prn = 1
sampling_frequency = 4MHz
sampled_code = gen_code(4000, gpsl1ca, prn, sampling_frequency, get_code_frequency(gpsl1ca), start_phase)
```

## Working with Different GNSS Signals

### GPS L1 C/A

GPS L1 C/A uses BPSK modulation with a 1023-chip C/A code:

```julia
gpsl1ca = GPSL1CA()
get_code_length(gpsl1ca)           # 1023
get_band(gpsl1ca)                  # L1()
get_center_frequency(gpsl1ca)      # 1575420000 Hz
get_code_frequency(gpsl1ca)        # 1023000 Hz
get_secondary_code_length(gpsl1ca) # 1 (no secondary code)
get_modulation(gpsl1ca)            # LOC()
get_signal_name(gpsl1ca)           # "GPS L1 C/A"
```

### GPS L5-I

GPS L5-I (the data-carrying component of GPS L5) uses BPSK modulation with a 10230-chip code and a 10-bit Neuman-Hofman secondary code:

```julia
gpsl5i = GPSL5I()
get_code_length(gpsl5i)           # 10230
get_band(gpsl5i)                  # L5()
get_center_frequency(gpsl5i)      # 1176450000 Hz
get_code_frequency(gpsl5i)        # 10230000 Hz
get_secondary_code_length(gpsl5i) # 10
get_secondary_code(gpsl5i)        # (1, 1, 1, 1, -1, -1, 1, -1, 1, -1)
```

### GPS L1C-D

GPS L1C-D is the data-carrying component of GPS L1C (IS-GPS-800G). It uses BOC(1,1) modulation with a 10230-chip Weil-based primary code and broadcasts the CNAV-2 message at 100 sps (post-LDPC channel symbol rate; the information bit rate after rate-½ decoding is 50 bps). It is broadcast by Block III/IIIF satellites and supports PRNs 1-63:

```julia
gpsl1c_d = GPSL1C_D()
get_code_length(gpsl1c_d)           # 10230
get_band(gpsl1c_d)                  # L1()
get_center_frequency(gpsl1c_d)      # 1575420000 Hz
get_code_frequency(gpsl1c_d)        # 1023000 Hz
get_data_frequency(gpsl1c_d)        # 100 Hz
get_secondary_code_length(gpsl1c_d) # 1 (no secondary code)
get_modulation(gpsl1c_d)            # BOCsin(1, 1)
```

### GPS L1C-P

GPS L1C-P is the pilot (dataless) component of GPS L1C. It carries 75% of the L1C power per IS-GPS-800G and uses TMBOC(6,1,4/33) modulation: every 33 primary chips, four positions (`{0, 4, 6, 29}`) use BOC(6,1) and the rest use BOC(1,1). It shares the same 10230-chip Weil construction as L1C-D (different per-PRN parameters) and adds an 18-second 1800-bit per-PRN overlay (exposed as a [`PerPRNSecondaryCode`](@ref GNSSSignals.PerPRNSecondaryCode)):

```julia
gpsl1c_p = GPSL1C_P()
get_code_length(gpsl1c_p)           # 10230
get_band(gpsl1c_p)                  # L1()
get_center_frequency(gpsl1c_p)      # 1575420000 Hz
get_code_frequency(gpsl1c_p)        # 1023000 Hz
get_data_frequency(gpsl1c_p)        # 0 Hz (dataless)
get_secondary_code_length(gpsl1c_p) # 1800
get_modulation(gpsl1c_p)            # TMBOC(BOCsin(1,1), BOCsin(6,1), …)
```

The overlay code is per-PRN, so [`get_secondary_code`](@ref) returns a [`PerPRNSecondaryCode`](@ref GNSSSignals.PerPRNSecondaryCode) wrapping the 1800 × 63 overlay matrix rather than a plain tuple.

### Galileo E1B

Galileo E1B (the data-carrying component of Galileo E1 OS) uses CBOC(6,1,1/11) modulation. It is transmitted on the same RF carrier as GPS L1 C/A — [`get_band`](@ref) returns [`L1`](@ref GNSSSignals.L1) for both:

```julia
gal_e1b = GalileoE1B()
get_code_length(gal_e1b)         # 4092
get_band(gal_e1b)                # L1()
get_center_frequency(gal_e1b)    # 1575420000 Hz
get_code_frequency(gal_e1b)      # 1023000 Hz
get_modulation(gal_e1b)          # CBOC(BOCsin(1,1), BOCsin(6,1), 10/11)
```

Note that the single-chip accessor [`get_code`](@ref) returns floating-point values for Galileo E1B (the CBOC subcarrier amplitudes are irrational), and `get_code_type(GalileoE1B())` is `Float32`. [`gen_code!`](@ref), however, always outputs `Int8`: for CBOC it bakes an `Int8` integer approximation of the two sub-carrier amplitudes (default `(19, 6) → ±25, ±13`), which is sign-identical to the spec with ~0 dB correlation loss.

### Galileo E1B (BOC(1,1) approximation)

Many software receivers substitute a BOC(1,1) replica for the full CBOC(6,1,1/11) E1B signal: BOC(6,1) carries only 1/11 of the signal power, so the correlation loss is roughly 0.45 dB and the replica needs only `fs ≥ 2 · 1.023` MHz instead of `fs ≥ 12 · 1.023` MHz. PocketSDR and other open-source GNSS-SDRs default to this approximation. [`GalileoE1B_BOC11`](@ref GNSSSignals.GalileoE1B_BOC11) is the typed variant: identical primary code, code length, code rate, data rate, and band as [`GalileoE1B`](@ref), but with `BOCsin(1, 1)` as the modulation:

```julia
e1b = GalileoE1B_BOC11()
get_code_length(e1b)         # 4092
get_modulation(e1b)          # BOCsin(1, 1)
get_code_type(e1b)           # Int16  (the get_code accessor type; gen_code! is Int8)
```

[`gen_code!`](@ref) outputs `Int8` for both variants. The reason to choose [`GalileoE1B_BOC11`](@ref GNSSSignals.GalileoE1B_BOC11) over [`GalileoE1B`](@ref) is the **much lower minimum sampling rate** (`fs ≥ 2 · 1.023` MHz vs `12 · 1.023` MHz) when the ~0.45 dB loss is acceptable.

## Bands

A [`Band`](@ref GNSSSignals.Band) represents a shared RF carrier frequency. Two signals with the same band can be driven by a single carrier NCO in a receiver — that is the architectural reason this abstraction exists.

```julia
get_band(GPSL1CA())            # L1()
get_band(GalileoE1B())         # L1()
get_band(GPSL5I())             # L5()

get_center_frequency(L1())     # 1575420000 Hz
get_center_frequency(L5())     # 1176450000 Hz
```

Band identity here is by RF frequency, not by ICD label: Galileo E1 returns `L1()` because it shares 1575.42 MHz with GPS L1.

## Basic Code Access

For accessing individual code values at specific phases (e.g., for analysis or custom resampling), use [`get_code`](@ref):

```julia
using GNSSSignals

gpsl1ca = GPSL1CA()
prn = 1

# Get a single code value at phase 0.0 for PRN = prn
code_value = get_code(gpsl1ca, 0.0, prn)  # Returns 1 or -1

# Get a full code period using broadcasting
code_phases = 0:1022
full_code = get_code.(gpsl1ca, code_phases, prn)
```

The phase is specified in chips and automatically wraps around the code length.

!!! note
    Prefer [`gen_code`](@ref) or [`gen_code!`](@ref) for generating sampled codes at a specific sampling frequency, as they are significantly faster.

## Accessing Raw Codes

To get the full code matrix directly:

```julia
gpsl1ca = GPSL1CA()
codes = get_codes(gpsl1ca)  # Matrix of size (code_length, num_prns)
size(codes)                 # (1023, 37)
```

Each column represents a different PRN.

## Signal Spectrum

To compute the power spectral density at a given frequency:

```julia
using GNSSSignals
using Unitful: kHz

gpsl1ca = GPSL1CA()
psd = get_code_spectrum(gpsl1ca, 0kHz)  # PSD at DC
```

For custom spectrum calculations:

```julia
using GNSSSignals
using Unitful: MHz, kHz

# BPSK spectrum
psd_bpsk = GNSSSignals.get_code_spectrum_BPSK(1.023MHz, 500kHz)

# BOC spectrum
psd_boc = GNSSSignals.get_code_spectrum_BOCsin(1.023MHz, 1.023MHz, 500kHz)
```

## Performance Tips

1. **Use `gen_code!`** instead of `gen_code` when generating codes repeatedly to avoid allocations.

2. **Pre-allocate buffers** for signal generation:
   ```julia
   num_iterations = 1000
   buffer = zeros(Int8, num_samples)
   for i in 1:num_iterations
       gen_code!(buffer, gpsl1ca, prn, sampling_frequency)
       # process buffer...
   end
   ```

   For repeated integrations that continue seamlessly across blocks (tracking), build a
   [`code_engine`](@ref)`(signal, prn, fs, fc)` once, seed `st = code_state(eng)`, and call
   `st = gen_code!(buffer, eng, st)` per block — the DDA setup is paid once, the per-call fill
   is allocation-free, and threading the returned state continues the stream seamlessly.
