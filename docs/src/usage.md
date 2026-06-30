# Usage

## Generating Sampled Code Signals (Recommended)

For most applications, use [`gen_code`](@ref) or [`gen_code!`](@ref) to generate sampled codes. These functions are highly optimized for real-time GNSS signal processing, using fixed-point arithmetic and minimizing memory access. They are significantly faster than calling [`get_code`](@ref) in a loop.

These functions exploit the fact that consecutive samples often map to the same code chip, avoiding per-sample floating-point operations and modulo calculations.

```julia
using GNSSSignals
using Unitful: Hz, kHz, MHz

gpsl1ca = GPSL1CA()
prn = 1
sampling_frequency = 4MHz
num_samples = 4000  # 1 ms at 4 MHz

# Allocating version
sampled_code = gen_code(num_samples, gpsl1ca, prn, sampling_frequency)

# In-place version (more efficient for repeated calls)
buffer = zeros(Int16, num_samples)
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

### GPS L2 CM

GPS L2 CM is the moderate-length, data-carrying component of the GPS L2 civil signal (L2C), on the L2 band (1227.6 MHz). It is a BPSK 10230-chip code at 511.5 kcps (a 20 ms period) carrying the CNAV message at 50 sps, with no secondary code (IS-GPS-200N §3.2.1.4). In the broadcast signal it is time-multiplexed chip-by-chip with the [`GPSL2CL`](@ref GNSSSignals.GPSL2CL) pilot; like GNSS-SDR and PocketSDR this implementation models the CM component on its own at its native chip rate. PRNs 1-63 are supported:

```julia
gpsl2cm = GPSL2CM()
get_code_length(gpsl2cm)            # 10230
get_band(gpsl2cm)                   # L2()
get_center_frequency(gpsl2cm)       # 1227600000 Hz
get_code_frequency(gpsl2cm)         # 511500 Hz
get_data_frequency(gpsl2cm)         # 50 Hz (CNAV symbol rate)
get_secondary_code_length(gpsl2cm)  # 1 (no secondary code)
get_modulation(gpsl2cm)             # LOC()
```

### GPS L2 CL

GPS L2 CL is the long, dataless pilot component of the GPS L2 civil signal. It shares the L2 CM code generator (different per-PRN initial state) but is short-cycled at a much longer 767250-chip period — a 1.5 s code at 511.5 kcps — and carries no data (IS-GPS-200N §3.2.1.5):

```julia
gpsl2cl = GPSL2CL()
get_code_length(gpsl2cl)            # 767250
get_band(gpsl2cl)                   # L2()
get_center_frequency(gpsl2cl)       # 1227600000 Hz
get_code_frequency(gpsl2cl)         # 511500 Hz
get_data_frequency(gpsl2cl)         # 0 Hz (dataless)
get_secondary_code_length(gpsl2cl)  # 1 (no secondary code)
get_modulation(gpsl2cl)             # LOC()
```

The full L2 CL code matrix is 767250 × 63; stored as `Int16` it occupies ~97 MB, so build `GPSL2CL()` once and reuse the instance.

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

Note that due to CBOC modulation, Galileo E1B code values are floating-point rather than integer.

### Galileo E1B (BOC(1,1) approximation)

Many software receivers substitute a BOC(1,1) replica for the full CBOC(6,1,1/11) E1B signal: BOC(6,1) carries only 1/11 of the signal power, so the correlation loss is roughly 0.45 dB and the replica needs only `fs ≥ 2 · 1.023` MHz instead of `fs ≥ 12 · 1.023` MHz. PocketSDR and other open-source GNSS-SDRs default to this approximation. [`GalileoE1B_BOC11`](@ref GNSSSignals.GalileoE1B_BOC11) is the typed variant: identical primary code, code length, code rate, data rate, and band as [`GalileoE1B`](@ref), but with `BOCsin(1, 1)` as the modulation and integer (`Int16`) output:

```julia
e1b = GalileoE1B_BOC11()
get_code_length(e1b)         # 4092
get_modulation(e1b)          # BOCsin(1, 1)
get_code_type(e1b)           # Int16
```

Use [`GalileoE1B`](@ref) for the full CBOC spec output (Float32); use [`GalileoE1B_BOC11`](@ref GNSSSignals.GalileoE1B_BOC11) when the 0.45 dB loss is acceptable and Int16 / lower sampling rate is preferable.

### Galileo E1C

Galileo E1C is the pilot (dataless) component of Galileo E1 OS. It shares the 4092-chip primary memory codes' construction with E1B but uses *different* codes (Galileo OS SIS ICD Annex C), and it carries a 25-chip secondary code (CS25) overlaid one chip per 4 ms primary period, giving a 100 ms tiered code. Its CBOC uses the BOC(6,1) component in anti-phase — CBOC(−) — where E1B uses CBOC(+) (Galileo OS SIS ICD §2.3.3):

```julia
gal_e1c = GalileoE1C()
get_code_length(gal_e1c)           # 4092
get_secondary_code_length(gal_e1c) # 25 (CS25)
get_band(gal_e1c)                  # L1()
get_data_frequency(gal_e1c)        # 0 Hz (pilot)
get_modulation(gal_e1c)            # CBOC(BOCsin(1,1), BOCsin(6,1), 10/11, -1)
```

The secondary code is shared across all PRNs, so [`get_secondary_code`](@ref) returns a [`SharedSecondaryCode`](@ref GNSSSignals.SharedSecondaryCode) of length 25. As with E1B, CBOC modulation makes the code values floating-point.

A [`GalileoE1C_BOC11`](@ref GNSSSignals.GalileoE1C_BOC11) variant provides the BOC(1,1) approximation (Int16 output, lower sampling rate) — the same substitution PocketSDR uses for E1C by default — with identical primary and CS25 secondary codes:

```julia
e1c = GalileoE1C_BOC11()
get_modulation(e1c)          # BOCsin(1, 1)
get_code_type(e1c)           # Int16
```

### Galileo E5a-I

Galileo E5a-I (the data-carrying component of Galileo E5a) uses a 10230-chip primary code at 10.23 Mcps with a 20-bit CS20 secondary code (shared by all SVIDs, giving a 20 ms tiered code). It is transmitted on the same RF carrier as GPS L5 — [`get_band`](@ref) returns [`L5`](@ref GNSSSignals.L5) for both. The wideband E5 signal is AltBOC(15,10), but like GNSS-SDR and PocketSDR this implementation models the E5a sideband on its own as BPSK(10):

```julia
e5a_i = GalileoE5aI()
get_code_length(e5a_i)           # 10230
get_band(e5a_i)                  # L5()
get_center_frequency(e5a_i)      # 1176450000 Hz
get_code_frequency(e5a_i)        # 10230000 Hz
get_data_frequency(e5a_i)        # 50 Hz (F/NAV symbol rate)
get_secondary_code_length(e5a_i) # 20
get_modulation(e5a_i)            # LOC()
```

### Galileo E5a-Q

Galileo E5a-Q is the pilot (dataless) component of Galileo E5a. It shares the E5a-I primary-code generator (different per-SVID register seeds) and overlays a 100-bit per-SVID CS100 secondary code (a 100 ms tiered code, exposed as a [`PerPRNSecondaryCode`](@ref GNSSSignals.PerPRNSecondaryCode)). PRNs 1-50 are supported:

```julia
e5a_q = GalileoE5aQ()
get_code_length(e5a_q)           # 10230
get_band(e5a_q)                  # L5()
get_center_frequency(e5a_q)      # 1176450000 Hz
get_code_frequency(e5a_q)        # 10230000 Hz
get_data_frequency(e5a_q)        # 0 Hz (dataless)
get_secondary_code_length(e5a_q) # 100
get_modulation(e5a_q)            # LOC()
```

## Bands

A [`Band`](@ref GNSSSignals.Band) represents a shared RF carrier frequency. Two signals with the same band can be driven by a single carrier NCO in a receiver — that is the architectural reason this abstraction exists.

```julia
get_band(GPSL1CA())            # L1()
get_band(GalileoE1B())         # L1()
get_band(GalileoE1C())         # L1()
get_band(GPSL2CM())            # L2()
get_band(GPSL5I())             # L5()

get_center_frequency(L1())     # 1575420000 Hz
get_center_frequency(L2())     # 1227600000 Hz
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
   buffer = zeros(Int16, num_samples)
   for i in 1:num_iterations
       gen_code!(buffer, gpsl1ca, prn, sampling_frequency)
       # process buffer...
   end
   ```
