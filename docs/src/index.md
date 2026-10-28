# GNSSSignals.jl

A Julia package for generating GNSS (Global Navigation Satellite System) spreading codes and signals.

## Features

### Supported GNSS Signals

| Signal | Band | Carrier Frequency | Chip Rate | Code Length | Secondary Code | Modulation |
|--------|------|-------------------|----------------|-------------|----------------|------------|
| GPS L1 C/A | L1 | 1575.42 MHz | 1.023 MHz | 1023 chips | — | BPSK |
| GPS L1C-D | L1 | 1575.42 MHz | 1.023 MHz | 10230 chips | — | BOC(1,1) |
| GPS L1C-P | L1 | 1575.42 MHz | 1.023 MHz | 10230 chips | 1800-chip overlay | TMBOC(6,1,4/33) |
| GPS L2 CM | L2 | 1227.6 MHz | 511.5 kHz | 10230 chips | — | BPSK |
| GPS L2 CL | L2 | 1227.6 MHz | 511.5 kHz | 767250 chips | — | BPSK (pilot) |
| GPS L5-I | L5 | 1176.45 MHz | 10.23 MHz | 10230 chips | 10-chip NH10 | BPSK |
| GPS L5-Q | L5 | 1176.45 MHz | 10.23 MHz | 10230 chips | 20-chip NH20 | BPSK (pilot) |
| Galileo E1B | L1 | 1575.42 MHz | 1.023 MHz | 4092 chips | — | CBOC(6,1,1/11) |
| Galileo E1B (BOC(1,1) approximation) | L1 | 1575.42 MHz | 1.023 MHz | 4092 chips | — | BOC(1,1) |
| Galileo E1C | L1 | 1575.42 MHz | 1.023 MHz | 4092 chips | 25-chip CS25 | CBOC(6,1,1/11), BOC(6,1) anti-phase (pilot) |
| Galileo E1C (BOC(1,1) approximation) | L1 | 1575.42 MHz | 1.023 MHz | 4092 chips | 25-chip CS25 | BOC(1,1) (pilot) |
| Galileo E5a-I | L5 | 1176.45 MHz | 10.23 MHz | 10230 chips | 20-chip CS20 | BPSK |
| Galileo E5a-Q | L5 | 1176.45 MHz | 10.23 MHz | 10230 chips | 100-chip CS100 | BPSK (pilot) |
| Galileo E5a-QP | L5 | 1176.45 MHz | 5.115 MHz | 330 chips | — | BPSK (quasi-pilot) |
| Galileo E5b-I | E5b | 1207.14 MHz | 10.23 MHz | 10230 chips | 4-chip CS4 | BPSK |
| Galileo E5b-Q | E5b | 1207.14 MHz | 10.23 MHz | 10230 chips | 100-chip CS100 | BPSK (pilot) |
| Galileo E6-B | E6 | 1278.75 MHz | 5.115 MHz | 5115 chips | — | BPSK |
| Galileo E6-C | E6 | 1278.75 MHz | 5.115 MHz | 5115 chips | 100-chip CS100 | BPSK (pilot) |
| BeiDou B1I | B1I | 1561.098 MHz | 2.046 MHz | 2046 chips | 20-chip NH20 (D1) | BPSK |
| BeiDou B3I | B3I | 1268.52 MHz | 10.23 MHz | 10230 chips | 20-chip NH20 (D1) | BPSK |
| BeiDou B2b-I | E5b | 1207.14 MHz | 10.23 MHz | 10230 chips | — | BPSK |
| BeiDou B2a data | L5 | 1176.45 MHz | 10.23 MHz | 10230 chips | 5-chip | BPSK |
| BeiDou B2a pilot | L5 | 1176.45 MHz | 10.23 MHz | 10230 chips | 100-chip | BPSK (pilot) |
| BeiDou B1C data | L1 | 1575.42 MHz | 1.023 MHz | 10230 chips | — | BOC(1,1) |
| BeiDou B1C pilot | L1 | 1575.42 MHz | 1.023 MHz | 10230 chips | 1800-chip | BOC(1,1) approximation of QMBOC(6,1,4/33) (pilot) |

Signals on the same RF carrier report the same [`Band`](@ref GNSSSignals.Band): GPS L1 C/A, GPS L1C, Galileo E1 and BeiDou B1C are all [`L1`](@ref GNSSSignals.L1) at 1575.42 MHz, so a receiver can drive them with the same carrier NCO. Likewise every signal reports the [`TimeSystem`](@ref GNSSSignals.TimeSystem) of its constellation — [`GPST`](@ref GNSSSignals.GPST), [`GST`](@ref GNSSSignals.GST) or [`BDT`](@ref GNSSSignals.BDT) — the level at which signals share a receiver clock bias.

### Modulation Types

- **LOC** — Linear Offset Carrier; the BPSK-like baseline with no subcarrier (GPS L1 C/A, L2C and L5, Galileo E5 and E6, BeiDou B1I, B2 and B3I)
- **BOCsin** — Sine-phased Binary Offset Carrier (GPS L1C-D, BeiDou B1C, the Galileo E1 BOC(1,1) approximations)
- **BOCcos** — Cosine-phased Binary Offset Carrier
- **CBOC** — Composite Binary Offset Carrier (Galileo E1B uses CBOC(+), Galileo E1C uses CBOC(−))
- **TMBOC** — Time-Multiplexed Binary Offset Carrier (GPS L1C-P)

## Installation

```julia
using Pkg
Pkg.add("GNSSSignals")
```

Or from the Julia REPL:

```julia-repl
julia> ]
pkg> add GNSSSignals
```

## Quick Start

For most applications, use [`gen_code`](@ref) or [`gen_code!`](@ref) to generate sampled codes. These functions are highly optimized and significantly faster than calling [`get_code`](@ref) in a loop.

```julia
using GNSSSignals
using Unitful: MHz

gpsl1ca = GPSL1CA()
prn = 1

# Generate 1 ms of sampled code at 4 MHz (recommended approach)
sampled_code = gen_code(4000, gpsl1ca, prn, 4MHz)

# For repeated calls, use the in-place version with a pre-allocated buffer.
# Output is Int8 (±1 for BPSK/BOC, an integer approximation for CBOC).
buffer = zeros(Int8, 4000)
gen_code!(buffer, gpsl1ca, prn, 4MHz)
```

For accessing individual code values at specific phases, use [`get_code`](@ref):

```julia
# Get a single code value
code_value = get_code(gpsl1ca, 0.0, prn)  # Returns 1 or -1

# Get a full code period (slower than gen_code for large arrays)
code_phases = 0:1022
full_code = get_code.(gpsl1ca, code_phases, prn)
```

For more detailed examples, see the [Usage](@ref) guide.

## Package Overview

GNSSSignals.jl provides functionality to:

- Generate spreading codes for every signal in the table above
- Sample codes at arbitrary frequencies with [`gen_code`](@ref) and [`gen_code!`](@ref)
- Access code values at specific phases with [`get_code`](@ref)
- Query signal parameters (code length, frequencies, modulation type, relative power, carrier phase offset)
- Identify signals, bands and constellations by `Symbol` id or display name
- Group signals by RF band via [`get_band`](@ref) and by time system via [`get_time_system`](@ref)
- Relate the GPS, Galileo and BeiDou time scales to each other and to TAI
- Compute signal spectra for analysis

The [`gen_code`](@ref) and [`gen_code!`](@ref) functions are highly optimized for real-time GNSS signal processing: each signal bakes its fully-modulated ±1 replica into an embedded `Int8` lookup table once, then resamples it with a drift-free fixed-point DDA and a single SIMD sliding-window permute (AVX-512 / AVX2 / NEON, scalar fallback). Output is `Int8`.

The package is designed to be used as a building block for GNSS receiver implementations and signal processing research.
