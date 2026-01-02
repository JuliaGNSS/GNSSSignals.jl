# GNSSSignals.jl

A Julia package for generating GNSS (Global Navigation Satellite System) spreading codes and signals.

## Features

### Supported GNSS Systems

| System | Carrier Frequency | Code Frequency | Code Length | Modulation |
|--------|-------------------|----------------|-------------|------------|
| GPS L1 | 1575.42 MHz | 1.023 MHz | 1023 chips | BPSK |
| GPS L5 | 1176.45 MHz | 10.23 MHz | 10230 chips | BPSK + Neuman-Hofman |
| Galileo E1B | 1575.42 MHz | 1.023 MHz | 4092 chips | CBOC(6,1,1/11) |

### Modulation Types

- **LOC** - Legacy/BPSK modulation (GPS L1, GPS L5)
- **BOCsin** - Sine-phased Binary Offset Carrier
- **BOCcos** - Cosine-phased Binary Offset Carrier
- **CBOC** - Composite Binary Offset Carrier (Galileo E1B)

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

gpsl1 = GPSL1()
prn = 1

# Generate 1 ms of sampled code at 4 MHz (recommended approach)
sampled_code = gen_code(4000, gpsl1, prn, 4MHz)

# For repeated calls, use the in-place version with a pre-allocated buffer
buffer = zeros(Int16, 4000)
gen_code!(buffer, gpsl1, prn, 4MHz)
```

For accessing individual code values at specific phases, use [`get_code`](@ref):

```julia
# Get a single code value
code_value = get_code(gpsl1, 0.0, prn)  # Returns 1 or -1

# Get a full code period (slower than gen_code for large arrays)
code_phases = 0:1022
full_code = get_code.(gpsl1, code_phases, prn)
```

For more detailed examples, see the [Usage](@ref) guide.

## Package Overview

GNSSSignals.jl provides functionality to:

- Generate spreading codes for GPS L1, GPS L5, and Galileo E1B
- Sample codes at arbitrary frequencies with [`gen_code`](@ref) and [`gen_code!`](@ref)
- Access code values at specific phases with [`get_code`](@ref)
- Query system parameters (code length, frequencies, modulation type)
- Compute signal spectra for analysis

The [`gen_code`](@ref) and [`gen_code!`](@ref) functions are highly optimized for real-time GNSS signal processing, using fixed-point arithmetic and minimizing memory access by exploiting the fact that consecutive samples often map to the same code chip.

The package is designed to be used as a building block for GNSS receiver implementations and signal processing research.
