[![Test](https://github.com/JuliaGNSS/GNSSSignals.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaGNSS/GNSSSignals.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaGNSS/GNSSSignals.jl/graph/badge.svg?token=QY2T178W3Z)](https://codecov.io/gh/JuliaGNSS/GNSSSignals.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/GNSSSignals.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/GNSSSignals.jl/dev)

# GNSSSignals.jl

A Julia package for generating GNSS spreading codes and signals.

## Features

* GPS L1 C/A (BPSK)
* GPS L1C-D (BOC(1,1))
* GPS L1C-P (TMBOC(6,1,4/33) with 1800-bit overlay)
* GPS L5-I (BPSK(10) with 10-bit Neuman-Hofman NH10 secondary code)
* GPS L5-Q (BPSK(10) pilot with 20-bit Neuman-Hofman NH20 secondary code)
* Galileo E1B (CBOC modulation)
* Galileo E1B BOC(1,1) approximation (lower sampling rate, Int16 output; common SDR substitute for full CBOC)
* Galileo E1C (CBOC(−) modulation, pilot component with 25-chip CS25 secondary code)
* Galileo E1C BOC(1,1) approximation (lower sampling rate, Int16 output; common SDR substitute for full CBOC)
* Galileo E5a-I (BPSK(10) with 20-bit CS20 secondary code)
* Galileo E5a-Q (BPSK(10) pilot with 100-bit per-SVID CS100 secondary code)
* Highly optimized code generation using fixed-point arithmetic for real-time signal processing

## Installation

```julia-repl
julia> ]
pkg> add GNSSSignals
```

## Quick Start

```julia
using GNSSSignals
using Unitful: MHz

gpsl1ca = GPSL1CA()
prn = 1

# Generate 1 ms of sampled code at 4 MHz
sampled_code = gen_code(4000, gpsl1ca, prn, 4MHz)

# For repeated calls, use the in-place version
buffer = zeros(Int16, 4000)
gen_code!(buffer, gpsl1ca, prn, 4MHz)
```

## Documentation

For detailed usage instructions and API reference, see the [documentation](https://JuliaGNSS.github.io/GNSSSignals.jl/stable).
