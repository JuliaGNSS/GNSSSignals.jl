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
* GPS L5-I (BPSK with Neuman-Hofman secondary code)
* Galileo E1B (CBOC modulation)
* Galileo E1B BOC(1,1) approximation (lower minimum sampling rate; common SDR substitute for full CBOC)
* Highly optimized code generation: each signal bakes its fully-modulated replica into an embedded `Int8` lookup table and resamples it with a drift-free fixed-point DDA + SIMD sliding-window permute (AVX-512 / AVX2 / NEON, scalar fallback)

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

# For repeated calls, use the in-place version (output is Int8)
buffer = zeros(Int8, 4000)
gen_code!(buffer, gpsl1ca, prn, 4MHz)
```

## Documentation

For detailed usage instructions and API reference, see the [documentation](https://JuliaGNSS.github.io/GNSSSignals.jl/stable).
