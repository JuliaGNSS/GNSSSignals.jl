[![Test](https://github.com/JuliaGNSS/GNSSSignals.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaGNSS/GNSSSignals.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaGNSS/GNSSSignals.jl/graph/badge.svg?token=QY2T178W3Z)](https://codecov.io/gh/JuliaGNSS/GNSSSignals.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/GNSSSignals.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/GNSSSignals.jl/dev)

# GNSSSignals.jl

A Julia package for generating GNSS spreading codes and signals.

## Features

* GPS L1 (BPSK)
* GPS L5 (BPSK with Neuman-Hofman secondary code)
* Galileo E1B (CBOC modulation)
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

gpsl1 = GPSL1()
prn = 1

# Generate 1 ms of sampled code at 4 MHz
sampled_code = gen_code(4000, gpsl1, prn, 4MHz)

# For repeated calls, use the in-place version
buffer = zeros(Int16, 4000)
gen_code!(buffer, gpsl1, prn, 4MHz)
```

## Documentation

For detailed usage instructions and API reference, see the [documentation](https://JuliaGNSS.github.io/GNSSSignals.jl/stable).
