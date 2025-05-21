[![Test](https://github.com/JuliaGNSS/GNSSSignals.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaGNSS/GNSSSignals.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaGNSS/GNSSSignals.jl/graph/badge.svg?token=QY2T178W3Z)](https://codecov.io/gh/JuliaGNSS/GNSSSignals.jl)

# Generate GNSS signals.

## Features

* GPS L1
* GPS L5
* Galileo E1B

## Getting started

Install:
```julia-repl
julia> ]
pkg> add GNSSSignals
```

## Usage

```julia
using GNSSSignals
code_phases = 0:1022
prn = 1
gpsl1 = GPSL1()
sampled_code = get_code.(gpsl1, code_phases, prn)
```
Output:
```julia
1023-element Array{Int8,1}:
  1
  1
  ⋮
 -1
 -1
```
In addition to that, there are some auxiliarly functions:

| Function                                                | Description                                                                        |
|---------------------------------------------------------|------------------------------------------------------------------------------------|
| `get_code_length(::AbstractGNSSSystem)`           | Get code length                                                                    |
| `get_secondary_code_length(::AbstractGNSSSystem)`  | Get secondary code length |
| `get_center_frequency(::AbstractGNSSSystem)`      | Get center frequency                                                               |
| `get_code_frequency(::AbstractGNSSSystem)`        | Get code frequency                                                                 |
| `get_data_frequency(::AbstractGNSSSystem)`        | Get data frequency                                                                 |
| `get_code(::AbstractGNSSSystem, phase, prn::Integer)` | Get code at phase `phase` from PRN `prn`                                           |
| `get_code_center_frequency_ratio(::AbstractGNSSSystem)` | Get code to center frequency ratio                                           |

#### Example

```julia-repl
julia> get_code_length(gpsl1)
1023
```
