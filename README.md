[![Build Status](https://travis-ci.org/JuliaGNSS/GNSSSignals.jl.svg?branch=master)](https://travis-ci.org/JuliaGNSS/GNSSSignals.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNSS/GNSSSignals.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGNSS/GNSSSignals.jl?branch=master)

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
sampled_code = get_code.(GPSL1, code_phases, prn)
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
| `get_code_length(::Type{<:AbstractGNSSSystem})`           | Get code length                                                                    |
| `get_secondary_code_length(::Type{<:AbstractGNSSSystem})`  | Get secondary code length |
| `get_center_frequency(::Type{<:AbstractGNSSSystem})`      | Get center frequency                                                               |
| `get_code_frequency(::Type{<:AbstractGNSSSystem})`        | Get code frequency                                                                 |
| `get_data_frequency(::Type{<:AbstractGNSSSystem})`        | Get data frequency                                                                 |
| `get_code(::Type{<:AbstractGNSSSystem}, phase, prn::Int)` | Get code at phase `phase` from PRN `prn`                                           |
| `get_code_center_frequency_ratio(::Type{<:AbstractGNSSSystem})` | Get code to center frequency ratio                                           |

#### Example

```julia-repl
julia> get_code_length(GPSL1)
1023
```
