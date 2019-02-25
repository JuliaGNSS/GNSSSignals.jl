[![Build Status](https://travis-ci.org/JuliaGNSS/GNSSSignals.jl.svg?branch=master)](https://travis-ci.org/JuliaGNSS/GNSSSignals.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNSS/GNSSSignals.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGNSS/GNSSSignals.jl?branch=master)

# Generate GNSS signals.

## Features

* GPS L1
* GPS L5

## Getting started

Install:
```julia-repl
pkg> add https://github.com/JuliaGNSS/GNSSSignals.jl.git
```

## Usage

```julia
using GNSSSignals
gpsl1 = GPSL1()
code_freq = 1023e3
code_phase = 4
sample_freq = 4e6
prn = 1
sampled_code = gen_code.(Ref(gpsl1), 1:4000, code_freq, code_phase, sample_freq, prn)
```
Output:
```julia
4000-element Array{Int8,1}:
  1
  1
  1
  ⋮
 -1
 -1
  1
```

## Todo

* Galileo signals

## License

MIT License
