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
gen_sampled_code, get_code_phase = init_gpsl1_codes()
sat_prn = 1
code = gen_sampled_code(1:4000, 1023e3, 20, 4e6, sat_prn)
code_phase = get_code_phase(4000, 1023e3, 20, 4e6)
carrier = gen_carrier(1:4000, 1e3, 20 * pi / 180, 4e6)
carrier_phase = get_carrier_phase(4000, 1e3, 20 * pi / 180, 4e6)

sat_prn = 4
gen_gps_code, get_code_phase = init_gpsl5_codes()
gen_gps_code(1:10230, 1023e4, 20, 4e6, sat_prn)
get_code_phase(10230, 1023e4, 20, 4e6)
```

## Todo

* Galileo signals

## License

MIT License
