[![pipeline status](https://git.rwth-aachen.de/nav/GNSSSignals.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/GNSSSignals.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/GNSSSignals.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/GNSSSignals.jl/commits/master)

# Generate GNSS signals.

## Features

* GPS L1

## Getting started

Install:
```julia
Pkg.clone("git@git.rwth-aachen.de:nav/GNSSSignals.jl.git")
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

* add GPS Q5 of the L5 signal
* Galileo signals

## License

MIT License
