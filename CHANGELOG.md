# Changelog

# [2.0.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v1.0.3...v2.0.0) (2026-05-12)


* refactor!: model secondary codes as their own type hierarchy ([d572f67](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/d572f6751c2d1725ffe6ac53466dd4e8a1be932a))
* refactor!: rename signal types and factor out Band abstraction ([66aa1e3](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/66aa1e34b5708fbfc7d06b8791ea6ba60a6a49da))


### Bug Fixes

* address Aqua unbound type parameter and doc [@ref](https://github.com/ref) errors ([731961b](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/731961be3f56390a6e8f03d9f18ec8e0dfad9203))


### Performance Improvements

* extend inner-loop padding ladder to {4, 8, 16} ([15a4174](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/15a4174bf905e003cf8ec41925c78ed2d94b1c9c))
* hoist the per-chip secondary multiply out of the L5-I hot loop ([0451807](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/0451807c7794477f5f3369cb098d453c090ac4ae))
* pad inner store loop to ≥4 stores per chip ([ed1d199](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/ed1d1996a2d3e5b29ef465c3960c1ebd0f023ee8))
* per-arch inner-loop padding ladder + dedup tail helper ([a25cc32](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/a25cc328b8a44608e69cb28ef1271ec315236853))
* skip padding for real_num_inner == 9 ([3f1df4b](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/3f1df4bd426f53d3903b55552ce8bdb3a791c9d2))


### BREAKING CHANGES

* `get_secondary_code(::AbstractGNSSSignal)` returns a
`SecondaryCode` instance, not `Integer` or `Tuple`. The 2-arg
`get_secondary_code(signal, phase)` and the `Integer`/`Tuple`
dispatch helpers are gone.
*   - AbstractGNSS -> AbstractGNSSSignal
  - GPSL1 -> GPSL1CA
  - GPSL5 -> GPSL5I
  - get_system_string -> get_signal_name (returns a human-readable name
    like "GPS L1 C/A", not the Julia type name)
  - get_center_frequency is no longer defined per concrete signal;
    concrete signals expose get_band(::AbstractGNSSSignal) :: Band and
    get_center_frequency dispatches through it
  - data/codes_gps_l1.bin renamed to data/codes_gps_l1ca.bin

No deprecation shims — downstream callers update names in one pass.

## [1.0.3](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v1.0.2...v1.0.3) (2026-04-27)


### Performance Improvements

* drop Val{MESF}/Val{MED} args and vectorize CBOC subcarrier ([e4ab80f](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/e4ab80f5fd33d468e314ff887d9338021f499286))

## [1.0.2](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v1.0.1...v1.0.2) (2026-02-03)


### Bug Fixes

* add instance methods for get_modulation to fix type instability ([76b3f13](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/76b3f13191e08f4e48ac33764f26d11c172ff5fe))

## [1.0.1](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v1.0.0...v1.0.1) (2026-01-11)


### Bug Fixes

* prevent integer overflow in sample_code! fixed-point arithmetic ([64a5e8d](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/64a5e8d6b29072f1ba954ec932d2b694ae9070c0))

## [0.17.4](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v0.17.3...v0.17.4) (2025-12-24)


### Bug Fixes

* correct sample assignment at chip boundaries in sample_code! ([#47](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/47)) ([d6d3b68](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/d6d3b68194d5b33fdb0779f6afe60a39b05c4372))
