# Changelog

# [2.4.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.3.0...v2.4.0) (2026-06-30)


### Features

* add Galileo E1C (pilot component of Galileo E1 OS) ([55210b9](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/55210b99d026da2a68b91e28ab2af491f8264937))

# [2.3.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.2.2...v2.3.0) (2026-06-29)


### Features

* add Galileo E5a signals (E5a-I and E5a-Q) ([23e22a6](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/23e22a60346b05696ad98d2c4fafd1c60fcaff5a))

## [2.2.2](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.2.1...v2.2.2) (2026-06-08)


### Bug Fixes

* **common:** prevent Int64 overflow in sample_code! for very short outputs ([8f472e6](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/8f472e67c75b06d40bd618f2a280bbbcd27a5967)), closes [#63](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/63)

## [2.2.1](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.2.0...v2.2.1) (2026-05-20)


### Bug Fixes

* **l1c_d:** report channel symbol rate (100 Hz), not info bit rate (50 Hz) ([0fee3ab](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/0fee3ab94e3b3517b19f39ea64553335f5c0e801))

# [2.2.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.1.1...v2.2.0) (2026-05-19)


### Features

* add GalileoE1B_BOC11 (BOC(1,1) approximation of E1B) ([141922e](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/141922e1e3dfdf77c2a0ee082852e14397a9b1c1))

## [2.1.1](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.1.0...v2.1.1) (2026-05-19)


### Bug Fixes

* **cboc:** integer buffer raises ArgumentError with hint, not InexactError ([2ee607f](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/2ee607f3684defe03f6223b6e67dd9e1789677a6))


### Performance Improvements

* **tmboc:** route contiguous-Int16 views through SIMD fast path ([6e3b7e8](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/6e3b7e80391900015d9839d328b84593980891ba))

# [2.1.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.0.0...v2.1.0) (2026-05-15)


### Bug Fixes

* **boc:** align BOC sub-carrier with primary chip transitions at chip-aligned rates ([47de944](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/47de9443ad5f218c1ed14557ef99781323c4adf2))
* **sample_code:** apply secondary_start_index in the tail loop ([fe5aaa1](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/fe5aaa1a6001bbee8b3f28a569121d1ede4be02a))


### Features

* add GPS L1C-D and L1C-P signals ([81581d1](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/81581d1f629efc152b919692881069171036c6b1))


### Performance Improvements

* **tmboc:** explicit 16-lane SIMD.jl path for Int16 and Float32 ([03cd623](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/03cd623492b2f0cbf6c6b3100cc9b04431040369))
* **tmboc:** restore fast path for single-transition 16-lane blocks ([d6b0ebf](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/d6b0ebf078f4cdad2c30cf82a7dafbd7491500c0))
* **tmboc:** two-pass BOC1 + selective BOC2 fix-up ([94c9e0a](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/94c9e0af6e02433ef6216c222898b5b063d75f5b))

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
