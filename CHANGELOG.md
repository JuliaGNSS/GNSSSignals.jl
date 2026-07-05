# Changelog

# [3.2.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v3.1.0...v3.2.0) (2026-07-05)


### Features

* introduce per-constellation signal supertypes (`AbstractGPSSignal`, `AbstractGalileoSignal`); collapse `get_time_system` to one method per GNSS ([0d4f23f](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/0d4f23f))

# [3.1.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v3.0.0...v3.1.0) (2026-07-04)


### Features

* expose signal constants at the type level + TimeSystem abstraction ([d9813a2](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/d9813a20bd55dcb29f019fc8e5f103c634a6d2d3))

# [3.0.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.6.0...v3.0.0) (2026-07-02)


* feat(code_lut)!: embedded-LUT gen_code!/code_engine; remove CodeReplicaLUT + legacy generator ([0b92291](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/0b922916c013aa0f1046b74d6dfa5ef3f86389b2))


### Bug Fixes

* avoid varargs-splat allocation in CodeGeneratorLUT4 iterate ([47c4130](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/47c4130e3893fa77e04e6f8b03fc88b96c9b1b6a))
* **bench:** make the LUT continuing-fill rows API-agnostic across revisions ([7610aa3](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/7610aa3710df221b34d241e1e0e7fbfe3544c4a8)), closes [#69](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/69)
* **code_lut:** [@static-gate](https://github.com/static-gate) the AVX-512 branch so it never compiles on non-VBMI hosts ([8451496](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/84514962db5c27ee5a62e5c8f0cb72832a2185fd)), closes [#104](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/104)
* **code_lut:** apply non-baked secondary on the DDA-correct sample ([087a936](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/087a936a48259636bdc05e95b5e6fff2782febfb))
* **code_lut:** AVX2 runtime feature re-check + dedup default_backend to _select_backend ([5a186f4](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/5a186f45602906b59bd088c369a6f846aa4d9bd7)), closes [#104](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/104) [#104](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/104) [#104](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/104) [#126](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/126) [#129](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/129)
* **code_lut:** dense-target fallback in permute fill_continue! (strided Int8 views) ([abfe3c5](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/abfe3c56da9f93a1eb9920a1f42cd9ba93cc0f07)), closes [#103](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/103) [#124](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/124)
* **code_lut:** derive CBOC Int8 amplitudes from boc1_power ([9289f60](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/9289f608c410d360fd88e8b1a65a32e08520ac8e)), closes [#112](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/112)
* **code_lut:** keep a steady (N-less) _use_runfill for the continuing generator ([2f2484c](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/2f2484c3a4cef81d45931d788877ae2e14a1437a))
* **code_lut:** keep the continuing run-fill generator 0-allocation ([e387252](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/e387252fb87d9d41c7701c9df71453980743c48e))
* **code_lut:** re-validate CPU features at runtime to avoid AVX-512 SIGILL ([b651399](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/b651399366ca90bfd660c8d654a5ac03777b77a1)), closes [#104](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/104)
* **code_lut:** reject step_denominator != 2^30 loudly ([b58857d](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/b58857d42d5936c0ad469655bdb647243915f3f5)), closes [#109](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/109)
* **code_lut:** require Unitful.Frequency at public entry points ([24ac970](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/24ac9705a081e17e6fa0b95ff86e2eb63915b055)), closes [#105](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/105)
* **code_lut:** route strided Int8 targets through an indexed fallback ([b2a1290](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/b2a1290a96d4e294c53ddb76223775018839f69f)), closes [#103](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/103)
* **code_lut:** validate K in code_engine(Val(K)) and widen _phase_type margin ([94a8ba6](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/94a8ba60bb8369603ca2efce1ad3742ded4c964a)), closes [#128](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/128) [#128](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/128)
* **code_lut:** widen _make_engine value-engine strides to Int64 (32-bit) ([4492022](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/4492022b9ecda4fbc03d2595e70bb7b7fff30172)), closes [#108](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/108) [#125](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/125)
* **code_lut:** widen/re-anchor secondary-code period-walk (Int64 hang + 32-bit OOB) ([394e314](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/394e314edbf8243ed112826f4449df606058bb6c)), closes [#102](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/102) [#108](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/108) [#102](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/102) [#108](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/108) [#102](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/102) [#108](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/108)


### Features

* **code_lut:** add SignalLUT, build_signal_lut, embed lut in signal structs ([eb9121f](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/eb9121fad58e7442da856ad89982115425044be7))
* **code_lut:** bake cosine-phased BOC into the embedded LUT ([21c2e4e](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/21c2e4e5c6a9a64d3e1cb94d873833754c74cfba)), closes [#106](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/106) [#106](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/106)
* **code_lut:** embed SignalLUT in the signals added on master ([c08a6f7](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/c08a6f776442fee57511f9e16b1c36b2bba2cba3))
* **code_lut:** replace run-fill with an exact boundary-fill kernel ([5d25556](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/5d2555600dbff3920653ca90736f74af89bb2075))
* **code_lut:** support fractional sub-chip start phase ([f04bbd5](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/f04bbd5de5795bdeb36de34ec5b0f206ced37984))
* **code_lut:** support Galileo E1B CBOC via Int8 integer-amplitude LUT ([a1d22db](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/a1d22dba0c0cbf6cef96f097f6f27b413fd3af19)), closes [hi#oversampling](https://github.com/hi/issues/oversampling)
* **code_lut:** value-based code engine/state API; remove iterator approach ([37407e8](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/37407e8e4c16416639cddd420b4735d68e6aa49c))
* CodeGeneratorLUT4 — 4-wide code iterator (pairs with CarrierIterator4) ([66cbef1](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/66cbef1dcefa0e51ef0e71af1f4938e37d89d7ee))
* NEON tbl1 code-gen backend for Apple Silicon ([1761331](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/1761331e37e93494a1db07caa94cf62fc59b7898))
* SIMD code resampler (CodeLUT) with gen_code! plan dispatch ([a433a87](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/a433a87deb5199907f5675fda550dbe6b376e6cd))


### Performance Improvements

* **code_lut:** accumulate the boundary product instead of a mulx per chip ([0a02a86](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/0a02a86bc8cdd7e0e9ee5902806d2e8650a94cda))
* **code_lut:** add SW=8 store rung + 64x/128x benchmark rows ([82fee63](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/82fee63cd017fdf38f78dc10cbdfd81b2ef08883))
* **code_lut:** broadcast run-fill for high-oversampling code generation ([a7974b5](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/a7974b598adcd64b1ccc5bdf74c175d110c9f3ae)), closes [hi#oversampling](https://github.com/hi/issues/oversampling)
* **code_lut:** const-fold the AVX2/Portable backend branches (no runtime Ref read) ([572b624](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/572b6247a70627672cebf3c6aef8f546e73ebe3e)), closes [#104](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/104)
* **code_lut:** drop the AVX-512 base-advance idiv (mod -> conditional subtract) ([d5f8394](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/d5f8394be7001b881468801afa8747c124b00850))
* **code_lut:** fix run-fill crossover firing 1–2× oversampling too late ([a3bec4f](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/a3bec4f9abc2d2568bb50beccc6fbe8e6115a532))
* **code_lut:** ISA-select the boundary step; AVX2 gate back to 3 ([71ec18e](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/71ec18ee7f8d9ece3465e2a3cf78c4f631dbb07a))
* **code_lut:** lower AVX-512 run-fill threshold 8→7 (measured crossover) ([39fcb37](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/39fcb37a71bf30372d07bc7fbe3858af2f7e380b))
* **code_lut:** lower AVX2 run-fill threshold 4→3 (measured crossover) ([b173bcc](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/b173bcce15ac70c0d30adbcf551d6b885ec06789))
* **code_lut:** lower NEON run-fill threshold to m=3 (CI-tuned) ([0e5dff0](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/0e5dff0feb838ab2702ea36b36fb6bdd2c64a9b6))
* **code_lut:** make AVX-512 run-fill threshold N-aware (lower for short fills) ([cb0e2d2](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/cb0e2d20a95d246a1a1256b5a2d5d886e26c2a6b))
* **code_lut:** make the one-shot run-fill path 0-allocation ([b556137](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/b556137c2d94a0f304f18db81449bc8bf7e19190))
* **code_lut:** parameterize padded field to avoid copying the LUT column ([1a37a46](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/1a37a46f595c2797f85d40d6244c871bac5c207b)), closes [#107](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/107)
* **code_lut:** retune run-fill inner-store padding for the Int8 buffer ([87d5d5b](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/87d5d5baa35cd0bf0c1df244bacb46c741073481))
* **code_lut:** revert run-fill thresholds to [#69](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/69) values (avoid regressions) ([01acb00](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/01acb00a611cbc50a52b25518b3b84304b781243)), closes [#69-validated](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/69-validated)
* **code_lut:** round run-fill store width up to a power of two ([a20a0ff](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/a20a0ff451c45b933175af3a55b0f73499ac6db4))
* **code_lut:** split-constant recompute permute kernels (all backends) ([d767802](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/d767802ffa7df17c915344a5f1f03b0d716d4733))
* eliminate LUT code-generation construction allocations ([7aac2bf](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/7aac2bfbab2fdab974c77e7498b8bdef7f2212e9))
* single-stream SIMD tail in the array code kernels ([51c9f02](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/51c9f024f365467d3be5c32516f9fa41869db1ef))
* vectorize the AVX2 code-gen init (263 ns -> 33 ns/call) ([cd544c6](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/cd544c6c946df62907ef795cf977d8b423670e66))


### Reverts

* Revert "perf(code_lut): fix run-fill crossover firing 1–2× oversampling too late" ([1b6e156](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/1b6e156d1a620869e4a24dae36be24cda82461f5))
* Revert "perf(code_lut): retune run-fill inner-store padding for the Int8 buffer" ([41de19f](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/41de19fc64c4813f57176d016604a498417f56fd))


### BREAKING CHANGES

* the public sampled-code API changes element type and drops a
type. `gen_code!` now fills `AbstractVector{Int8}` only — it no longer accepts
the `Int16`/`Float32` buffers that `get_code_type` used to select (including the
`Float32` buffer previously required for CBOC signals) — and
`gen_code(n, signal, prn, …)` returns a `Vector{Int8}` instead of a
`Vector{get_code_type(signal)}`. Callers that pre-allocate their own output
buffer must change its element type to `Int8`, and callers that relied on the
`gen_code` return eltype must update accordingly. The `CodeReplicaLUT` plan type
is removed and no longer exported; continuing/tracking generation moves from the
plan/iterator objects to the value-threaded
`code_engine`/`code_state`/`gen_code!(out, eng, st)` API. That is why this is
breaking: existing code that allocated non-Int8 buffers, depended on the old
return type, or referenced `CodeReplicaLUT` will no longer compile or run
unchanged.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>

# [2.6.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.5.0...v2.6.0) (2026-06-30)


### Features

* add GPS L2C civil signals (L2 CM and L2 CL) ([9cd0e2a](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/9cd0e2a96af28e688a25ee3a6f29a29d52dbb758))

# [2.5.0](https://github.com/JuliaGNSS/GNSSSignals.jl/compare/v2.4.0...v2.5.0) (2026-06-30)


### Features

* add GPS L5-Q pilot signal ([21c1e1b](https://github.com/JuliaGNSS/GNSSSignals.jl/commit/21c1e1b73944eb10a84776397fb9d6b686b6ee60))

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
