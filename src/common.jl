"""
$(SIGNATURES)

Get the full code matrix for a GNSS signal.

Returns the codes as a matrix where each column represents a PRN.

# Arguments
- `signal`: A GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B()`)

# Returns
- `Matrix`: Code matrix of size `(code_length, num_prns)`

# Examples
```julia-repl
julia> codes = get_codes(GPSL1CA())
julia> size(codes)
(1023, 37)
```
"""
function get_codes(signal::AbstractGNSSSignal)
    signal.codes
end

"""
$(SIGNATURES)

Widen a primary-code matrix from its on-disk / LFSR-generated `Int8`
representation to `Int16`.

Chip values are ±1 and would fit in `Int8`, but on x86_64 / AVX2 hardware
storing as `Int16` is materially faster for `gen_code!`: the inner store
loop emits a clean `vpbroadcastw` + `vmovq` pattern, while Int8 storage
triggers an `shl 8 + or` byte-packing antipattern (3 extra μops per chip)
when the buffer is also `Int8`. Storing Int8 chips into an Int16 buffer
recovers the clean codegen but still loses ~14 % because the `movsx`
load chain runs slower than `movzx`.

Probed alternatives that did **not** recover Int8 perf on AVX2:
- `@simd ivdep` annotation on the inner store loop (no measurable
  effect — LLVM has already made its codegen choice for `Val`-known
  fixed-trip loops)
- `unsafe_store!` of a replicated `UInt32` to bypass LLVM's
  pattern-matcher (helped Int8/Int8 by ~5 % but still ~12 % slower
  than Int16/Int16, and a regression for Int16/Int16)
- `SIMD.jl` `Vec{N,T}` broadcast-store (slower at every NUM_INNER
  measured because the abstraction forces an xmm materialization for
  small N)

Memory cost of widening is small in absolute terms — the largest
current matrix (GPS L5-I, 10230 × 37) is 757 KB at Int16. The trade-off
may invert on AVX-512 or non-x86 hardware; revisit if you have access
to such platforms.
"""
widen_codes_to_storage(codes::AbstractMatrix) = Int16.(codes)

"""
$(SIGNATURES)

Get the `Symbol` identifier of a GNSS signal — the machine-readable per-signal
key, e.g. `:GPSL1CA`, `:GalileoE1B`, `:GalileoE5aI`.

This is the finest identity level: the same PRN measured on two bands of one
constellation (e.g. GPS L1 and L5) has two distinct signal ids.

Defaults to `nameof` of the signal type; override `get_signal_id(::Type{MySignal})`
to pin a specific symbol (e.g. so an approximation type reports the id of the
signal it approximates). Works on an instance or a type and folds to a
compile-time constant.

Distinct from [`get_signal_name`](@ref), which is a human-readable display
string (`"GPS L1 C/A"`), and coarser than [`get_band_id`](@ref) (`:L1`).

```julia-repl
julia> get_signal_id(GPSL1CA())
:GPSL1CA

julia> get_signal_id(GalileoE1B)
:GalileoE1B
```
"""
@inline get_signal_id(::Type{S}) where {S<:AbstractGNSSSignal} = nameof(S)
@inline get_signal_id(s::AbstractGNSSSignal) = get_signal_id(typeof(s))

"""
$(SIGNATURES)

Get the `Symbol` identifier of the GNSS constellation a signal belongs to —
`:GPS` or `:Galileo`.

This is the coarsest identity level: every signal of one constellation, across
all bands and PRNs, maps to the same id (GPS L1 C/A, GPS L5I and GPS L1C are all
`:GPS`). It's the natural key for logging, dictionaries and constellation-level
branching; for a type-level test prefer `signal isa AbstractGPSSignal`.

Defined once per constellation on the abstract signal type
(`get_constellation_id(::Type{<:AbstractGPSSignal}) = :GPS`), so every concrete
signal inherits it without constructing a value. Works on an instance or a type
and folds to a compile-time constant.

Distinct from [`get_time_system`](@ref) (the constellation's reference time
scale, `GPST()`/`GST()`), and coarser than [`get_band_id`](@ref) (`:L1`) and
[`get_signal_id`](@ref) (`:GPSL1CA`).

```julia-repl
julia> get_constellation_id(GPSL1CA())
:GPS

julia> get_constellation_id(GalileoE1B)
:Galileo
```
"""
@inline get_constellation_id(::Type{<:AbstractGPSSignal}) = :GPS
@inline get_constellation_id(::Type{<:AbstractGalileoSignal}) = :Galileo
@inline get_constellation_id(s::AbstractGNSSSignal) = get_constellation_id(typeof(s))

# Each concrete signal defines these on `::Type{<:Signal}` (per-signal files);
# these forward an instance to its type.
@inline get_code_length(s::AbstractGNSSSignal) = get_code_length(typeof(s))
@inline get_code_frequency(s::AbstractGNSSSignal) = get_code_frequency(typeof(s))
@inline get_data_frequency(s::AbstractGNSSSignal) = get_data_frequency(typeof(s))

# NOTE: the legacy fixed-point `gen_code!` and its `sample_code!` / `dispatch_sample_code_worker!`
# / `sample_code_worker!` / `sample_code_worker_generic!` / `sample_code_tail!` /
# `_pad_inner_iterations` machinery (plus the `SAMPLE_CODE_INNER_THRESHOLD` / `HAS_AVX512`
# constants), and the `_select_codes_for` secondary-matrix selector (with its GPSL5I/GPSL1C_P
# specializations and their cached `negated_codes` fields), have all been retired: the embedded
# SIMD LUT in `code_lut.jl` is now THE `gen_code!`. The subcarrier multiply
# (`multiply_with_subcarrier!`, the `calc_subcarrier_*` helpers, and the `_tmboc_*` /
# `_pack_tmboc_pattern` TMBOC kernels) has now been retired too — it had no production caller
# once the LUT became `gen_code!`. What remains here is just `get_codes` /
# `widen_codes_to_storage`, the allocating `gen_code` wrapper, and the `get_code_spectrum_*`
# helpers. `get_code` / `get_code_unsafe` (modulation.jl) never used any of the deleted helpers.

"""
$(SIGNATURES)

Generate a sampled code signal for a given PRN.

Allocates a new `Vector{Int8}` and returns the spreading code sampled at the specified
sampling frequency. Thin allocating wrapper over the in-place [`gen_code!`](@ref); see it for
the Int8 output semantics (±1, or the CBOC integer approximation), fractional-phase support,
and the `sampling_frequency ≥ code_frequency · subchip_factor` requirement.

# Arguments
- `num_samples`: Number of samples to generate
- `signal`: GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B()`)
- `prn`: PRN number of the satellite
- `sampling_frequency`: Sampling frequency (must be ≥ code frequency · subchip_factor)
- `code_frequency`: Code chipping rate (default: signal's nominal code frequency)
- `start_phase`: Initial code phase in chips (default: 0.0)
- `start_index`: Index offset (default: 0)

# Returns
- `Vector{Int8}`: Sampled code signal

# Examples
```julia-repl
julia> using Unitful: MHz
julia> sampled_code = gen_code(4000, GPSL1CA(), 1, 4MHz)
julia> length(sampled_code)
4000
```
"""
function gen_code(
    num_samples::Integer,
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(signal),
    start_phase = 0.0,
    start_index::Integer = 0,
)
    code = zeros(Int8, num_samples)
    gen_code!(code, signal, prn, sampling_frequency, code_frequency, start_phase, start_index)
end

"""
$(SIGNATURES)

Get the ratio of code frequency to center frequency.

This ratio is used to compute the code Doppler from the carrier Doppler.

# Arguments
- `signal`: A GNSS signal instance

# Returns
- `Float64`: The code-to-center frequency ratio

# Examples
```julia-repl
julia> get_code_center_frequency_ratio(GPSL1CA())
0.0006493506493506494
```
"""
@inline function get_code_center_frequency_ratio(signal::AbstractGNSSSignal)
    get_code_frequency(signal) / get_center_frequency(signal)
end

"""
$(SIGNATURES)

Get the minimum number of bits needed to represent the code length.

Calculates the number of bits required to represent the full code length,
including secondary code if present.

# Arguments
- `signal`: A GNSS signal instance

# Returns
- `Int`: Number of bits needed

# Examples
```julia-repl
julia> min_bits_for_code_length(GPSL1CA())
10
julia> min_bits_for_code_length(GPSL5I())
17
```
"""
@inline function min_bits_for_code_length(signal::AbstractGNSSSignal)
    ndigits(get_code_length(signal) * get_secondary_code_length(signal); base = 2)
end

"""
$(SIGNATURES)

Get the length of the secondary code.

# Arguments
- `signal`: A GNSS signal instance

# Returns
- `Int`: Secondary code length (1 if no secondary code)

# Examples
```julia-repl
julia> get_secondary_code_length(GPSL1CA())
1
julia> get_secondary_code_length(GPSL5I())
10
```
"""
@inline function get_secondary_code_length(signal::AbstractGNSSSignal)
    secondary_code_length(get_secondary_code(signal))
end

"""
$(SIGNATURES)

Calculate the spectral power density of a BPSK modulated signal.

Computes the power spectral density at baseband frequency `f` for a BPSK
signal with chip rate `fc`.

# Arguments
- `fc`: Code chip rate
- `f`: Baseband frequency at which to evaluate the spectrum

# Returns
- Spectral power density value

# Examples
```julia-repl
julia> using Unitful: MHz, kHz
julia> get_code_spectrum_BPSK(1.023MHz, 0kHz)
9.775171065493646e-7
```
"""
function get_code_spectrum_BPSK(fc::Frequency, f)
    return get_code_spectrum_BPSK(fc / 1Hz, f)
end
function get_code_spectrum_BPSK(fc, f::Frequency)
    return get_code_spectrum_BPSK(fc, f / 1Hz)
end
function get_code_spectrum_BPSK(fc::Frequency, f::Frequency)
    return get_code_spectrum_BPSK(fc / 1Hz, f / 1Hz)
end
function get_code_spectrum_BPSK(fc, f)
    return sinc(f / fc)^2 / fc
end

"""
$(SIGNATURES)

Calculate the spectral power density of a sine-phased BOC modulated signal.

Computes the power spectral density at baseband frequency `f` for a BOC(sin)
signal with chip rate `fc` and subcarrier frequency `fs`.

# Arguments
- `fc`: Code chip rate
- `fs`: Subcarrier frequency
- `f`: Baseband frequency at which to evaluate the spectrum

# Returns
- Spectral power density value
"""
function get_code_spectrum_BOCsin(fc::Frequency, fs::Frequency, f)
    return get_code_spectrum_BOCsin(fc / 1Hz, fs / 1Hz, f)
end
function get_code_spectrum_BOCsin(fc, fs, f::Frequency)
    return get_code_spectrum_BOCsin(fc, fs, f / 1Hz)
end
function get_code_spectrum_BOCsin(fc::Frequency, fs::Frequency, f::Frequency)
    return get_code_spectrum_BOCsin(fc / 1Hz, fs / 1Hz, f / 1Hz)
end
function get_code_spectrum_BOCsin(fc, fs, f)
    return ((sinc(f / fc) * tan(pi * f / (2 * fs)))^2 / fc)
end

"""
$(SIGNATURES)

Calculate the spectral power density of a cosine-phased BOC modulated signal.

Computes the power spectral density at baseband frequency `f` for a BOC(cos)
signal with chip rate `fc` and subcarrier frequency `fs`.

# Arguments
- `fc`: Code chip rate
- `fs`: Subcarrier frequency
- `f`: Baseband frequency at which to evaluate the spectrum

# Returns
- Spectral power density value
"""
function get_code_spectrum_BOCcos(fc::Frequency, fs::Frequency, f)
    return get_code_spectrum_BOCcos(fc / 1Hz, fs / 1Hz, f)
end
function get_code_spectrum_BOCcos(fc, fs, f::Frequency)
    return get_code_spectrum_BOCcos(fc, fs, f / 1Hz)
end
function get_code_spectrum_BOCcos(fc::Frequency, fs::Frequency, f::Frequency)
    return get_code_spectrum_BOCcos(fc / 1Hz, fs / 1Hz, f / 1Hz)
end
function get_code_spectrum_BOCcos(fc, fs, f)
    return (2 * sinc(f / fc) * sinpi(f / 4fs)^2 / cospi(f / 2fs))^2 / fc
end
