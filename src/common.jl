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

# Maximum num_inner_iterations for which we generate a Val-specialized variant.
# Covers oversampling ratios up to 64, which is enough for virtually all GNSS
# receiver sampling rates. Above this, a @simd fallback is used; it remains
# within a few percent of the specialized version.
const SAMPLE_CODE_INNER_THRESHOLD = 64

# Runtime-detected at module load. `+avx512f` is the base AVX-512 feature
# (foundation); every other AVX-512 subset implies it. Used by
# `_pad_inner_iterations` to pick an arch-tuned padding ladder.
const HAS_AVX512 = let
    try
        f = ccall(:jl_get_cpu_features, Any, ())
        f isa AbstractString && occursin("+avx512f", f)
    catch
        false
    end
end

# Pad the inner store loop to a length that LLVM emits a clean wide SIMD
# store for. Per-chip "extra" writes get overwritten by the next chip's
# writes (overwrite-tolerance), so over-padding only costs slightly extra
# store bandwidth. The very last main-loop chip's extras must fit in the
# buffer, so the caller holds back `num_inner_iterations - real_num_inner`
# chips for the tail loop to handle with bounds-respecting writes.
#
# Two ladders, picked per-arch from the runtime feature string:
#
# - AVX2: ladder steps at {4, 8, 16}, with a `real == 9` exception
#   (pad=16 there loses ~19% on Zen 3 / Intel client AVX2). Tuned on
#   Zen 3 (this laptop's reference). Confirmed: real ∈ {10..16} all
#   benefit from pad=16; real=9 stays unpadded.
#
# - AVX-512: ladder steps at {4, 8, 12, 16, 20, 24}. Additions vs AVX2:
#     * pad to 12 for real ∈ {10, 11, 12}: beats pad=16 by 13–23% on
#       Zen 5.
#     * real ∈ {17, 18}: no padding — both regress >20 % at pad=24
#       and >8 % at pad=20.
#     * pad to 20 for real == 19: −15 % vs natural; pad=24 loses here.
#     * pad to 24 for real ∈ {21, 22, 23}: beats natural by 9–24 %.
#   The `real == 9` exception is preserved (same outlier on both
#   archs). Above 24, no padding pays.
#
# Other architectures (ARM NEON, AVX-512 with VBMI2 tweaks, Apple
# Silicon) fall through to the AVX2 ladder by default. If you re-tune
# for another arch, sweep `real ∈ 2..32` over `benchmark/benchmarks.jl`
# and add a new branch here.
@inline function _pad_inner_iterations(real_num_inner::Int)
    if HAS_AVX512
        real_num_inner <= 4  ? 4  :
        real_num_inner <= 8  ? 8  :
        real_num_inner == 9  ? 9  :
        real_num_inner <= 12 ? 12 :
        real_num_inner <= 16 ? 16 :
        real_num_inner <= 18 ? real_num_inner :
        real_num_inner <= 20 ? 20 :
        real_num_inner <= 23 ? 24 : real_num_inner
    else
        real_num_inner <= 4  ? 4  :
        real_num_inner <= 8  ? 8  :
        real_num_inner == 9  ? 9  :
        real_num_inner <= 16 ? 16 : real_num_inner
    end
end

"""
$(SIGNATURES)

Generate the code signal in-place for a given PRN.

This is a highly optimized function for generating sampled spreading codes at arbitrary
sampling rates. It uses fixed-point arithmetic and minimizes memory access by exploiting
the fact that consecutive samples often map to the same code chip. The algorithm avoids
per-sample floating-point operations and modulo calculations, making it suitable for
real-time GNSS signal processing applications.

Samples the spreading code at the specified sampling frequency and stores the result
in the provided buffer. Includes subcarrier modulation for BOC-type signals.

# Arguments
- `sampled_code`: Pre-allocated output buffer
- `signal`: GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B()`)
- `prn`: PRN number of the satellite
- `sampling_frequency`: Sampling frequency (must be larger than code frequency)
- `code_frequency`: Code chipping rate (default: signal's nominal code frequency)
- `start_phase`: Initial code phase in chips (default: 0.0)
- `start_index_shift`: Index offset for the output buffer (default: 0)
- `PHASET`: Integer type for phase calculations (default: `Int32`)

# Returns
- The modified `sampled_code` buffer

# Examples
```julia-repl
julia> using Unitful: MHz
julia> buffer = zeros(Int16, 4000)
julia> gen_code!(buffer, GPSL1CA(), 1, 4MHz)
```
"""
function gen_code!(
    sampled_code::AbstractVector,
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(signal),
    start_phase = 0.0,
    start_index_shift::Integer = 0,
    PHASET = Int32,
)
    sample_code!(
        sampled_code,
        signal,
        prn,
        sampling_frequency,
        code_frequency,
        start_phase,
        start_index_shift,
    )
    multiply_with_subcarrier!(
        sampled_code,
        get_modulation(signal),
        sampling_frequency,
        code_frequency,
        start_phase,
        start_index_shift,
        PHASET,
    )
    sampled_code
end

function sample_code!(
    sampled_code::AbstractVector,
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(signal),
    start_phase = 0.0,
    start_index_shift::Integer = 0,
)
    modulated_code_frequency = code_frequency * get_code_factor(signal)
    frequency_ratio = sampling_frequency / modulated_code_frequency
    modulated_code_frequency > sampling_frequency && error(
        "The sampling frequency must be larger than the code frequency multiplied by code factor (larger than $modulated_code_frequency, it is $sampling_frequency).",
    )

    # The -2 (instead of -1) reserves an extra bit for overflow headroom.
    # delta_sum accumulates (num_samples * frequency_ratio) in fixed-point representation.
    # With -1, delta_sum uses nearly the full Int range and can overflow when combined
    # with the (1 << fixed_point) offset in the initial delta_sum calculation.
    fixed_point = sizeof(Int) * 8 - 2 - ndigits(length(sampled_code); base = 2)

    frequency_ratio_fixed_point = round(Int, frequency_ratio * 1 << fixed_point)

    start_phase_including_shift =
        start_phase + start_index_shift * code_frequency / sampling_frequency

    sec = get_secondary_code(signal)
    primary_length = get_code_length(signal)
    secondary_length = secondary_code_length(sec)
    full_cycle_length = primary_length * secondary_length
    # Compute fractional part before any normalization to preserve exact floating-point value
    floored_code_phase = floor(Int, start_phase_including_shift)
    frac_part = floored_code_phase - start_phase_including_shift  # Always in (-1, 0]
    # Split the absolute phase into (primary chip index, secondary chip index)
    # within one full secondary cycle.
    absolute_chip = mod(floored_code_phase, full_cycle_length)
    code_start_index = mod(absolute_chip, primary_length)
    secondary_start_index = div(absolute_chip, primary_length)
    # The -256 offset handles boundary cases where a sample falls exactly on a chip
    # boundary (e.g., phase = 2.0). Without this offset, floor(2.0) = 2 would assign
    # two samples to the first chip when only one belongs there. The offset ensures
    # samples exactly on boundaries are assigned to the next chip.
    delta_sum =
        floor(Int, frac_part * frequency_ratio_fixed_point) + (1 << fixed_point) - 256
    raw_num_code_samples =
        Int(fld(modulated_code_frequency * length(sampled_code), sampling_frequency))
    real_num_inner = ceil(Int, frequency_ratio)
    num_inner_iterations = _pad_inner_iterations(real_num_inner)
    tail_slack = num_inner_iterations - real_num_inner
    num_code_samples_to_iterate = max(0, raw_num_code_samples - tail_slack)
    dispatch_sample_code_worker!(
        sampled_code,
        signal,
        sec,
        prn,
        frequency_ratio_fixed_point,
        fixed_point,
        code_start_index,
        secondary_start_index,
        delta_sum,
        num_code_samples_to_iterate,
        primary_length,
        secondary_length,
        num_inner_iterations,
        tail_slack,
    )
    return sampled_code
end

"""
$(SIGNATURES)

Select which code matrix the inner store loop should read from, given
the current primary-period's secondary-chip value.

Returns a `(codes, multiplier)` tuple. The inner loop reads
`codes[i, prn] * multiplier` per chip.

The default implementation returns `(signal.codes, sec_val)` — the inner
loop multiplies each chip by the secondary value as usual. Per-signal
specializations can pre-negate the codes matrix at construction time
and return `(positive_or_negated_codes, true)`; the `* true` then
elides at compile time, hoisting the per-chip multiply out of the hot
loop. See [`GPSL5I`](@ref) for an example.

This is an internal helper for the `sample_code_worker!` /
`sample_code_worker_generic!` hot path. `get_code` does its own
straightforward lookup and does not use this dispatch.
"""
@inline _select_codes_for(signal::AbstractGNSSSignal, sec_val) = (signal.codes, sec_val)

# Shared tail loop for both `sample_code_worker!` and
# `sample_code_worker_generic!`. The tail handles the chips held back
# from the main loop (`tail_slack`) plus a 3-chip safety margin; each
# iteration clamps writes to the buffer end via `min(...)` and breaks
# once `prev` has filled the buffer. Inlined so the workers don't pay
# a call overhead per `gen_code!`.
@inline function sample_code_tail!(
    sampled_code,
    signal,
    sec,
    prn,
    frequency_ratio_fixed_point,
    fixed_point,
    code_start_index,
    secondary_start_index,
    delta_sum,
    prev,
    processed_code_samples,
    primary_length,
    secondary_length,
    tail_slack,
)
    @inbounds for i = 0:(tail_slack + 2)
        # Offset (in chips) from the start of the buffer (which already
        # sits at primary-chip `code_start_index` of secondary-chip
        # `secondary_start_index`).
        chip_offset = processed_code_samples + code_start_index + i
        # Absolute chip index within one full secondary cycle. Must add
        # `secondary_start_index * primary_length` so the tail picks up
        # the right secondary index for buffers that fit entirely in
        # the tail loop (small N) and start past the first primary
        # period.
        absolute_chip = mod(
            secondary_start_index * primary_length + chip_offset,
            primary_length * secondary_length,
        )
        next_code_idx = mod(absolute_chip, primary_length) + 1
        sec_idx = div(absolute_chip, primary_length)
        sec_val = secondary_value(sec, prn, sec_idx)
        codes_view, code_mul = _select_codes_for(signal, sec_val)
        next_code = codes_view[next_code_idx, prn] * code_mul
        delta_sum += frequency_ratio_fixed_point
        next_prev = delta_sum >> fixed_point
        num_iterations = min(next_prev - prev, length(sampled_code) - prev)
        for j = 1:num_iterations
            sampled_code[prev+j] = next_code
        end
        prev = next_prev
        prev >= length(sampled_code) && break
    end
    return sampled_code
end

# Inner worker parameterized on `Val{NUM_INNER}` so the fixed-trip inner loop
# gets fully unrolled and vectorized by LLVM.
#
# `sec` is a SecondaryCode; for `NoSecondaryCode` the per-chip multiply
# folds to a no-op at compile time. For `SharedSecondaryCode` /
# `PerPRNSecondaryCode` we hoist the lookup outside the per-chip loop —
# the secondary value is constant within one primary-code period.
function sample_code_worker!(
    sampled_code::AbstractVector,
    signal::AbstractGNSSSignal,
    sec::SecondaryCode,
    prn::Integer,
    frequency_ratio_fixed_point::Int,
    fixed_point::Int,
    code_start_index::Int,
    secondary_start_index::Int,
    delta_sum::Int,
    num_code_samples_to_iterate::Int,
    primary_length::Int,
    secondary_length::Int,
    ::Val{NUM_INNER},
    tail_slack::Int,
) where {NUM_INNER}
    prev = 0
    num_code_iterations =
        cld(num_code_samples_to_iterate + code_start_index, primary_length)
    processed_code_samples = 0
    @inbounds for k = 0:(num_code_iterations-1)
        iteration_begin = (k == 0 ? code_start_index : 0) + 1
        iteration_end = min(
            num_code_samples_to_iterate + (k == 0 ? code_start_index : 0) -
            processed_code_samples - 1,
            primary_length,
        )
        iterations = iteration_begin:iteration_end
        processed_code_samples += length(iterations)
        sec_val = secondary_value(
            sec, prn, mod(secondary_start_index + k, secondary_length))
        # `_select_codes_for` lets signals with a ±1 SharedSecondaryCode
        # swap between a precomputed positive/negated code matrix here,
        # hoisting the per-chip multiply out of the inner loop. The
        # default returns (signal.codes, sec_val), preserving the original
        # behaviour for signals without that specialization.
        codes_view, code_mul = _select_codes_for(signal, sec_val)
        for i in iterations
            next_code = codes_view[i, prn] * code_mul
            for j = 1:NUM_INNER
                sampled_code[prev+j] = next_code
            end
            delta_sum += frequency_ratio_fixed_point
            prev = delta_sum >> fixed_point
        end
    end
    return sample_code_tail!(
        sampled_code, signal, sec, prn,
        frequency_ratio_fixed_point, fixed_point,
        code_start_index, secondary_start_index,
        delta_sum, prev, processed_code_samples,
        primary_length, secondary_length, tail_slack,
    )
end

# Fallback for oversampling ratios above SAMPLE_CODE_INNER_THRESHOLD.
#
# Kept structurally in sync with `sample_code_worker!` above — the outer
# `k` loop and the tail are identical, only the inner store differs
# (`@simd ivdep` over a runtime length here vs. a `Val`-known fixed trip
# in the specialized worker). If you fix a bug in one outer/tail
# loop, fix it in the other.
function sample_code_worker_generic!(
    sampled_code::AbstractVector,
    signal::AbstractGNSSSignal,
    sec::SecondaryCode,
    prn::Integer,
    frequency_ratio_fixed_point::Int,
    fixed_point::Int,
    code_start_index::Int,
    secondary_start_index::Int,
    delta_sum::Int,
    num_code_samples_to_iterate::Int,
    primary_length::Int,
    secondary_length::Int,
    num_inner_iterations::Int,
    tail_slack::Int,
)
    prev = 0
    num_code_iterations =
        cld(num_code_samples_to_iterate + code_start_index, primary_length)
    processed_code_samples = 0
    @inbounds for k = 0:(num_code_iterations-1)
        iteration_begin = (k == 0 ? code_start_index : 0) + 1
        iteration_end = min(
            num_code_samples_to_iterate + (k == 0 ? code_start_index : 0) -
            processed_code_samples - 1,
            primary_length,
        )
        iterations = iteration_begin:iteration_end
        processed_code_samples += length(iterations)
        sec_val = secondary_value(
            sec, prn, mod(secondary_start_index + k, secondary_length))
        codes_view, code_mul = _select_codes_for(signal, sec_val)
        for i in iterations
            next_code = codes_view[i, prn] * code_mul
            @simd ivdep for j = 1:num_inner_iterations
                sampled_code[prev+j] = next_code
            end
            delta_sum += frequency_ratio_fixed_point
            prev = delta_sum >> fixed_point
        end
    end
    return sample_code_tail!(
        sampled_code, signal, sec, prn,
        frequency_ratio_fixed_point, fixed_point,
        code_start_index, secondary_start_index,
        delta_sum, prev, processed_code_samples,
        primary_length, secondary_length, tail_slack,
    )
end

@generated function dispatch_sample_code_worker!(
    sampled_code,
    signal,
    sec,
    prn,
    frequency_ratio_fixed_point,
    fixed_point,
    code_start_index,
    secondary_start_index,
    delta_sum,
    num_code_samples_to_iterate,
    primary_length,
    secondary_length,
    num_inner_iterations,
    tail_slack,
)
    branches = [
        :(
            num_inner_iterations == $i && return sample_code_worker!(
                sampled_code,
                signal,
                sec,
                prn,
                frequency_ratio_fixed_point,
                fixed_point,
                code_start_index,
                secondary_start_index,
                delta_sum,
                num_code_samples_to_iterate,
                primary_length,
                secondary_length,
                Val($i),
                tail_slack,
            )
        # `sample_code!` pads `num_inner_iterations` to a minimum of 4, so
        # the dispatcher never sees values 1, 2, or 3.
        ) for i = 4:SAMPLE_CODE_INNER_THRESHOLD
    ]
    quote
        $(branches...)
        return sample_code_worker_generic!(
            sampled_code,
            signal,
            sec,
            prn,
            frequency_ratio_fixed_point,
            fixed_point,
            code_start_index,
            secondary_start_index,
            delta_sum,
            num_code_samples_to_iterate,
            primary_length,
            secondary_length,
            num_inner_iterations,
            tail_slack,
        )
    end
end

@inline function calc_subcarrier_phase_and_delta(
    modulation::BOC,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase,
    PHASET,
)
    subcarrier_frequency = code_frequency * 2 * modulation.m
    subcarrier_frequency > sampling_frequency && error(
        "The sampling frequency should be larger than twice the code frequency multiplied by the sub-carrier factor (larger than $subcarrier_frequency, it is $sampling_frequency).",
    )

    fixed_point = sizeof(PHASET) * 8 - 1

    subcarrier_phase =
        floor(
            unsigned(PHASET),
            mod(start_phase, 1) * unsigned(one(PHASET)) << fixed_point,
        ) << 1 * unsigned(PHASET(modulation.m))

    # Ceil-round the phase increment so that, at sampling rates where
    # the BOC half-cycle boundary lands on an exact sample boundary
    # (for example, `sampling_frequency = 12 × code_frequency` with
    # BOC(1,1), where each chip is 12 samples and the half-cycle is 6
    # samples), the integer phase accumulator reaches the sign-bit
    # threshold on the spec-aligned sample. Each BOC chip then splits
    # cleanly into two halves of equal length with the same primary
    # chip sign — the first half positive (relative to the chip),
    # then negated — matching the IS-GPS-800G / IS-GPS-200 convention
    # used by every other reference implementation we checked
    # (GNSS-SDR, PocketSDR). With `floor` the accumulator missed the
    # threshold by ≤ 1 unit in the last place at exact alignment and
    # the sign-flip landed one sample late, drifting the BOC
    # sub-carrier against the primary chip rate over the buffer.
    delta_subcarrier_phase = ceil(
        unsigned(PHASET),
        subcarrier_frequency * unsigned(one(PHASET)) << fixed_point / sampling_frequency,
    )
    fixed_point, subcarrier_phase, delta_subcarrier_phase
end

function calc_subcarrier_bit(
    i,
    fixed_point,
    subcarrier_phase,
    delta_subcarrier_phase,
    PHASET,
)
    (reinterpret(PHASET, delta_subcarrier_phase * i + subcarrier_phase) >> fixed_point) <<
    true + true
end

# Compute the subcarrier bit directly as ±1.0f0 by splicing the sign bit of the
# phase accumulator into the bit pattern of +1.0f0 (0x3F800000). This lets LLVM
# vectorize the CBOC multiply loop, which the integer-shift form blocks due to
# the reinterpret chain.
@inline function calc_subcarrier_bit_f32(i::UInt32, subcarrier_phase::UInt32, delta_subcarrier_phase::UInt32)
    phase = delta_subcarrier_phase * i + subcarrier_phase
    reinterpret(Float32, (phase & 0x80000000) | 0x3F800000)
end

# Integer version of the same trick: produce ±1 from the sign bit of the
# phase accumulator using only operations LLVM happily vectorizes
# (multiply, add, arithmetic shift, bitwise OR). Returns Int32 because
# Int32 vector ops are widely supported; downstream code converts to
# the buffer's element type.
@inline function calc_subcarrier_bit_i32(i::UInt32, subcarrier_phase::UInt32, delta_subcarrier_phase::UInt32)
    phase = delta_subcarrier_phase * i + subcarrier_phase
    # Arithmetic shift on a signed integer fills with the sign bit:
    # phase ≥ 0 → 0x00000000, phase < 0 → 0xFFFFFFFF. OR with 1 gives 1
    # or -1 in two's complement.
    (reinterpret(Int32, phase) >> 31) | 1
end

function multiply_with_subcarrier!(
    sampled_code::AbstractVector{T},
    modulation::BOCsin,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
    PHASET = Int32,
) where {T}
    fixed_point, subcarrier_phase, delta_subcarrier_phase = calc_subcarrier_phase_and_delta(
        modulation,
        sampling_frequency,
        code_frequency,
        start_phase,
        PHASET,
    )

    @inbounds for (index, i) in enumerate(
        PHASET(start_index):PHASET(length(sampled_code) - 1 + start_index),
    )
        sampled_code[index] *= T(
            calc_subcarrier_bit(
                i,
                fixed_point,
                subcarrier_phase,
                delta_subcarrier_phase,
                PHASET,
            ),
        )
    end
    sampled_code
end

function multiply_with_subcarrier!(
    sampled_code::AbstractVector{Float32},
    modulation::CBOC,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
    PHASET = Int32,
)
    _, subcarrier_phase_boc1, delta_subcarrier_phase_boc1 =
        calc_subcarrier_phase_and_delta(
            modulation.boc1,
            sampling_frequency,
            code_frequency,
            start_phase,
            Int32,
        )

    _, subcarrier_phase_boc2, delta_subcarrier_phase_boc2 =
        calc_subcarrier_phase_and_delta(
            modulation.boc2,
            sampling_frequency,
            code_frequency,
            start_phase,
            Int32,
        )

    boc1_amplitude = Float32(sqrt(modulation.boc1_power))
    boc2_amplitude = Float32(sqrt(1 - modulation.boc1_power))

    N = length(sampled_code)
    # Reinterpret as UInt32 so negative start_index values wrap correctly; the
    # integer phase accumulator is modular.
    si = reinterpret(UInt32, Int32(start_index))
    @inbounds @simd ivdep for index = 1:N
        i = UInt32(index - 1) + si
        b1 = calc_subcarrier_bit_f32(i, subcarrier_phase_boc1, delta_subcarrier_phase_boc1)
        b2 = calc_subcarrier_bit_f32(i, subcarrier_phase_boc2, delta_subcarrier_phase_boc2)
        sampled_code[index] *= boc1_amplitude * b1 + boc2_amplitude * b2
    end
    sampled_code
end

# Fallback for non-Float32 output types: same algorithm, integer subcarrier bits.
function multiply_with_subcarrier!(
    sampled_code::AbstractVector{T},
    modulation::CBOC,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
    PHASET = Int32,
) where {T}
    fixed_point_boc1, subcarrier_phase_boc1, delta_subcarrier_phase_boc1 =
        calc_subcarrier_phase_and_delta(
            modulation.boc1,
            sampling_frequency,
            code_frequency,
            start_phase,
            PHASET,
        )

    fixed_point_boc2, subcarrier_phase_boc2, delta_subcarrier_phase_boc2 =
        calc_subcarrier_phase_and_delta(
            modulation.boc2,
            sampling_frequency,
            code_frequency,
            start_phase,
            PHASET,
        )

    boc1_amplitude = T(sqrt(modulation.boc1_power))
    boc2_amplitude = T(sqrt(1 - modulation.boc1_power))

    @inbounds for (index, i) in enumerate(
        PHASET(start_index):PHASET(length(sampled_code) - 1 + start_index),
    )
        sampled_code[index] *= (
            boc1_amplitude * calc_subcarrier_bit(
                i,
                fixed_point_boc1,
                subcarrier_phase_boc1,
                delta_subcarrier_phase_boc1,
                PHASET,
            ) +
            boc2_amplitude * calc_subcarrier_bit(
                i,
                fixed_point_boc2,
                subcarrier_phase_boc2,
                delta_subcarrier_phase_boc2,
                PHASET,
            )
        )
    end
    sampled_code
end

# Pack a TMBOC pattern tuple into a UInt64 bitmask. NPAT ≤ 64.
# This lets the hot loop do a single `(bits >> chip_pos) & 1` to
# look up the BOC choice rather than indexing a tuple.
@inline function _pack_tmboc_pattern(pattern::NTuple{NPAT, Bool}) where {NPAT}
    @assert NPAT <= 64
    bits = UInt64(0)
    @inbounds for k = 1:NPAT
        if pattern[k]
            bits |= UInt64(1) << (k - 1)
        end
    end
    return bits
end

# TMBOC subcarrier multiply (Int16 buffers only).
#
# Strategy: explicit 16-lane SIMD via SIMD.jl. Per SIMD iteration we
# process 16 samples:
#   1. Compute chip_pos for each lane via a 16-lane SIMD ramp (scalar
#      base, per-lane bump on chip transition).
#   2. `pattern_bits >> chip_pos_lane` (`vpsrlvd`) selects BOC(1,1) vs
#      BOC(6,1) per lane in a single vector instruction — no gather,
#      no per-lane modulo.
#   3. Compute both BOC sign-bits per lane via the sign-bit splice
#      trick (one arithmetic shift + or, vectorized), `vifelse` to
#      select.
#   4. Load 16 buffer samples, multiply, store.
#
# The packed pattern must fit in 32 bits (so `vpsrlvd` is sufficient).
# For L1C-P — TMBOC(6,1,4/33) with BOC(6,1) at chip positions
# {0, 4, 6, 29} — this holds. The Int16 dispatcher below falls back
# to the auto-vectorized two-pass for non-contiguous Int16 buffers and
# for patterns whose packed bits don't fit in `UInt32`. Buffers with
# any other element type (`Float32`, `Float64`, …) raise a
# `MethodError`; convert to `Int16` first.

# Compute the last 1-based sample index in absolute chip count `chip`.
# Returns N if the buffer ends within this chip. Used by the two-pass
# fallback for non-contiguous Int16 buffers.
@inline function _tmboc_chip_end(
    chip::Int,
    chip_acc_fp0::UInt64,
    chip_delta_fp::UInt64,
    N::Int,
)
    threshold = UInt64(chip + 1) << 32
    threshold <= chip_acc_fp0 && return N
    need = threshold - chip_acc_fp0
    samples_into_buffer = div(need + chip_delta_fp - UInt64(1), chip_delta_fp)
    return min(N, Int(samples_into_buffer))
end

# Build a 16-lane Vec{16, UInt32} of chip_pos values for the 16 samples
# starting at `chip_acc_base`. Chip count varies by 0 or 1 across the
# lanes for our typical operating range (~14.6 samples/chip), so the
# common case is a constant broadcast.
@inline function _tmboc_chip_pos_16(
    chip_acc_base::UInt64,
    chip_delta_fp::UInt64,
    int_chip_offset::Int,
    ::Val{NPAT},
) where {NPAT}
    L = 16
    base_count = chip_acc_base >> 32
    chip_step = chip_delta_fp * UInt64(L)
    chip_pos_base = (UInt64(int_chip_offset) + base_count) % UInt64(NPAT)

    # Fast path 1: no chip transition within the 16-lane block. All
    # lanes share `chip_pos_base`. Dominant case whenever the block
    # spans less than one chip — i.e. when the samples-per-chip ratio
    # is larger than the SIMD lane count.
    next_threshold = (base_count + UInt64(1)) << 32
    if chip_acc_base + chip_step <= next_threshold
        return Vec{16, UInt32}(UInt32(chip_pos_base))
    end

    lane_idx = Vec{16, UInt32}((
        UInt32(0), UInt32(1), UInt32(2), UInt32(3),
        UInt32(4), UInt32(5), UInt32(6), UInt32(7),
        UInt32(8), UInt32(9), UInt32(10), UInt32(11),
        UInt32(12), UInt32(13), UInt32(14), UInt32(15),
    ))
    npat_v = Vec{16, UInt32}(UInt32(NPAT))

    # Fast path 2: exactly one chip transition within the block. This
    # is the common case when the samples-per-chip ratio is close to
    # the SIMD lane count. For example at sampling_frequency = 15 MHz
    # and code_frequency = 1.023 MHz the block spans about 1.09 chips,
    # so roughly 91 % of blocks have a single chip transition and 9 %
    # have two. Find the first lane that belongs to the next chip and
    # bump the lanes at or after it by one chip-pos with a single
    # comparison-plus-select.
    second_threshold = (base_count + UInt64(2)) << 32
    if chip_acc_base + chip_step <= second_threshold
        first_lane_in_next_chip = if chip_acc_base >= next_threshold
            0
        else
            Int(cld(next_threshold - chip_acc_base, chip_delta_fp))
        end
        delta_v = vifelse(
            lane_idx < Vec{16, UInt32}(UInt32(first_lane_in_next_chip)),
            Vec{16, UInt32}(0),
            Vec{16, UInt32}(1),
        )
        chip_pos_v = Vec{16, UInt32}(UInt32(chip_pos_base)) + delta_v
        return vifelse(chip_pos_v >= npat_v, chip_pos_v - npat_v, chip_pos_v)
    end

    # General case: two transitions inside the block. The per-lane
    # chip-count delta is 0, 1, or 2 depending on which thresholds the
    # lane's accumulator value has crossed.
    first_lane_in_next_chip = if chip_acc_base >= next_threshold
        0
    else
        Int(cld(next_threshold - chip_acc_base, chip_delta_fp))
    end
    first_lane_two_chips_ahead = if chip_acc_base >= second_threshold
        0
    else
        Int(cld(second_threshold - chip_acc_base, chip_delta_fp))
    end
    one_v = Vec{16, UInt32}(1)
    zero_v = Vec{16, UInt32}(0)
    delta_v =
        vifelse(
            lane_idx < Vec{16, UInt32}(UInt32(first_lane_in_next_chip)),
            zero_v,
            one_v,
        ) +
        vifelse(
            lane_idx < Vec{16, UInt32}(UInt32(first_lane_two_chips_ahead)),
            zero_v,
            one_v,
        )
    chip_pos_v = Vec{16, UInt32}(UInt32(chip_pos_base)) + delta_v
    # Branch-free modulo-NPAT. After adding the per-lane delta,
    # `chip_pos_v` is in `[0, NPAT + 1]`, so at most two successive
    # subtractions of NPAT bring every lane back into the range
    # `[0, NPAT)`.
    chip_pos_v = vifelse(chip_pos_v >= npat_v, chip_pos_v - npat_v, chip_pos_v)
    return vifelse(chip_pos_v >= npat_v, chip_pos_v - npat_v, chip_pos_v)
end

# Two-pass fallback for buffers that don't qualify for the SIMD path.
# Hit when either the buffer isn't a contiguous `Vector{Int16}`
# (SIMD.jl `vload` requires it — e.g. a `view` lands here) or when
# the packed `pattern_bits` doesn't fit in `UInt32`.
#
# Pass 1 multiplies the whole buffer by the BOC(1,1) sign-bit
# (pure auto-vectorized loop). Pass 2 walks chips and, for each
# BOC(6,1) chip, re-multiplies that chip's samples by
# `BOC(1,1) · BOC(6,1)` (which is in {-1, +1}), so the net effect at
# those samples is BOC(6,1) instead of BOC(1,1).
function _tmboc_two_pass_i16!(
    sampled_code::AbstractVector{Int16},
    sc_phase_boc1::UInt32, sc_delta_boc1::UInt32,
    sc_phase_boc2::UInt32, sc_delta_boc2::UInt32,
    pattern_bits::UInt64,
    chip_acc_fp0::UInt64, chip_delta_fp::UInt64,
    int_chip_offset::Int, si::UInt32,
    ::Val{NPAT},
) where {NPAT}
    N = length(sampled_code)
    @inbounds @simd ivdep for index = 1:N
        i = UInt32(index - 1) + si
        b1 = calc_subcarrier_bit_i32(i, sc_phase_boc1, sc_delta_boc1)
        sampled_code[index] = (Int32(sampled_code[index]) * b1) % Int16
    end
    sample_start = 1
    chip = 0
    while sample_start <= N
        sample_end = _tmboc_chip_end(chip, chip_acc_fp0, chip_delta_fp, N)
        chip_pos = (int_chip_offset + chip) % NPAT
        use_boc2 = ((pattern_bits >> chip_pos) & UInt64(1)) != UInt64(0)
        if use_boc2
            @inbounds @simd ivdep for index = sample_start:sample_end
                i = UInt32(index - 1) + si
                b1 = calc_subcarrier_bit_i32(i, sc_phase_boc1, sc_delta_boc1)
                b2 = calc_subcarrier_bit_i32(i, sc_phase_boc2, sc_delta_boc2)
                sampled_code[index] = (Int32(sampled_code[index]) * b1 * b2) % Int16
            end
        end
        sample_start = sample_end + 1
        chip += 1
    end
    sampled_code
end

# 16-lane SIMD.jl Int16 fast path. Requires a contiguous `Vector{Int16}`
# and a pattern that fits in UInt32.
function _tmboc_simd_i16!(
    sampled_code::Vector{Int16},
    sc_phase_boc1::UInt32, sc_delta_boc1::UInt32,
    sc_phase_boc2::UInt32, sc_delta_boc2::UInt32,
    pattern_bits::UInt32,
    chip_acc_fp0::UInt64, chip_delta_fp::UInt64,
    int_chip_offset::Int, si::UInt32,
    ::Val{NPAT},
) where {NPAT}
    N = length(sampled_code)
    L = 16
    Nblk = (N ÷ L) * L

    lane_idx = Vec{16, UInt32}((
        UInt32(0), UInt32(1), UInt32(2), UInt32(3),
        UInt32(4), UInt32(5), UInt32(6), UInt32(7),
        UInt32(8), UInt32(9), UInt32(10), UInt32(11),
        UInt32(12), UInt32(13), UInt32(14), UInt32(15),
    ))
    sc_phase_boc1_v = Vec{16, UInt32}(sc_phase_boc1 + sc_delta_boc1 * si) +
                       lane_idx * Vec{16, UInt32}(sc_delta_boc1)
    sc_phase_boc2_v = Vec{16, UInt32}(sc_phase_boc2 + sc_delta_boc2 * si) +
                       lane_idx * Vec{16, UInt32}(sc_delta_boc2)
    sc_step_boc1 = Vec{16, UInt32}(sc_delta_boc1 * UInt32(L))
    sc_step_boc2 = Vec{16, UInt32}(sc_delta_boc2 * UInt32(L))

    chip_step = chip_delta_fp * UInt64(L)
    chip_acc_base = chip_acc_fp0
    pattern_v = Vec{16, UInt32}(pattern_bits)

    @inbounds for blk = 0:L:(Nblk - 1)
        chip_pos_v = _tmboc_chip_pos_16(chip_acc_base, chip_delta_fp, int_chip_offset, Val(NPAT))
        bits_v = pattern_v >> chip_pos_v
        use_boc2_v = (bits_v & Vec{16, UInt32}(1)) != Vec{16, UInt32}(0)

        phase1 = reinterpret(Vec{16, Int32}, sc_phase_boc1_v)
        phase2 = reinterpret(Vec{16, Int32}, sc_phase_boc2_v)
        b1 = (phase1 >> Int32(31)) | Vec{16, Int32}(1)
        b2 = (phase2 >> Int32(31)) | Vec{16, Int32}(1)
        b = vifelse(use_boc2_v, b2, b1)
        b_i16 = convert(Vec{16, Int16}, b)

        v = vload(Vec{16, Int16}, sampled_code, blk + 1)
        vstore(v * b_i16, sampled_code, blk + 1)

        sc_phase_boc1_v += sc_step_boc1
        sc_phase_boc2_v += sc_step_boc2
        chip_acc_base += chip_step
    end

    # Scalar tail for samples after the last full SIMD block.
    @inbounds for index = (Nblk + 1):N
        i = UInt32(index - 1) + si
        b1 = calc_subcarrier_bit_i32(i, sc_phase_boc1, sc_delta_boc1)
        b2 = calc_subcarrier_bit_i32(i, sc_phase_boc2, sc_delta_boc2)
        acc = chip_acc_fp0 + UInt64(index - 1) * chip_delta_fp
        chip = Int(acc >> 32)
        chip_pos = (int_chip_offset + chip) % NPAT
        bit = (pattern_bits >> chip_pos) & UInt32(1)
        b = bit != 0 ? b2 : b1
        sampled_code[index] = (Int32(sampled_code[index]) * b) % Int16
    end
    sampled_code
end

function multiply_with_subcarrier!(
    sampled_code::AbstractVector{Int16},
    modulation::TMBOC{B1, B2, NPAT},
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
    PHASET = Int32,
) where {B1, B2, NPAT}
    _, sc_phase_boc1, sc_delta_boc1 = calc_subcarrier_phase_and_delta(
        modulation.boc1, sampling_frequency, code_frequency, start_phase, Int32,
    )
    _, sc_phase_boc2, sc_delta_boc2 = calc_subcarrier_phase_and_delta(
        modulation.boc2, sampling_frequency, code_frequency, start_phase, Int32,
    )
    pattern_bits64 = _pack_tmboc_pattern(modulation.pattern)
    chip_delta_fp = ceil(UInt64, Float64(code_frequency / sampling_frequency) * (UInt64(1) << 32))
    int_chip_offset = mod(floor(Int, start_phase), NPAT)
    chip_acc_fp0 = UInt64(floor(UInt64, mod(start_phase, 1.0) * (UInt64(1) << 32)))
    si = reinterpret(UInt32, Int32(start_index))

    if sampled_code isa Vector{Int16} && pattern_bits64 <= UInt64(typemax(UInt32))
        return _tmboc_simd_i16!(
            sampled_code,
            sc_phase_boc1, sc_delta_boc1,
            sc_phase_boc2, sc_delta_boc2,
            UInt32(pattern_bits64),
            chip_acc_fp0, chip_delta_fp,
            int_chip_offset, si, Val(NPAT),
        )
    end
    return _tmboc_two_pass_i16!(
        sampled_code,
        sc_phase_boc1, sc_delta_boc1,
        sc_phase_boc2, sc_delta_boc2,
        pattern_bits64,
        chip_acc_fp0, chip_delta_fp,
        int_chip_offset, si, Val(NPAT),
    )
end

function multiply_with_subcarrier!(
    sampled_code::AbstractVector{T},
    modulation::BOCcos,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
    PHASET = Int32,
) where {T}
    multiply_with_subcarrier!(
        sampled_code,
        BOCsin(modulation.m, modulation.n),
        sampling_frequency,
        code_frequency,
        start_phase + 0.25,
        start_index,
    )
end

function multiply_with_subcarrier!(
    sampled_code::AbstractVector{T},
    modulation::LOC,
    sampling_frequency::Frequency,
    code_frequency::Frequency,
    start_phase = 0.0,
    start_index::Integer = 0,
    PHASET = Int32,
) where {T}
    sampled_code
end

"""
$(SIGNATURES)

Generate a sampled code signal for a given PRN.

Allocates and returns a new buffer containing the spreading code sampled at the
specified sampling frequency. For in-place operation, use [`gen_code!`](@ref).

# Arguments
- `num_samples`: Number of samples to generate
- `signal`: GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B()`)
- `prn`: PRN number of the satellite
- `sampling_frequency`: Sampling frequency (must be larger than code frequency)
- `code_frequency`: Code chipping rate (default: signal's nominal code frequency)
- `start_phase`: Initial code phase in chips (default: 0.0)
- `start_index`: Index offset (default: 0)

# Returns
- `Vector`: Sampled code signal

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
    code = zeros(get_code_type(signal), num_samples)
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
