"""
$(SIGNATURES)

Get the full code matrix for a GNSS system.

Returns the codes as a matrix where each column represents a PRN.

# Arguments
- `gnss`: A GNSS system instance (e.g., `GPSL1()`, `GPSL5()`, `GalileoE1B()`)

# Returns
- `Matrix`: Code matrix of size `(code_length, num_prns)`

# Examples
```julia-repl
julia> codes = get_codes(GPSL1())
julia> size(codes)
(1023, 37)
```
"""
function get_codes(gnss::AbstractGNSS)
    gnss.codes
end

function calculate_num_inner_iterations(
    gnss::AbstractGNSS,
    maximum_expected_sampling_frequency::Val{MESF},
    maximum_expected_doppler::Val{MED} = Val(8000Hz),
) where {MESF,MED}
    ceil(
        Int,
        MESF / get_code_factor(gnss) /
        (get_code_frequency(gnss) - MED * get_code_center_frequency_ratio(gnss)),
    )
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
- `gnss`: GNSS system instance (e.g., `GPSL1()`, `GPSL5()`, `GalileoE1B()`)
- `prn`: PRN number of the satellite
- `sampling_frequency`: Sampling frequency (must be larger than code frequency)
- `code_frequency`: Code chipping rate (default: system's nominal code frequency)
- `start_phase`: Initial code phase in chips (default: 0.0)
- `start_index_shift`: Index offset for the output buffer (default: 0)
- `maximum_expected_sampling_frequency`: Maximum expected sampling frequency for optimization
- `maximum_expected_doppler`: Maximum expected Doppler frequency (default: 8000 Hz)
- `PHASET`: Integer type for phase calculations (default: `Int32`)

# Returns
- The modified `sampled_code` buffer

# Examples
```julia-repl
julia> using Unitful: MHz
julia> buffer = zeros(Int16, 4000)
julia> gen_code!(buffer, GPSL1(), 1, 4MHz)
```
"""
function gen_code!(
    sampled_code::AbstractVector,
    gnss::AbstractGNSS,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(gnss),
    start_phase = 0.0,
    start_index_shift::Integer = 0,
    maximum_expected_sampling_frequency::Val{MESF} = Val(sampling_frequency),
    maximum_expected_doppler::Val{MED} = Val(8000Hz),
    PHASET = Int32,
) where {MED,MESF}
    sample_code!(
        sampled_code,
        gnss,
        prn,
        sampling_frequency,
        code_frequency,
        start_phase,
        start_index_shift,
        maximum_expected_sampling_frequency,
        maximum_expected_doppler,
    )
    multiply_with_subcarrier!(
        sampled_code,
        get_modulation(gnss),
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
    gnss::AbstractGNSS,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(gnss),
    start_phase = 0.0,
    start_index_shift::Integer = 0,
    maximum_expected_sampling_frequency::Val{MESF} = Val(sampling_frequency),
    maximum_expected_doppler::Val{MED} = Val(8000Hz),
) where {MED,MESF}
    modulated_code_frequency = code_frequency * get_code_factor(gnss)
    frequency_ratio = sampling_frequency / modulated_code_frequency
    modulated_code_frequency > sampling_frequency && error(
        "The sampling frequency must be larger than the code frequency multiplied by code factor (larger than $modulated_code_frequency, it is $sampling_frequency).",
    )

    code_frequency > (get_code_frequency(gnss) + MED) && error(
        "The code frequency $code_frequency is larger than expected ($(get_code_frequency(gnss) + MED))). Please increase the expected maximum Doppler frequency $MED",
    )

    # The -2 (instead of -1) reserves an extra bit for overflow headroom.
    # delta_sum accumulates (num_samples * frequency_ratio) in fixed-point representation.
    # With -1, delta_sum uses nearly the full Int range and can overflow when combined
    # with the (1 << fixed_point) offset in the initial delta_sum calculation.
    fixed_point = sizeof(Int) * 8 - 2 - ndigits(length(sampled_code); base = 2)

    frequency_ratio_fixed_point = round(Int, frequency_ratio * 1 << fixed_point)

    start_phase_including_shift =
        start_phase + start_index_shift * code_frequency / sampling_frequency

    code_length = get_code_length(gnss) * get_secondary_code_length(gnss)
    # Compute fractional part before any normalization to preserve exact floating-point value
    floored_code_phase = floor(Int, start_phase_including_shift)
    frac_part = floored_code_phase - start_phase_including_shift  # Always in (-1, 0]
    # Normalize only the integer part for array indexing
    code_start_index = mod(floored_code_phase, code_length)
    # The -256 offset handles boundary cases where a sample falls exactly on a chip
    # boundary (e.g., phase = 2.0). Without this offset, floor(2.0) = 2 would assign
    # two samples to the first chip when only one belongs there. The offset ensures
    # samples exactly on boundaries are assigned to the next chip.
    delta_sum =
        floor(Int, frac_part * frequency_ratio_fixed_point) + (1 << fixed_point) - 256
    prev = 0
    num_code_samples_to_iterate =
        Int(fld(modulated_code_frequency * length(sampled_code), sampling_frequency))
    num_code_iterations = cld(num_code_samples_to_iterate + code_start_index, code_length)
    num_inner_iterations = calculate_num_inner_iterations(
        gnss,
        maximum_expected_sampling_frequency,
        maximum_expected_doppler,
    )
    processed_code_samples = 0
    @inbounds for k = 0:(num_code_iterations-1)
        iteration_begin = (k == 0 ? code_start_index : 0) + 1
        iteration_end = min(
            num_code_samples_to_iterate + (k == 0 ? code_start_index : 0) -
            processed_code_samples - 1,
            code_length,
        )
        iterations = iteration_begin:iteration_end
        processed_code_samples += length(iterations)
        for i in iterations
            next_code = gnss.codes[i, prn]
            for j = 1:num_inner_iterations
                sampled_code[prev+j] = next_code
            end
            delta_sum += frequency_ratio_fixed_point
            prev = delta_sum >> fixed_point
        end
    end
    @inbounds for i = 0:2
        next_code_idx = mod(processed_code_samples + code_start_index + i, code_length) + 1
        next_code = gnss.codes[next_code_idx, prn]
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

    delta_subcarrier_phase = floor(
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
- `gnss`: GNSS system instance (e.g., `GPSL1()`, `GPSL5()`, `GalileoE1B()`)
- `prn`: PRN number of the satellite
- `sampling_frequency`: Sampling frequency (must be larger than code frequency)
- `code_frequency`: Code chipping rate (default: system's nominal code frequency)
- `start_phase`: Initial code phase in chips (default: 0.0)
- `start_index`: Index offset (default: 0)

# Returns
- `Vector`: Sampled code signal

# Examples
```julia-repl
julia> using Unitful: MHz
julia> sampled_code = gen_code(4000, GPSL1(), 1, 4MHz)
julia> length(sampled_code)
4000
```
"""
function gen_code(
    num_samples::Integer,
    gnss::AbstractGNSS,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(gnss),
    start_phase = 0.0,
    start_index::Integer = 0,
)
    code = zeros(get_code_type(gnss), num_samples)
    gen_code!(code, gnss, prn, sampling_frequency, code_frequency, start_phase, start_index)
end

"""
$(SIGNATURES)

Get the ratio of code frequency to center frequency.

This ratio is used to compute the code Doppler from the carrier Doppler.

# Arguments
- `gnss`: A GNSS system instance

# Returns
- `Float64`: The code-to-center frequency ratio

# Examples
```julia-repl
julia> get_code_center_frequency_ratio(GPSL1())
0.0006493506493506494
```
"""
@inline function get_code_center_frequency_ratio(gnss::AbstractGNSS)
    get_code_frequency(gnss) / get_center_frequency(gnss)
end

"""
$(SIGNATURES)

Get the minimum number of bits needed to represent the code length.

Calculates the number of bits required to represent the full code length,
including secondary code if present.

# Arguments
- `gnss`: A GNSS system instance

# Returns
- `Int`: Number of bits needed

# Examples
```julia-repl
julia> min_bits_for_code_length(GPSL1())
10
julia> min_bits_for_code_length(GPSL5())
17
```
"""
@inline function min_bits_for_code_length(gnss::AbstractGNSS)
    ndigits(get_code_length(gnss) * get_secondary_code_length(gnss); base = 2)
end

"""
$(SIGNATURES)

Get the length of the secondary code.

# Arguments
- `gnss`: A GNSS system instance

# Returns
- `Int`: Secondary code length (1 if no secondary code)

# Examples
```julia-repl
julia> get_secondary_code_length(GPSL1())
1
julia> get_secondary_code_length(GPSL5())
10
```
"""
@inline function get_secondary_code_length(gnss::AbstractGNSS)
    length(get_secondary_code(gnss))
end

"""
$(SIGNATURES)

Get the secondary code value at a given phase.

# Arguments
- `gnss`: A GNSS system instance
- `phase`: Code phase in chips

# Returns
- Secondary code value at the given phase

# Examples
```julia-repl
julia> get_secondary_code(GPSL5(), 10230.0)  # Start of second code period
1
```
"""
@inline function get_secondary_code(gnss::AbstractGNSS, phase)
    get_secondary_code(gnss, get_secondary_code(gnss), phase)
end

"""
$(SIGNATURES)

Get secondary code value when code is a single integer (no secondary code).

Returns the code value unchanged.
"""
@inline function get_secondary_code(gnss::AbstractGNSS, code::Integer, phase)
    code
end

"""
$(SIGNATURES)

Get secondary code value at phase when code is a tuple (has secondary code).

Computes the secondary code index from the phase and returns the corresponding value.
"""
@inline function get_secondary_code(gnss::AbstractGNSS, code::Tuple, phase)
    code[mod(floor(Int, phase / get_code_length(gnss)), get_secondary_code_length(gnss))+1]
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
