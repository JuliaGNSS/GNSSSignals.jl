"""
$(SIGNATURES)

Get codes of GNSS system as a Matrix where each column
represents a PRN.
```julia-repl
julia> get_code(GPSL1())
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

Generate the code signal inplace for PRN-Number `prn` of system `gnss` at chip rate
`code_frequency`, sampled at sampling rate `sampling_frequency`. Make sure, that
`sampling_frequency` is larger than `code_frequency` to avoid overflows with the
modulo calculation.
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

    fixed_point = sizeof(Int) * 8 - 1 - ndigits(length(sampled_code); base = 2)

    frequency_ratio_fixed_point = round(Int, frequency_ratio * 1 << fixed_point)

    start_phase_including_shift =
        start_phase + start_index_shift * code_frequency / sampling_frequency

    code_length = get_code_length(gnss) * get_secondary_code_length(gnss)
    floored_code_phase = floor(Int, start_phase_including_shift)
    code_start_index = mod(floored_code_phase, code_length)
    delta_sum =
        round(
            Int,
            (floored_code_phase - start_phase_including_shift) *
            frequency_ratio_fixed_point,
        ) + 1 << fixed_point
    prev = 0
    num_total_iterations =
        Int(fld(modulated_code_frequency * length(sampled_code), sampling_frequency))
    total_iteration_end = num_total_iterations + code_start_index
    num_code_iterations = cld(total_iteration_end, code_length)
    num_inner_iterations = calculate_num_inner_iterations(
        gnss,
        maximum_expected_sampling_frequency,
        maximum_expected_doppler,
    )
    @inbounds for k = 0:(num_code_iterations-1)
        iteration_begin = (k == 0 ? code_start_index : 0) + 1
        iteration_end = min(total_iteration_end - code_length * k - 1, code_length)
        for i = iteration_begin:iteration_end
            next_code = gnss.codes[i, prn]
            for j = 1:num_inner_iterations
                sampled_code[prev+j] = next_code
            end
            delta_sum += frequency_ratio_fixed_point
            prev = delta_sum >> fixed_point
        end
    end
    @inbounds for i = 0:2
        next_code_idx =
            mod(
                total_iteration_end - code_length * (num_code_iterations - 1) + i - 1,
                code_length,
            ) + 1
        next_code = gnss.codes[next_code_idx, prn]
        num_iterations = min(num_inner_iterations, length(sampled_code) - prev)
        for j = 1:num_iterations
            sampled_code[prev+j] = next_code
        end
        prev + num_iterations == length(sampled_code) && break
        delta_sum += frequency_ratio_fixed_point
        prev = delta_sum >> fixed_point
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

Generate the code signal for PRN-Number `prn` of system `gnss` at chip rate
`code_frequency`, sampled at sampling rate `sampling_frequency`. Make sure, that
`sampling_frequency` is larger than `code_frequency` to avoid overflows with the
modulo calculation.
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

Get code to center frequency ratio
```julia-repl
julia> get_code_center_frequency_ratio(GPSL1())
```
"""
@inline function get_code_center_frequency_ratio(gnss::AbstractGNSS)
    get_code_frequency(gnss) / get_center_frequency(gnss)
end

"""
$(SIGNATURES)

Get the minimum number of bits that are needed to represent the code length
"""
@inline function min_bits_for_code_length(gnss::AbstractGNSS)
    ndigits(get_code_length(gnss) * get_secondary_code_length(gnss); base = 2)
end

"""
$(SIGNATURES)

Get secondary code length
"""
@inline function get_secondary_code_length(gnss::AbstractGNSS)
    length(get_secondary_code(gnss))
end

"""
$(SIGNATURES)

Get secondary code at phase
"""
@inline function get_secondary_code(gnss::AbstractGNSS, phase)
    get_secondary_code(gnss, get_secondary_code(gnss), phase)
end

"""
$(SIGNATURES)

Get secondary code at phase
"""
@inline function get_secondary_code(gnss::AbstractGNSS, code::Integer, phase)
    code
end

"""
$(SIGNATURES)

Get secondary code at phase
"""
@inline function get_secondary_code(gnss::AbstractGNSS, code::Tuple, phase)
    code[mod(floor(Int, phase / get_code_length(gnss)), get_secondary_code_length(gnss))+1]
end

"""
$(SIGNATURES)

Calculate the spectral power of a BPSK modulated signal with chiprate `fc`
at baseband frequency `f`
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
Calculate the spectral power of a sine phased BOC modulated signal with chiprate
`fc` and subcarrier frequency `fs` at baseband frequency `f`
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
Calculate the spectral power of a cosine phased BOC modulated signal with chiprate
`fc` and subcarrier frequency `fs` at baseband frequency `f`
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
