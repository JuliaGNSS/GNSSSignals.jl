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

"""
$(SIGNATURES)

Generate the code signal inplace for PRN-Number `prn` of system `gnss` at chip rate
`code_frequency`, sampled at sampling rate `sampling_frequency`. Make sure, that
`sampling_frequency` is larger than `code_frequency` to avoid overflows with the
modulo calculation.
"""
function gen_code!(
    code::AbstractVector,
    gnss::AbstractGNSS,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(gnss),
    start_phase = 0.0,
    start_index::Integer = 0,
)
    code_frequency > sampling_frequency &&
        error("The code freqeuncy must not be larger than the sampling frequency.")
    num_samples = length(code)
    fixed_point = sizeof(Int) * 8 - 1 - min_bits_for_code_length(gnss)
    FP = Fixed{Int,fixed_point}
    total_code_length = get_code_length(gnss) * get_secondary_code_length(gnss)
    fp_total_code_length = FP(total_code_length / get_code_factor(gnss))
    delta = FP(upreferred(code_frequency / sampling_frequency))
    code_phase = FP(
        mod(
            FP(mod(start_phase, total_code_length)) + start_index * delta,
            total_code_length,
        ),
    )
    @inbounds for i ∈ 1:num_samples
        code[i] = get_code_unsafe(gnss, code_phase, prn)
        code_phase += delta
        code_phase -= (code_phase >= fp_total_code_length) * fp_total_code_length
    end
    return code
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
