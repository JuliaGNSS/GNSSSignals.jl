"""
$(SIGNATURES)

Get codes of GNSS system as a Matrix where each column
represents a PRN.
```julia-repl
julia> get_code(gpsl1)
```
"""
function get_codes(gnss::AbstractGNSS{T}) where T <: AbstractMatrix
    @view gnss.codes[get_code_length(gnss) + 1:2 * get_code_length(gnss), :]
end
function get_codes(gnss::AbstractGNSS{T}) where T <: CuMatrix
    gnss.codes
end

"""
$(SIGNATURES)

Get code of type <: `AbstractGNSS` at phase `phase` of PRN `prn`.
```julia-repl
julia> get_code(GPSL1, 1200.3, 1)
```
"""
Base.@propagate_inbounds function get_code(
    gnss::AbstractGNSS,
    phase,
    prn::Integer
)
    floored_phase = floor(Int, phase)
    get_code_unsafe(
        gnss,
        mod(floored_phase, get_code_length(gnss) * get_secondary_code_length(gnss)),
        prn
    )
end

"""
$(SIGNATURES)

Get code of type <: `AbstractGNSS` at phase `phase` of PRN `prn`.
The phase will not be wrapped by the code length. The phase has to smaller
than the code length incl. secondary code.
```julia-repl
julia> get_code_unsafe(GPSL1, 10.3, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(
    gnss::AbstractGNSS,
    phase,
    prn::Integer
)
    get_code_unsafe(gnss, floor(Int, phase), prn)
end

"""
$(SIGNATURES)

Get code of GNSS system at phase `phase` of prn `prn`.
The phase will not be wrapped by the code length. The phase has to smaller
than the code length and must be an integer.
```julia-repl
julia> get_code_unsafe(gpsl1, 10, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(
    gnss::AbstractGNSS,
    phase::Integer,
    prn::Integer
)
    gnss.codes[get_code_length(gnss) + phase + 1, prn]
end
Base.@propagate_inbounds function get_code_unsafe(
    gnss::AbstractGNSS{C},
    phase::Integer,
    prn::Integer
) where {C <: CuMatrix}
    gnss.codes[phase + 1, prn]
end


"""
$(SIGNATURES)

Get code to center frequency ratio
```julia-repl
julia> get_code_unsafe(GPSL1, 10.3, 1)
```
"""
@inline function get_code_center_frequency_ratio(gnss::AbstractGNSS)
    get_code_frequency(gnss) / get_center_frequency(gnss)
end

"""
$(SIGNATURES)

Minimum bits that are needed to represent the code length
"""
@inline function min_bits_for_code_length(gnss::AbstractGNSS)
    ndigits(get_code_length(gnss) * get_secondary_code_length(gnss); base=2)
end


"""
$(SIGNATURES)
Calculate the spectral power of a BPSK modulated signal with chiprate `fc`
at baseband frequency `f`
"""
function get_code_spectrum_BPSK(fc::Frequency, f)
    return get_code_spectrum_BPSK(fc/1Hz, f)
end
function get_code_spectrum_BPSK(fc, f::Frequency)
    return get_code_spectrum_BPSK(fc, f/1Hz)
end
function get_code_spectrum_BPSK(fc::Frequency, f::Frequency)
    return get_code_spectrum_BPSK(fc/1Hz, f/1Hz)
end
function get_code_spectrum_BPSK(fc, f)
    return sinc(f/fc)^2 / fc
end


"""
$(SIGNATURES)
Calculate the spectral power of a sine phased BOC modulated signal with chiprate
`fc` and subcarrier frequency `fs` at baseband frequency `f`
"""
function get_code_spectrum_BOCsin(fc::Frequency, fs::Frequency, f)
    return get_code_spectrum_BOCsin(fc/1Hz, fs/1Hz, f)
end
function get_code_spectrum_BOCsin(fc, fs, f::Frequency)
    return get_code_spectrum_BOCsin(fc, fs, f/1Hz)
end
function get_code_spectrum_BOCsin(fc::Frequency, fs::Frequency, f::Frequency)
    return get_code_spectrum_BOCsin(fc/1Hz, fs/1Hz, f/1Hz)
end
function get_code_spectrum_BOCsin(fc, fs, f)
    return ((sinc(f/fc) * tan(pi*f/(2*fs)))^2/ fc)
end


"""
$(SIGNATURES)
Calculate the spectral power of a cosine phased BOC modulated signal with chiprate
`fc` and subcarrier frequency `fs` at baseband frequency `f`
"""
function get_code_spectrum_BOCcos(fc::Frequency, fs::Frequency, f)
    return get_code_spectrum_BOCcos(fc/1Hz, fs/1Hz, f)
end
function get_code_spectrum_BOCcos(fc, fs, f::Frequency)
    return get_code_spectrum_BOCcos(fc, fs, f/1Hz)
end
function get_code_spectrum_BOCcos(fc::Frequency, fs::Frequency, f::Frequency)
    return get_code_spectrum_BOCcos(fc/1Hz, fs/1Hz, f/1Hz)
end
function get_code_spectrum_BOCcos(fc, fs, f)
    return (2 * sinc(f/fc) * sinpi(f/4fs)^2 / cospi(f/2fs))^2 / fc
end
