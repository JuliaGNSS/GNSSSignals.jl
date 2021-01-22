"""
$(SIGNATURES)

Get code of type <: `AbstractGNSSSystem` at phase `phase` of PRN `prn`.
```julia-repl
julia> get_code(GPSL1, 1200.3, 1)
```
"""
Base.@propagate_inbounds function get_code(
    ::Type{T},
    phase,
    prn::Integer
) where T <: AbstractGNSSSystem
    floored_phase = floor(Int, phase)
    get_code_unsafe(
        T,
        mod(floored_phase, get_code_length(T) * get_secondary_code_length(T)),
        prn
    )
end

"""
$(SIGNATURES)

Get code of type <: `AbstractGNSSSystem` at phase `phase` of PRN `prn`.
The phase will not be wrapped by the code length. The phase has to smaller
than the code length incl. secondary code.
```julia-repl
julia> get_code_unsafe(GPSL1, 10.3, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(
    ::Type{T},
    phase,
    prn::Integer
) where T <: AbstractGNSSSystem
    get_code_unsafe(T, floor(Int, phase), prn)
end

"""
$(SIGNATURES)

Get code to center frequency ratio
```julia-repl
julia> get_code_unsafe(GPSL1, 10.3, 1)
```
"""
@inline function get_code_center_frequency_ratio(::Type{T}) where T <: AbstractGNSSSystem
    get_code_frequency(T) / get_center_frequency(T)
end

"""
$(SIGNATURES)

Minimum bits that are needed to represent the code length
"""
function min_bits_for_code_length(::Type{S}) where S <: AbstractGNSSSystem
    for i = 1:32
        if get_code_length(S) * get_secondary_code_length(S) <= 1 << i
            return i
        end
    end
    return 0
end

"""
$(SIGNATURES)

Calculate the baseband spectrum of an BPSK modulated signal
"""
function get_code_spectrum_BPSK(code_frequency::Frequency, frequencies)
    return get_code_spectrum_BPSK(code_frequency/1Hz, frequencies)
end
function get_code_spectrum_BPSK(code_frequency, frequencies::Frequency)
    return get_code_spectrum_BPSK(code_frequency, frequencies/1Hz)
end
function get_code_spectrum_BPSK(code_frequency::Frequency, frequencies::Frequency)
    return get_code_spectrum_BPSK(code_frequency/1Hz, frequencies/1Hz)
end
function get_code_spectrum_BPSK(code_frequency, frequencies)
    return sinc.(frequencies./code_frequency).^2 ./ code_frequency
end


"""
$(SIGNATURES)

Calculate the baseband spectrum of a sine phased BOC modulated signal
"""
function get_code_spectrum_BOCsin(
    code_frequency::Frequency, 
    subcarrier_frequency::Frequency, 
    frequencies
)
    return get_code_spectrum_BOCsin(code_frequency/1Hz, subcarrier_frequency/1Hz, frequencies)
end
function get_code_spectrum_BOCsin(
    code_frequency, 
    subcarrier_frequency, 
    frequencies::Frequency
)
    return get_code_spectrum_BOCsin(code_frequency, subcarrier_frequency, frequencies/1Hz)
end
function get_code_spectrum_BOCsin(
    code_frequency::Frequency, 
    subcarrier_frequency::Frequency, 
    frequencies::Frequency
)
    return get_code_spectrum_BOCsin(code_frequency/1Hz, subcarrier_frequency/1Hz, frequencies/1Hz)
end
function get_code_spectrum_BOCsin(code_frequency, subcarrier_frequency, frequencies)
    return ((sinc.(frequencies./code_frequency) 
        .* tan.(pi.*frequencies./(2*subcarrier_frequency))).^2 
        ./ code_frequency)
end

"""
$(SIGNATURES)

Calculate the baseband spectrum of a cosine phased BOC modulated signal
"""
function get_code_spectrum_BOCcos(
    code_frequency::Frequency, 
    subcarrier_frequency::Frequency, 
    frequencies
)
    return get_code_spectrum_BOCcos(code_frequency/1Hz, subcarrier_frequency/1Hz, frequencies)
end
function get_code_spectrum_BOCcos(
    code_frequency, 
    subcarrier_frequency, 
    frequencies::Frequency
)
    return get_code_spectrum_BOCcos(code_frequency, subcarrier_frequency, frequencies/1Hz)
end
function get_code_spectrum_BOCcos(
    code_frequency::Frequency, 
    subcarrier_frequency::Frequency, 
    frequencies::Frequency
)
    return get_code_spectrum_BOCcos(code_frequency/1Hz, subcarrier_frequency/1Hz, frequencies/1Hz)
end
function get_code_spectrum_BOCcos(code_frequency, subcarrier_frequency, frequencies)
    fs = subcarrier_frequency
    fc = code_frequency
    f  = frequencies
    
    @. return (2 * sinc(f/fc) * sinpi(f/4fs)^2 / cospi(f/2fs))^2 / fc
end