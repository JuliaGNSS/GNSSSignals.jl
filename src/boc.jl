"""
$(SIGNATURES)

Get codes of base GNSS system for BOC as a Matrix, where each 
column represents a PRN.
"""
function get_codes(::Type{<:BOC{T}}) where {T <: AbstractGNSSSystem}
    get_codes(T)
end

"""
$(SIGNATURES)

Get code length of base GNSS system for BOC.
"""
@inline function get_code_length(::Type{<:BOC{T}}) where T <: AbstractGNSSSystem
    get_code_length(T)
end

"""
$(SIGNATURES)

Get secondary code length of base GNSS system for BOC.
"""
@inline function get_secondary_code_length(::Type{<:BOC{T}}) where T <: AbstractGNSSSystem
    get_secondary_code_length(T)
end

"""
$(SIGNATURES)

Get center frequency of base GNSS system for BOC.
"""
@inline function get_center_frequency(::Type{<:BOC{T}}) where T <: AbstractGNSSSystem
    get_center_frequency(T)
end

"""
$(SIGNATURES)

Get code frequency of base GNSS system for BOC.
"""
function get_code_frequency(::Type{BOC{T,M,N}}) where {T<:AbstractGNSSSystem, M, N}
    N * 1_023_000Hz
end

"""
$(SIGNATURES)

Get data frequency of base GNSS system for BOC.
"""
function get_data_frequency(::Type{<:BOC{T}}) where T <: AbstractGNSSSystem
    get_data_frequency(T)
end

"""
$(SIGNATURES)

Get code of type `BOC{T}(m,n) where T<:AbstractGNSSSystem` at 
phase `phase` of PRN `prn`.
"""
Base.@propagate_inbounds function get_code(
    ::Type{BOC{T,M,N}},
    phase,
    prn::Integer
) where {T <: AbstractGNSSSystem, M, N}
    floored_phase = floor(Int, phase)
    floored_boc_phase = floor(Int, phase * 2 * M / N)
    get_code_unsafe(
        T,
        mod(
            floored_phase,
            get_code_length(T) * get_secondary_code_length(T)
        ),
        prn
    ) * (iseven(floored_boc_phase) << 1 - 1)
end

"""
$(SIGNATURES)

Get code of type `BOC{T}(m,n) where T <: AbstractGNSSSystem` at
phase `phase` of PRN `prn`. The phase will not be wrapped by the code 
length. The phase has to be smaller than the code length. 
"""
Base.@propagate_inbounds function get_code_unsafe(
    ::Type{BOC{T,M,N}},
    phase,
    prn::Integer
) where {T <: AbstractGNSSSystem, M, N}
    floored_phase = floor(Int, phase)
    floored_boc_phase = floor(Int, phase * 2 * M / N)
    get_code_unsafe(T, floored_phase, prn) * 
        (iseven(floored_boc_phase)<<1 - 1)
end