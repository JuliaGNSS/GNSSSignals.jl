"""
$(SIGNATURES)

Get codes of base GNSS system for generic BOC as a Matrix, where each 
column represents a PRN.
"""
function get_codes(::Type{<:GenericBOC{T}}) where {T<:AbstractGNSSSystem}
    get_codes(T)
end

"""
$(SIGNATURES)

Get code length of base GNSS system for generic BOC.
"""
@inline function get_code_length(::Type{<:GenericBOC{T}}) where T<:AbstractGNSSSystem
    get_code_length(T)
end

"""
$(SIGNATURES)

Get secondary code length of base GNSS system for generic BOC.
"""
@inline function get_secondary_code_length(::Type{<:GenericBOC{T}}) where T<:AbstractGNSSSystem
    get_secondary_code_length(T)
end

"""
$(SIGNATURES)

Get center frequency of base GNSS system for generic BOC.
"""
@inline function get_center_frequency(::Type{<:GenericBOC{T}}) where T<:AbstractGNSSSystem
    get_center_frequency(T)
end

"""
$(SIGNATURES)

Get code frequency of base GNSS system for generic BOC.
"""
function get_code_frequency(::Type{GenericBOC{T,m,n}}) where {T<:AbstractGNSSSystem, m, n}
    n * 1_023_000Hz
end

"""
$(SIGNATURES)

Get data frequency of base GNSS system for generic BOC.
"""
function get_data_frequency(::Type{<:GenericBOC{T}}) where T<:AbstractGNSSSystem
    get_data_frequency(T)
end

"""
$(SIGNATURES)

Get code of type `GenericBOC{T}(m,n) where T<:AbstractGNSSSystem` at 
phase `phase` of PRN `prn`.
"""
Base.@propagate_inbounds function get_code(
    ::Type{GenericBOC{T,m,n}},
    phase,
    prn::Integer
) where {T <: AbstractGNSSSystem, m, n}
    floored_phase = floor(Int, phase)
    floored_boc_phase = floor(Int, phase * 2 * m/n)
    get_code_unsafe(
        T,
        mod(
            floored_phase,
            get_code_length(T) * get_secondary_code_length(T)
        ),
        prn
    ) * (iseven(floored_boc_phase)<<1 - 1)
end

"""
$(SIGNATURES)

Get code of type `GenericBOC{T}(m,n) where T<:AbstractGNSSSystem` at
phase `phase` of PRN `prn`. The phase will not be wrapped by the code 
length. The phase has to be smaller than the code length. 
"""
Base.@propagate_inbounds function get_code_unsafe(
    ::Type{GenericBOC{T,m,n}},
    phase,
    prn::Integer
) where {T <: AbstractGNSSSystem,m,n}
    floored_phase = floor(Int, phase)
    floored_boc_phase = floor(Int, phase*2*m/n)
    get_code_unsafe(T, floored_phase, prn) * 
        (iseven(floored_boc_phase)<<1 - 1)
end