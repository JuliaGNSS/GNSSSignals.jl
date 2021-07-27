struct BOCcos{T <: AbstractGNSS, M, N} <: AbstractGNSSBOCcos{T, M, N}
    system::T
end

function BOCcos(system::T, m, n) where T <: AbstractGNSS
    BOCcos{T, m, n}(system)
end

"""
$(SIGNATURES)

Get codes of base GNSS system for BOCcos as a Matrix, where each
column represents a PRN.
"""
function get_codes(boc::AbstractGNSSBOCcos)
    get_codes(boc.system)
end

"""
$(SIGNATURES)

Get code length of base GNSS system for BOCcos.
"""
@inline function get_code_length(boc::AbstractGNSSBOCcos)
    get_code_length(boc.system)
end

"""
$(SIGNATURES)

Get secondary code length of base GNSS system for BOCcos.
"""
@inline function get_secondary_code_length(boc::AbstractGNSSBOCcos)
    get_secondary_code_length(boc.system)
end

"""
$(SIGNATURES)

Get center frequency of base GNSS system for BOCcos.
"""
@inline function get_center_frequency(boc::AbstractGNSSBOCcos)
    get_center_frequency(boc.system)
end

"""
$(SIGNATURES)

Get code frequency of base GNSS system for BOCcos.
"""
function get_code_frequency(boc::AbstractGNSSBOCcos{M, N}) where {M, N}
    N * get_code_frequency(boc.system)
end

"""
$(SIGNATURES)

Get data frequency of base GNSS system for BOCcos.
"""
function get_data_frequency(boc::AbstractGNSSBOCcos)
    get_data_frequency(boc.system)
end

"""
$(SIGNATURES)

Get code of BOC at
phase `phase` of PRN `prn`.
"""
Base.@propagate_inbounds function get_code(
    boc::AbstractGNSSBOCcos{M, N},
    phase,
    prn::Integer
) where {M, N}
    floored_phase = floor(Int, phase)
    floored_BOC_phase = floor(Int, phase * 2 * M / N)
    modded_floored_phase = mod(
        floored_phase,
        get_code_length(boc.system) * get_secondary_code_length(boc.system)
    )
    get_code_unsafe(boc.system, modded_floored_phase, prn) *
        (iseven(floored_BOC_phase) << 1 - 1)
end

"""
$(SIGNATURES)

Get code of type BOC at
phase `phase` of PRN `prn`. The phase will not be wrapped by the code
length. The phase has to be smaller than the code length.
"""
Base.@propagate_inbounds function get_code_unsafe(
    boc::AbstractGNSSBOCcos{M, N},
    phase,
    prn::Integer
) where {M, N}
    floored_phase = floor(Int, phase)
    floored_BOC_phase = floor(Int, phase * 2 * M / N)
    get_code_unsafe(boc.system, floored_phase, prn) *
        (iseven(floored_BOC_phase) << 1 - 1)
end
Base.@propagate_inbounds function get_code_unsafe(
    boc::AbstractGNSSBOCcos{M, N},
    phase::Integer,
    prn::Integer
) where {M, N}
    floored_BOC_phase = floor(Integer, phase * 2 * M / N)
    get_code_unsafe(boc.system, phase, prn) *
        (iseven(floored_BOC_phase) << 1 - 1)
end
