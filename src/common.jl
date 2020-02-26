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
