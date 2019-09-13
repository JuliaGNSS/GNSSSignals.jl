"""
$(SIGNATURES)

Get code of type <: `AbstractGNSSSystem` at phase `phase` of prn `prn`.
```julia-repl
julia> get_code(GPSL1, 1200.3, 1)
```
"""
function get_code(::Type{T}, phase, prn::Int) where T <: AbstractGNSSSystem
    get_code_unsafe(T, mod(phase, get_code_length(T)), prn)
end

"""
$(SIGNATURES)

Get code of type <: `AbstractGNSSSystem` at phase `phase` of prn `prn`.
The phase will not be wrapped by the code length. The phase has to smaller
than the code length.
```julia-repl
julia> get_code_unsafe(GPSL1, 10.3, 1)
```
"""
function get_code_unsafe(::Type{T}, phase, prn::Int) where T <: AbstractGNSSSystem
    get_code_unsafe(T, floor(Int, phase), prn::Int)
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
