#const CIS_LUT = SVector{64}(cis.((0:63) / 64 * 2π))
#function cis_fast(x)
#    @inbounds CIS_LUT[(floor(Int, x / 2π * 64) & 63) + 1]
#end

"""
$(SIGNATURES)

Get sine approximation very fast. It is valid between -π and π.
```julia-repl
julia> sin_vfast(π / 4)
```
"""
@inline function sin_vfast(x::T) where T <: Real
    B = T(4 / π)
    C = T(-4 / (π * π))
    B * x + C * x * abs(x)
end

"""
$(SIGNATURES)

Get cosine approximation very fast. It is valid between -π and π.
```julia-repl
julia> cos_vfast(π / 4)
```
"""
@inline function cos_vfast(x::T) where T <: Real
    x += T(π / 2)
    x -= (x > T(π)) * T(2π)
    sin_vfast(x)
end

"""
$(SIGNATURES)

Get sine approximation fast. It is valid between -π and π.
```julia-repl
julia> sin_fast(π / 4)
```
"""
@inline function sin_fast(x::T) where T <: Real
    P = T(0.225)
    y = sin_vfast(x)
    P * (y * abs(y) - y) + y
end

"""
$(SIGNATURES)

Get cosine approximation fast. It is valid between -π and π.
```julia-repl
julia> cos_fast(π / 4)
```
"""
@inline function cos_fast(x::T) where T <: Real
    x += T(π / 2)
    x -= (x > T(π)) * T(2π)
    sin_fast(x)
end

"""
$(SIGNATURES)

Get \exp(iz) approximation fast. It is valid between -π and π.
```julia-repl
julia> cos_fast(π / 4)
```
"""
function cis_fast(x)
    complex(cos_fast(x), sin_fast(x))
end

"""
$(SIGNATURES)

Get \exp(iz) approximation very fast. It is valid between -π and π.
```julia-repl
julia> cis_vfast(π / 4)
```
"""
function cis_vfast(x)
    complex(cos_vfast(x), sin_vfast(x))
end

"""
$(SIGNATURES)

Get \exp(iz) approximation fast. It is valid between -π and π.
```julia-repl
julia> get_carrier_fast_unsafe(π / 4)
```
"""
@inline function get_carrier_fast_unsafe(x)
    cis_fast(x)
end

"""
$(SIGNATURES)

Get \exp(iz) approximation very fast. It is valid between -π and π.
```julia-repl
julia> get_carrier_vfast_unsafe(π / 4)
```
"""
@inline function get_carrier_vfast_unsafe(x)
    cis_vfast(x)
end
