#const CIS_LUT = SVector{32}(complex.(floor.(Int8, cos.(2π .* (0:31) ./ 32) .* (1 << 4)), floor.(Int8, sin.(2π .* (0:31) ./ 32) .* (1 << 4))))
#cis_fast(x) = @inbounds CIS_LUT[x & 31 + 1]
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


get_quadrant_size_power(x::Int16) = 7
get_carrier_amplitude_power(x::Int16) = 7
@inline function calc_A(x::Int16)
    p = 15; r = 1; A = Int16(23170); C = Int16(-425)
    n = get_quadrant_size_power(x)
    a = get_carrier_amplitude_power(x)
    x² = (x * x) >> n
    rounding = one(x) << (p - a - 1)
    (A + rounding + (x² * C) >> r) >> (p - a)
end
@inline function calc_B(x::Int16)
    p = 14; r = 3; B = Int16(-17790); D = Int16(351)
    n = get_quadrant_size_power(x)
    a = get_carrier_amplitude_power(x)
    x² = (x * x) >> n
    rounding = one(x) << (p - a - 1)
    (rounding + x * (B + (x² * D) >> r) >> n) >> (p - a)
end
@inline function get_first_bit_sign(x)
    n = get_quadrant_size_power(x)
    mysign(x << (sizeof(x) * 8 - n - 1))
end
@inline function get_second_bit_sign(x)
    n = get_quadrant_size_power(x)
    mysign(x << (sizeof(x) * 8 - n - 2))
end
@inline function get_angle(x)
    n = get_quadrant_size_power(x)
    x & (one(x) << n - one(x)) - one(x) << (n - 1)
end
@inline mysign(x) = x >= zero(x) ? one(x) : -one(x)

"""
$(SIGNATURES)

Fixed point sin and cos
"""
@inline function fpsincos(x::Union{Int16, Int32, Int64})
    first_bit_sign = get_first_bit_sign(x)
    second_bit_sign = get_second_bit_sign(x)
    angle = get_angle(x)
    A = calc_A(angle)
    B = calc_B(angle)

    cos_approx = second_bit_sign * (first_bit_sign * A + B)
    sin_approx = second_bit_sign * (A - first_bit_sign * B)
    sin_approx, cos_approx
end
