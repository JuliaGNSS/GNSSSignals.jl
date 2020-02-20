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

get_quadrant_size_power(x::VInt16) = 7
get_carrier_amplitude_power(x::VInt16) = 7
@inline function calc_A(x::VInt16)
    p = 15; r = 1; A = Int16(23170); C = Int16(-425)
    n = get_quadrant_size_power(x)
    a = get_carrier_amplitude_power(x)
    x² = (x * x) >> n
    (A + (x² * C) >> r) >> (p - a)
end
@inline function calc_B(x::VInt16)
    p = 14; r = 3; B = Int16(-17790); D = Int16(351)
    n = get_quadrant_size_power(x)
    a = get_carrier_amplitude_power(x)
    x² = (x * x) >> n
    (x * (B + (x² * D) >> r) >> n) >> (p - a)
end
@inline function get_first_bit_sign(x::VInt16)
    n = get_quadrant_size_power(x)
    mysign(x << (16 - n - 1))
end
@inline function get_second_bit_sign(x::VInt16)
    n = get_quadrant_size_power(x)
    mysign(x << (16 - n - 2))
end
@inline function get_quarter_angle(x)
    n = get_quadrant_size_power(x)
    x & (one(x) << n - one(x)) - one(x) << (n - 1)
end
@inline mysign(x) = vifelse(x >= zero(x), one(x), -one(x))

"""
$(SIGNATURES)

Fixed point cos
"""
@inline function fpcos(phase)
    first_bit_sign = get_first_bit_sign(phase)
    second_bit_sign = get_second_bit_sign(phase)
    quarter_angle = get_quarter_angle(phase)
    A = calc_A(quarter_angle)
    B = calc_B(quarter_angle)

    second_bit_sign * (first_bit_sign * A + B)
end

"""
$(SIGNATURES)

Fixed point sin
"""
@inline function fpsin(phase)
    first_bit_sign = get_first_bit_sign(phase)
    second_bit_sign = get_second_bit_sign(phase)
    quarter_angle = get_quarter_angle(phase)
    A = calc_A(quarter_angle)
    B = calc_B(quarter_angle)

    second_bit_sign * (A - first_bit_sign * B)
end

"""
$(SIGNATURES)

Fixed point sin and cos
"""
@inline function fpsincos(x)
    (fpsin(x), fpcos(x))
end

"""
$(SIGNATURES)

Fixed point carrier phases
"""
function fpcarrier_phases!(
    phases::Vector{T},
    carrier_frequency,
    sample_frequency,
    start_phase::AbstractFloat;
    start_sample::Integer = 1,
    num_samples::Integer = length(phases)
) where T <: Integer
    n = get_quadrant_size_power(zero(T)) + 2
    fixed_point = 32 - n - 2
    delta = floor(Int32, carrier_frequency * 1 << (fixed_point + n) / sample_frequency)
    fixed_point_start_phase = floor(Int32, start_phase * 1 << (fixed_point + n))
    fixed_point_phase = fixed_point_start_phase - delta
    @inbounds for i = start_sample:num_samples + start_sample - 1
        fixed_point_phase = fixed_point_phase + delta
        phases[i] = T(fixed_point_phase >> fixed_point)
    end
end

"""
$(SIGNATURES)

Fixed point carrier
"""
function fpcarrier!(
    carrier_sin::VT,
    carrier_cos::VT,
    phases::VT;
    start_sample::Integer = 1,
    num_samples::Integer = length(phases)
) where VT <: Vector{Int16}
    @avx unroll = 8 for i = start_sample:num_samples + start_sample - 1
        carrier_sin[i] = fpsin(phases[i])
        carrier_cos[i] = fpcos(phases[i])
    end
end

"""
$(SIGNATURES)

Fixed point carrier
"""
function fpcarrier!(
    carrier::StructArray{Complex{Int16}},
    carrier_frequency,
    sample_frequency,
    start_phase::AbstractFloat;
    start_sample::Integer = 1,
    num_samples::Integer = length(carrier)
)
    fpcarrier_phases!(
        carrier.re,
        carrier_frequency,
        sample_frequency,
        start_phase,
        start_sample = start_sample,
        num_samples = num_samples
    )
    fpcarrier!(
        carrier.im,
        carrier.re,
        carrier.re,
        start_sample = start_sample,
        num_samples = num_samples
    )
end
