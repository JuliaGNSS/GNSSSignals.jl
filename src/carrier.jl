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
    num_samples::Integer = length(phases),
    bits::Val{N} = Val(5)
) where {T <: Integer, N}
    n = N + 2
    fixed_point = 32 - n - 2
    delta = floor(Int32, carrier_frequency * 1 << (fixed_point + n) / sample_frequency)
    fixed_point_start_phase = floor(Int32, start_phase * 1 << (fixed_point + n))
    fixed_point_phase = fixed_point_start_phase - delta
    @inbounds for i = start_sample:num_samples + start_sample - 1
        fixed_point_phase = fixed_point_phase + delta
        phases[i] = T(fixed_point_phase >> fixed_point)
    end
    phases
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
    num_samples::Integer = length(phases),
    bits::Val{N} = Val(5)
) where {VT <: Vector{Int16}, N}
    @avx unroll = 6 for i = start_sample:num_samples + start_sample - 1
        carrier_sin[i] = fpsin(phases[i], bits)
        carrier_cos[i] = fpcos(phases[i], bits)
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
    num_samples::Integer = length(carrier),
    bits::Val{N} = Val(5)
) where N
    fpcarrier_phases!(
        carrier.re,
        carrier_frequency,
        sample_frequency,
        start_phase,
        start_sample = start_sample,
        num_samples = num_samples,
        bits = bits
    )
    fpcarrier!(
        carrier.im,
        carrier.re,
        carrier.re,
        start_sample = start_sample,
        num_samples = num_samples,
        bits = bits
    )
    carrier
end
