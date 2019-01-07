const CIS_LUT = SVector{64}(cis.((0:63) / 64 * 2π))
function cis_fast(x)
    @inbounds CIS_LUT[(floor(Int, x / 2π * 64) & 63) + 1]
end
"""
$(SIGNATURES)

Generate carrier at sample points `samples` with frequency `f`, phase `φ₀` and sampling frequency `f_s`.
# Examples
```julia-repl
julia> gen_carrier.(1:4000, 200Hz, 10 * π / 180, 4e6Hz)
```
"""
function gen_carrier(sample, f, φ₀, f_s)
    cis(2π * f / f_s * sample + φ₀)
end
function gen_carrier_fast(sample, f, φ₀, f_s)
    cis_fast(2π * f / f_s * sample + φ₀)
end

"""
$(SIGNATURES)

Calculate carrier phase at sample point `sample` with frequency `f`, phase `φ₀` and sampling frequency `f_s`.
# Examples
```julia-repl
julia> get_carrier_phase(4000, 200Hz, 10 * π / 180, 4e6Hz)
```
"""
function calc_carrier_phase(sample, f, φ₀, f_s)
    mod2pi(2π * f / f_s * sample + φ₀)
end

"""
$(SIGNATURES)

Generate sampled code at sample points `samples` with the code frequency `f`, code phase `φ₀` and sampling
frequency `f_s`. The code is provided by `code`.
# Examples
```julia-repl
julia> gen_code.(1:4000, 1023e3Hz, 2, 4e6Hz, Ref([1, -1, 1, 1, 1]))
```
"""
function gen_code(sample, f, φ₀, f_s, codes, prn)
    codes[1 + mod(floor(Int, f / f_s * sample + φ₀), size(codes, 1)), prn]
end

function gen_code(gnss_system::T, sample, f, φ₀, f_s, prn) where T <: AbstractGNSSSystem
    gen_code(sample, f, φ₀, f_s, gnss_system.codes, prn)
end

function gen_code(gnss_system::GPSL1, sample::T, f, φ₀, f_s, prn) where T <: Union{UInt16,UInt32}
    gnss_system.codes[1 + mod_1023(floor(T, f / f_s * sample + φ₀)), prn]
end

function mod_1023(x::UInt16)
    x = (x & 1023) + (x >> 10)
    (x + ((x + 1) >> 10)) & 1023
end

function mod_1023(x::UInt32)
    x = (x & 1023) + (x >> 10) + (x >> 20) + (x >> 30)
    (x + ((x + 1) >> 10)) & 1023
end

"""
$(SIGNATURES)

Calculates the code phase at sample point `sample` with the code frequency `f`, code phase `φ₀`, sampling
frequency `f_s` and code length `code_length`.
# Examples
```julia-repl
julia> calc_code_phase(4000, 1023e3Hz, 2, 4e6Hz, 1023)
```
"""
function calc_code_phase(sample, f, φ₀, f_s, code_length)
    mod(f / f_s * sample + φ₀, code_length)
end

function calc_code_phase_unsafe(sample, f, φ₀, f_s)
    f / f_s * sample + φ₀
end
