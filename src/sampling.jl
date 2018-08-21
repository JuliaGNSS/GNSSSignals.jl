"""
$(SIGNATURES)

Generate carrier at sample points `samples` with frequency `f`, phase `φ₀` and sampling frequency `f_s`.
# Examples
```julia-repl
julia> gen_carrier(1:4000, 200, 10 * π / 180, 4e6)
```
"""
function gen_carrier(samples, f, φ₀, f_s)
    arg = (2 * π * f / f_s) .* samples .+ φ₀
    sin_sig, cos_sig = Yeppp.sin(arg), Yeppp.cos(arg) # use Yeppp for better performance
    complex.(cos_sig, sin_sig) # or cis.(arg)
end

"""
$(SIGNATURES)

Calculate carrier phase at sample point `sample` with frequency `f`, phase `φ₀` and sampling frequency `f_s`.
# Examples
```julia-repl
julia> get_carrier_phase(4000, 200, 10 * π / 180, 4e6)
```
"""
function get_carrier_phase(sample, f, φ₀, f_s)
    mod2pi((2 * π * f / f_s) * sample + φ₀)
end

"""
$(SIGNATURES)

Generate sampled code at sample points `samples` with the code frequency `f`, code phase `φ₀` and sampling 
frequency `f_s`. The code is provided by `code`.
# Examples
```julia-repl
julia> gen_sat_code(1:4000, 1023e3, 2, 4e6, [1, -1, 1, 1, 1])
```
"""
function gen_sat_code(samples, f, φ₀, f_s, code)
    code_indices = floor.(Int, f ./ f_s .* samples .+ φ₀)
    code_indices .= 1.+ mod.(code_indices, length(code))
        code[code_indices]
end

"""
$(SIGNATURES)

Calculates the code phase at sample point `sample` with the code frequency `f`, code phase `φ₀`, sampling 
frequency `f_s` and code length `code_length`.
# Examples
```julia-repl
julia> get_sat_code_phase(4000, 1023e3, 2, 4e6, 1023)
```
"""
function get_sat_code_phase(sample, f, φ₀, f_s, code_length)
    mod(f / f_s * sample + φ₀ + code_length / 2, code_length) - code_length / 2
end