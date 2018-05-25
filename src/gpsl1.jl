"""
$(SIGNATURES)

Returns functions to generate sampled code and code phase for the GPS L1 signal.
# Examples
```julia-repl
julia> gen_gpsl1_code, get_gpsl1_code_phase = init_gpsl1_codes()
julia> gen_gpsl1_code(samples, f, φ₀, f_s, sat)
julia> get_code_phase(sample, f, φ₀, f_s)
```
"""
function init_gpsl1_codes()
    code_length = 1023
    codes = read_in_codes(joinpath(Base.Pkg.dir("GNSSSignals"), "data/codes_gps_l1.bin"), code_length)
    gen_code(samples, f, φ₀, f_s, sat) = gen_sat_code(samples, f, φ₀, f_s, codes[:,sat])
    get_code_phase(sample, f, φ₀, f_s) = get_sat_code_phase(sample, f, φ₀, f_s, code_length)
    gen_code, get_code_phase
end