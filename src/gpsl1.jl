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
function GPSL1()
    code_length = 1023
    codes = read_in_codes(joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_gps_l1.bin"), code_length)
    GPSL1(codes, code_length, 1ms, 1023e3Hz, 1.57542e9Hz, 20)
end
