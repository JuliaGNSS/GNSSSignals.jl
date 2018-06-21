"""
$(SIGNATURES)

Reads codes from a file with filename `filename` (including the path). The code length is provided 
by `code_length`.
# Examples
```julia-repl
julia> read_in_codes("/data/gpsl1codes.bin", 1023)
```
"""
function read_in_codes(filename, code_length)
    file_stats = stat(filename)
    num_prn_codes = floor(Int, file_stats.size / code_length)
    codes = open(filename) do file_stream
        read(file_stream, Int8, code_length, num_prn_codes)
    end
end 


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