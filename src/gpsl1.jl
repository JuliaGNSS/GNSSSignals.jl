
function init_gpsl1_codes()
    code_length = 1023
    codes = read_in_codes(joinpath(Base.Pkg.dir("GNSSSignals"), "data/codes_gps_l1.bin"), code_length)
    gen_code(t, f, φ₀, f_s, sat) = gen_sat_code(t, f, φ₀, f_s, codes[:,sat])
    get_code_phase(t, f, φ₀, f_s) = get_sat_code_phase(t, f, φ₀, f_s, code_length)
    gen_code, get_code_phase
end