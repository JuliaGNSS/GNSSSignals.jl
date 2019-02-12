"""
$(SIGNATURES)

Returns a `GPSL1 <: AbstractGNSSSystem` type which holds information about the
GPSL5. It can e.g. be used to generate the PRN code.
"""
function GPSL1()
    code_length = 1023
    codes = read_in_codes(joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_gps_l1.bin"), code_length)
    GPSL1(codes, code_length, 1023e3Hz, 1.57542e9Hz, 20)
end
