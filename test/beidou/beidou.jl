# BeiDou Open Service signal tests.
#
# The `ICD_*` verification vectors in `icd_vectors.jl` are the first/last-24-chip
# octal values from the official BDS-SIS-ICD PDFs; every generated primary and
# secondary code is pinned to them below. B1I (whose ICD publishes no chip
# vectors) is checked structurally, and B3I additionally via its published
# register-shift table.
include("icd_vectors.jl")

# Pack a ±1 chip vector into octal, MSB first, matching the ICD chip convention
# (a code chip of +1 → bit 1, -1 → bit 0). Length must be a multiple of 3.
function _bd_octal(v)
    b = map(x -> x > 0 ? 1 : 0, v)
    join(string(4b[i] + 2b[i+1] + b[i+2]) for i = 1:3:length(b))
end

# Assert that column `prn` of `mat` reproduces the ICD (prn, first24, last24)
# octal vectors, for every row in `table`.
function _bd_verify_octal(mat, table)
    for (prn, first24, last24) in table
        col = mat[:, prn]
        @test _bd_octal(col[1:24]) == first24
        @test _bd_octal(col[end-23:end]) == last24
    end
end

# All BeiDou primary codes are ±1 Int8/Int16 with the expected length.
function _bd_check_pm1(codes, len)
    @test size(codes, 1) == len
    @test all(x -> x == 1 || x == -1, codes)
end

include("b1i.jl")
include("b3i.jl")
include("b2b.jl")
include("b2a.jl")
include("b1c.jl")
