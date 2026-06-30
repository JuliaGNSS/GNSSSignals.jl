using Test, GNSSSignals, Statistics, Aqua
import Unitful: Hz, MHz, Frequency
import GNSSSignals: BOCsin, BOCcos, CBOC
using CodecZlib: GzipDecompressorStream

# Decode a hex-packed ±1 sample fixture (LSB-first packing into hex
# nibbles, so 4 samples per nibble). Bit 0 of nibble `k` holds sample
# `4k+1`, bit 3 holds sample `4k+4`. Bit value 0 → -1, 1 → +1.
function _load_packed_hex_fixture(filename::AbstractString, n_samples::Integer)
    hex = open(filename) do io
        strip(read(GzipDecompressorStream(io), String))
    end
    out = Vector{Int16}(undef, n_samples)
    @inbounds for k = 1:n_samples
        nibble = parse(UInt8, hex[(k - 1) ÷ 4 + 1]; base = 16)
        bit = (nibble >> ((k - 1) % 4)) & UInt8(1)
        out[k] = bit == 0 ? Int16(-1) : Int16(1)
    end
    out
end

# The AVX-512 SDE-emulation CI job sets `GNSS_TEST_SIMD_ONLY` and exists only to exercise the
# SIMD backends (`code_lut.jl` forces every host backend). The rest of the suite is CPU-agnostic
# and fully covered by the native CI matrix, so scope the emulated run to the SIMD tests: under
# emulation the full suite is the slowest job *and* flakes Aqua's timeout-based
# `persistent_tasks` check. Native runs (no env var) run everything.
const SIMD_ONLY = haskey(ENV, "GNSS_TEST_SIMD_ONLY")

if !SIMD_ONLY
    @testset "Aqua.jl" begin
        Aqua.test_all(GNSSSignals; deps_compat=(check_extras=false,))
    end

    include("test_codes.jl")
    include("bands.jl")
    include("modulation.jl")
    include("gps/l1ca.jl")
    include("gps/l5i.jl")
    include("gps/l1c_codes.jl")
    include("gps/l1c_d.jl")
    include("gps/l1c_p.jl")
    include("galileo/e1b.jl")
    include("common.jl")
end

include("code_lut.jl")
