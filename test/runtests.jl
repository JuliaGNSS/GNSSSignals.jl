using Test, GNSSSignals, Statistics, Aqua
import Unitful: Hz, MHz, Frequency
import GNSSSignals: BOCsin, BOCcos, CBOC
using CodecZlib: GzipDecompressorStream

# Every concrete signal, band and time system the package defines. A type added
# later must be appended here to be swept by the accessor-layer invariant tests
# ("every signal answers every identity accessor" and friends). Test-suite
# stand-ins (e.g. `TestOnlyBand`) stay out: the invariants are the package's to
# keep.
const ALL_SIGNALS = (
    GPSL1CA,
    GPSL1C_D,
    GPSL1C_P,
    GPSL2CM,
    GPSL2CL,
    GPSL5I,
    GPSL5Q,
    GalileoE1B,
    GalileoE1B_BOC11,
    GalileoE1C,
    GalileoE1C_BOC11,
    GalileoE5aI,
    GalileoE5aQ,
    GalileoE6B,
    GalileoE6C,
    BeiDouB1I,
    BeiDouB3I,
    BeiDouB2bI,
    BeiDouB2aI,
    BeiDouB2aQ,
    BeiDouB1C_D,
    BeiDouB1C_P,
)
const ALL_BANDS = (L1, L2, L5, B1I, B3I, B2b, E6)
const ALL_TIME_SYSTEMS = (GPST, GST, BDT)

# Decode a hex-packed ±1 sample fixture (LSB-first packing into hex
# nibbles, so 4 samples per nibble). Bit 0 of nibble `k` holds sample
# `4k+1`, bit 3 holds sample `4k+4`. Bit value 0 → -1, 1 → +1.
function _load_packed_hex_fixture(filename::AbstractString, n_samples::Integer)
    hex = open(filename) do io
        strip(read(GzipDecompressorStream(io), String))
    end
    out = Vector{Int16}(undef, n_samples)
    @inbounds for k = 1:n_samples
        nibble = parse(UInt8, hex[(k-1)÷4+1]; base = 16)
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
    include("time_systems.jl")
    include("modulation.jl")
    include("gps/l1ca.jl")
    include("gps/l5.jl")
    include("gps/l1c_codes.jl")
    include("gps/l1c_d.jl")
    include("gps/l1c_p.jl")
    include("gps/l2c.jl")
    include("galileo/e1b.jl")
    include("galileo/e1c.jl")
    include("galileo/e5a.jl")
    include("galileo/e6.jl")
    include("beidou/beidou.jl")
    include("common.jl")
end

include("code_lut.jl")
