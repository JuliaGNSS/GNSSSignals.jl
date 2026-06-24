module GNSSSignals

using DocStringExtensions
using FixedPointNumbers
using SIMD: Vec, vload, vstore, vifelse, shufflevector
using Statistics
using Unitful: Frequency, Hz, upreferred, ustrip, @u_str

import Base.show

export AbstractGNSSSignal,
    Band,
    L1,
    L2,
    L5,
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
    TMBOC,
    gen_code!,
    gen_code,
    CodeReplicaLUT,
    CodeGeneratorLUT,
    get_codes,
    get_code_type,
    get_code_length,
    get_secondary_code_length,
    get_center_frequency,
    get_code_frequency,
    get_code,
    get_code_unsafe,
    get_data_frequency,
    get_code_center_frequency_ratio,
    get_code_spectrum,
    get_band,
    get_signal_name,
    min_bits_for_code_length,
    get_modulation,
    get_secondary_code,
    SecondaryCode,
    NoSecondaryCode,
    SharedSecondaryCode,
    PerPRNSecondaryCode

"""
    AbstractGNSSSignal{C}

Abstract supertype for a GNSS signal.

A *signal* here means a specific transmission such as GPS L1 C/A, GPS L5-I,
or Galileo E1B — the thing with one spreading code, one modulation, one
nominal chip rate, and one RF carrier. This is the level the SIS-ICDs
define: e.g. "GPS L1 C/A Signal Specification".

Concrete subtypes include [`GPSL1CA`](@ref), [`GPSL5I`](@ref), and
[`GalileoE1B`](@ref). The type parameter `C` is the code matrix type.

Signals that share an RF carrier expose the same [`Band`](@ref) via
[`get_band`](@ref); this is what lets a receiver share a carrier NCO
between them (e.g. GPS L1 C/A and Galileo E1B both on 1575.42 MHz).
"""
abstract type AbstractGNSSSignal{C} end

Base.Broadcast.broadcastable(s::AbstractGNSSSignal) = Ref(s)

Base.show(io::IO, x::AbstractGNSSSignal) = print(io, "$(typeof(x))()")

"""
$(SIGNATURES)

Read codes from a binary file.

Reads codes encoded in the specified type from a file. The code length must be provided
by `code_length` and the number of PRNs by `num_prns`.

# Arguments
- `type`: The data type of the codes (e.g., `Int8`)
- `filename`: Path to the binary file containing the codes
- `num_prns`: Number of PRN codes in the file
- `code_length`: Length of each code sequence

# Returns
- `Matrix{type}`: A matrix of size `(code_length, num_prns)` containing the codes

# Examples
```julia-repl
julia> read_in_codes(Int8, "/data/gpsl1cacodes.bin", 32, 1023)
```
"""
function read_in_codes(type, filename, num_prns, code_length)
    open(filename) do file_stream
        read!(file_stream, Array{type}(undef, code_length, num_prns))
    end
end

include("modulation.jl")
include("bands.jl")
include("secondary_codes.jl")
include("gps/l1ca.jl")
include("gps/l5.jl")
include("gps/l1c_constants.jl")
include("gps/l1c_codes.jl")
include("gps/l1c_d.jl")
include("gps/l1c_p.jl")
include("gps/l2c_constants.jl")
include("gps/l2c.jl")
include("galileo/e1b.jl")
include("galileo/e1c.jl")
include("galileo/e5a.jl")
include("common.jl")
include("code_lut.jl")
end
