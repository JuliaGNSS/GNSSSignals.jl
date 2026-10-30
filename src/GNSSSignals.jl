module GNSSSignals

using Dates: Dates, DateTime
using DocStringExtensions
using FixedPointNumbers
using SIMD: Vec, vload, vstore, vifelse, shufflevector
using Statistics
using Unitful: Frequency, Hz, s, upreferred, ustrip, @u_str

import Base.show

export AbstractGNSSSignal,
    AbstractGPSSignal,
    AbstractGalileoSignal,
    AbstractBeiDouSignal,
    Band,
    L1,
    L2,
    L5,
    B1I,
    B3I,
    E5b,
    E6,
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
    GalileoE5aQP,
    GalileoE5bI,
    GalileoE5bQ,
    GalileoE6B,
    GalileoE6C,
    BeiDouB1I,
    BeiDouB3I,
    BeiDouB2bI,
    BeiDouB2aI,
    BeiDouB2aQ,
    BeiDouB1C_D,
    BeiDouB1C_P,
    gen_code!,
    gen_code,
    get_codes,
    get_code_length,
    get_secondary_code_length,
    get_center_frequency,
    get_chip_rate,
    get_code,
    get_code_amplitude,
    get_symbol_rate,
    get_code_periods_per_symbol,
    get_carrier_phase_offset,
    get_chip_rate_to_center_frequency_ratio,
    get_code_spectrum,
    get_band,
    get_band_id,
    get_band_name,
    get_signal_id,
    get_constellation_id,
    get_constellation_name,
    get_signal_name,
    get_relative_power,
    TimeSystem,
    GPST,
    GST,
    BDT,
    get_time_system,
    get_time_system_id,
    get_time_system_name,
    get_system_start_time,
    get_tai_system_start_time,
    get_epoch_offset,
    get_time_scale_offset,
    get_week_start_epoch,
    get_week_and_time_of_week,
    get_week,
    get_time_of_week,
    TAITime,
    get_tai_offset,
    get_modulation,
    get_secondary_code

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
between components of one satellite (e.g. GPS L1 C/A and GPS L1C). Signals
from different satellites still have independent carrier phases and Dopplers.
"""
abstract type AbstractGNSSSignal{C} end

"""
    AbstractGPSSignal{C} <: AbstractGNSSSignal{C}

Abstract supertype for a signal transmitted by the GPS constellation, e.g.
[`GPSL1CA`](@ref), [`GPSL5I`](@ref).

Its purpose is to carry the constellation-level facts that every GPS signal
shares, so they can be stated once instead of per signal: the time system
([`get_time_system`](@ref), and with it [`get_time_system_id`](@ref) /
[`get_time_system_name`](@ref)) and the constellation identity
([`get_constellation_id`](@ref) / [`get_constellation_name`](@ref)). Each is one
method on this type — `get_time_system(::Type{<:AbstractGPSSignal}) = GPST()` —
covering all GPS signals through subtype dispatch. Genuinely per-signal facts
([`get_band`](@ref), [`get_modulation`](@ref), the chip rates …) stay defined
on the concrete signal types.
"""
abstract type AbstractGPSSignal{C} <: AbstractGNSSSignal{C} end

"""
    AbstractGalileoSignal{C} <: AbstractGNSSSignal{C}

Abstract supertype for a signal transmitted by the Galileo constellation, e.g.
[`GalileoE1B`](@ref), [`GalileoE5aI`](@ref).

The Galileo counterpart to [`AbstractGPSSignal`](@ref): it carries the facts
every Galileo signal shares, e.g. `get_time_system(::Type{<:AbstractGalileoSignal})
= GST()`.
"""
abstract type AbstractGalileoSignal{C} <: AbstractGNSSSignal{C} end

"""
    AbstractBeiDouSignal{C} <: AbstractGNSSSignal{C}

Abstract supertype for a signal transmitted by the BeiDou constellation, e.g.
[`BeiDouB1I`](@ref), [`BeiDouB2aI`](@ref), [`BeiDouB1C_D`](@ref).

The BeiDou counterpart to [`AbstractGPSSignal`](@ref): it carries the facts
every BeiDou signal shares, e.g. `get_time_system(::Type{<:AbstractBeiDouSignal})
= BDT()` (BeiDou Time).
"""
abstract type AbstractBeiDouSignal{C} <: AbstractGNSSSignal{C} end

Base.Broadcast.broadcastable(s::AbstractGNSSSignal) = Ref(s)

Base.show(io::IO, x::AbstractGNSSSignal) = print(io, "$(typeof(x))()")

include("signal_lut.jl")

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
function read_in_codes(::Type{T}, filename, num_prns, code_length) where {T}
    # The two-argument `open` and an explicit `close`, not the `do`-block form:
    # the closure form's keyword plumbing is the one call Julia's `--trim=safe`
    # verifier cannot resolve, and a statically compiled receiver (the
    # hardware-correlator loop process of GNSSReceiver.jl) constructs its signals
    # through this function.
    file_stream = open(filename, "r")
    try
        read!(file_stream, Array{T}(undef, code_length, num_prns))
    finally
        close(file_stream)
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
include("galileo/codes.jl")
include("galileo/e1b.jl")
include("galileo/e1c.jl")
include("galileo/e5a.jl")
include("galileo/e5a_qp.jl")
include("galileo/e5b.jl")
include("galileo/e6.jl")
include("beidou/codes.jl")
include("beidou/b1i.jl")
include("beidou/b3i.jl")
include("beidou/b2b.jl")
include("beidou/b2a.jl")
include("beidou/b1c.jl")
include("tai_time.jl")
include("time_systems.jl")
include("common.jl")
include("spectrum.jl")
include("code_lut.jl")
end
