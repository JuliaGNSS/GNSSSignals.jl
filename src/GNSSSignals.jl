module GNSSSignals

using Dates: DateTime
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
    B2b,
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
    BeiDouB1I,
    BeiDouB3I,
    BeiDouB2bI,
    BeiDouB2aI,
    BeiDouB2aQ,
    BeiDouB1C_D,
    BeiDouB1C_P,
    TMBOC,
    gen_code!,
    gen_code,
    code_engine,
    code_state,
    code_lookup,
    code_advance,
    code_width,
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
    get_band_id,
    get_signal_id,
    get_signal_name,
    TimeSystem,
    GPST,
    GST,
    BDT,
    get_time_system,
    get_system_start_time,
    get_tai_offset,
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

"""
    AbstractGPSSignal{C} <: AbstractGNSSSignal{C}

Abstract supertype for a signal transmitted by the GPS constellation, e.g.
[`GPSL1CA`](@ref), [`GPSL5I`](@ref).

Its purpose is to carry the constellation-level facts that every GPS signal
shares, so they can be stated once instead of per signal. The time system is
the current example: `get_time_system(::Type{<:AbstractGPSSignal}) = GPST()`
covers all GPS signals through subtype dispatch. Genuinely per-signal facts
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

# ─────────────────────────────────────────────────────────────────────────────
# `SignalLUT`: per-SIGNAL embedded LUT holding ALL PRNs' baked code tables as one matrix.
#
# Defined here (before the signal structs that embed it as a `lut::SignalLUT` field, always
# populated) so the struct field types resolve; built and consumed by the LUT adapter in
# `code_lut.jl` (`build_signal_lut`, `gen_code!`). The metadata (subchip_factor P,
# table_length, period_subchips) is identical across the PRNs of a signal — only the baked
# Int8 column differs — so we store one `padded` matrix whose column `prn` is the per-PRN
# baked padded table, plus the shared metadata once. Resampled zero-copy by `gen_code!`
# (a transient view-backed `CodeLUT.CodeTable` over a matrix column). See `build_signal_lut`.
# ─────────────────────────────────────────────────────────────────────────────
struct SignalLUT
    padded::Matrix{Int8}      # (table_length + WINDOW_PAD) × num_prns; column prn = that PRN's padded baked table
    subchip_factor::Int       # P
    # Residual NON-baked secondary applied per primary period, as a (Ls × n) matrix. `Int8` with
    # a single row of 1 when none/baked. Stored per-PRN (column prn) so a long PER-PRN overlay
    # (GPS L1C-P's 1800-chip code) is applied at runtime WITHOUT baking it (which would blow the
    # matrix up ~1800×). A SHARED secondary (GPS L5I's NH10) stores one column reused for all PRNs.
    secondary::Matrix{Int8}
    table_length::Int         # L·P (length of one PRN's baked table, i.e. column length minus WINDOW_PAD)
    period_subchips::Int      # sub-chips per primary period (Lp·P), for secondary application
end

# Residual secondary column for PRN `prn` (the shared-secondary matrix has one column reused
# for every PRN; the per-PRN matrix has one column each).
@inline _signal_lut_secondary(lut::SignalLUT, prn::Int) =
    @view lut.secondary[:, size(lut.secondary, 2) == 1 ? 1 : prn]

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
include("beidou/codes.jl")
include("beidou/b1i.jl")
include("beidou/b3i.jl")
include("beidou/b2b.jl")
include("beidou/b2a.jl")
include("beidou/b1c.jl")
include("time_systems.jl")
include("common.jl")
include("code_lut.jl")
end
