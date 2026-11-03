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
    B3,
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
    GalileoE1C,
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
    BOCsin,
    CBOC,
    TMBOC,
    gen_code!,
    gen_code,
    get_codes,
    get_supported_prns,
    get_code_length,
    get_secondary_code_length,
    get_center_frequency,
    get_chip_rate,
    get_code,
    get_code_amplitude,
    get_bit_rate,
    get_code_periods_per_bit,
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
    get_time_system_offset,
    get_week_start_epoch,
    get_week,
    get_time_of_week,
    fold_week_crossover,
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

Carries the constellation-level facts every GPS signal shares: the time system
([`get_time_system`](@ref), with [`get_time_system_id`](@ref) /
[`get_time_system_name`](@ref)) and the constellation identity
([`get_constellation_id`](@ref) / [`get_constellation_name`](@ref)), each one method
on this type, e.g. `get_time_system(::Type{<:AbstractGPSSignal}) = GPST()`. Per-signal
facts ([`get_band`](@ref), [`get_modulation`](@ref), the chip rates …) stay on the
concrete signal types.
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
# `SignalLUT`: per-signal embedded LUT holding all PRNs' baked code tables as one matrix.
#
# Defined before the signal structs that embed it (`lut::SignalLUT`); built and consumed by
# `code_lut.jl` (`build_signal_lut`, `gen_code!`). The metadata is identical across a
# signal's PRNs, so one `padded` matrix (column `prn` = that PRN's baked table) plus the
# shared metadata is stored. `gen_code!` resamples it zero-copy.
# ─────────────────────────────────────────────────────────────────────────────
struct SignalLUT
    padded::Matrix{Int8}      # (table_length + WINDOW_PAD) × num_prns; column prn = that PRN's padded baked table
    subchip_factor::Int       # P
    # Residual non-baked secondary applied per primary period, as a (Ls × n) matrix; a single
    # row of 1 when none/baked. A long per-PRN overlay (GPS L1C-P's 1800 chips) is stored per
    # column and applied at runtime, since baking it would grow the matrix ~1800×. A shared
    # secondary (GPS L5I's NH10) stores one column reused for all PRNs.
    secondary::Matrix{Int8}
    table_length::Int         # L·P (length of one PRN's baked table, i.e. column length minus WINDOW_PAD)
    period_subchips::Int      # sub-chips per primary period (Lp·P), for secondary application
    # RMS amplitude of the baked code, `sqrt(mean(code^2))` over one primary table. Exactly
    # `1` for ±1 codes (BPSK/BOC/TMBOC); for the multi-level CBOC Int8 approximation it is
    # `sqrt(a1^2 + a2^2)` (≈ 19.92 for E1B's baked (19, 6)) and tracks any `cboc_amplitudes`
    # override. PRN-independent. Read via `get_code_amplitude`; divide correlation results
    # by it to get a unit-power scale.
    code_amplitude::Float64
end

# Residual secondary column for PRN `prn` (a shared secondary has one column for all PRNs).
@inline _signal_lut_secondary(lut::SignalLUT, prn::Int) =
    @view lut.secondary[:, size(lut.secondary, 2) == 1 ? 1 : prn]

"""
$(SIGNATURES)

Read codes of the given type from a binary file, returning a `Matrix{type}` of size
`(code_length, num_prns)`.

# Examples
```julia-repl
julia> read_in_codes(Int8, "/data/gpsl1cacodes.bin", 32, 1023)
```
"""
function read_in_codes(::Type{T}, filename, num_prns, code_length) where {T}
    # Two-argument `open` and explicit `close`, not the `do`-block form: its keyword
    # plumbing is the one call Julia's `--trim=safe` verifier cannot resolve.
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
