"""
    GPSL1C_P{C, M} <: AbstractGPSSignal{C}

GPS L1C pilot signal (the dataless component of L1C, carrying 75% of
the L1C power per IS-GPS-800G; broadcast by Block III/IIIF satellites).

The broadcast signal is TMBOC(6,1,4/33) modulated on the L1 band (1575.42 MHz): every 33
primary chips, four positions `{0, 4, 6, 29}` use BOC(6,1) and the
rest use BOC(1,1). This type provides a continuous BOC(1,1) replica;
use [`GPSL1C_P_TMBOC`](@ref) for the full TMBOC replica.
The 10230-chip primary code at 1.023 Mcps uses the same Weil
construction as L1C-D with different per-PRN parameters. It is modulo-2 added
to an 18 s, 1800-bit per-PRN overlay code (IS-GPS-800G §3.2.2.1.2)
exposed here as a [`PerPRNSecondaryCode`](@ref).

PRNs 1-63 supported.

# Example
```julia
gpsl1c_p = GPSL1C_P()
get_code_length(gpsl1c_p)            # 10230
get_secondary_code_length(gpsl1c_p)  # 1800
get_band(gpsl1c_p)                   # L1()
```
"""
struct GPSL1C_P{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    overlay_codes::M    # 1800 × 63 Int8 ±1 matrix, exposed via PerPRNSecondaryCode
    lut::SignalLUT      # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GPSL1C_P_TMBOC{C, M} <: AbstractGPSSignal{C}

Full TMBOC(6,1,4/33) replica of the GPS L1C pilot. It shares its primary
code, per-PRN overlay, and signal metadata with [`GPSL1C_P`](@ref).
"""
struct GPSL1C_P_TMBOC{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    overlay_codes::M
    lut::SignalLUT
end

const GPSL1CPReplica = Union{GPSL1C_P, GPSL1C_P_TMBOC}

get_modulation(::Type{<:GPSL1C_P}) = BOCsin(1, 1)
get_modulation(::Type{<:GPSL1C_P_TMBOC}) =
    TMBOC(BOCsin(1, 1), BOCsin(6, 1),
          ntuple(k -> (k - 1) ∈ L1C_TMBOC_BOC6_POSITIONS, Val(L1C_TMBOC_PERIOD)))

@inline get_band(::Type{<:GPSL1CPReplica}) = L1()

# 75 % of the L1C composite (IS-GPS-800J, Table 3.2-1). See [`GPSL1C_D`](@ref)'s
# note and [`get_relative_power`](@ref).
@inline get_relative_power(::Type{<:GPSL1CPReplica}) = 0.75

@inline get_signal_name(::Type{<:GPSL1C_P}) = "GPS L1C-P"
@inline get_replica_name(::Type{<:GPSL1C_P}) = "GPS L1C-P (BOC(1,1) approximation)"
@inline get_signal_name(::Type{<:GPSL1C_P_TMBOC}) = "GPS L1C-P"
@inline get_replica_name(::Type{<:GPSL1C_P_TMBOC}) = "GPS L1C-P (TMBOC)"
@inline get_signal_id(::Type{<:GPSL1C_P_TMBOC}) = :GPSL1C_P

function read_gpsl1c_p_codes()
    _l1c_build_primary_codes(L1C_P_WEIL_INDEX, L1C_P_INSERTION_INDEX)
end

function _gpsl1c_p_replica(::Type{T}) where {T<:GPSL1CPReplica}
    codes = read_gpsl1c_p_codes()
    overlay = _l1c_build_overlay_codes()
    # The 1800-chip per-PRN overlay is far too long to bake (1800·10230·12 ≫ typemax(Int16)),
    # so it stays residual in the SignalLUT and is applied per primary period at gen time.
    lut = build_signal_lut(get_modulation(T), codes, PerPRNSecondaryCode(overlay))
    T(codes, overlay, lut)
end

GPSL1C_P() = _gpsl1c_p_replica(GPSL1C_P)
GPSL1C_P_TMBOC() = _gpsl1c_p_replica(GPSL1C_P_TMBOC)

"""
$(SIGNATURES)

Get the primary code length for GPS L1C-P.

# Returns
- `Int`: 10230 chips
"""
@inline get_code_length(::Type{<:GPSL1CPReplica}) = L1C_PRIMARY_LENGTH

"""
$(SIGNATURES)

Get the secondary code (overlay) for GPS L1C-P.

Per IS-GPS-800G §3.2.2.1.2, every primary period (10 ms) is XOR'd with
one chip of a 1800-bit per-PRN LFSR-generated overlay code, giving an
18 s total cycle.

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 1800 × 63 overlay matrix
"""
@inline get_secondary_code(s::GPSL1CPReplica) = PerPRNSecondaryCode(s.overlay_codes)

"""
$(SIGNATURES)

Get the nominal chip rate of GPS L1C-P.

# Returns
- `Frequency`: 1.023 MHz
"""
@inline get_chip_rate(::Type{<:GPSL1CPReplica}) = 1_023_000Hz

"""
$(SIGNATURES)

Get the bit rate, as tracked, of GPS L1C-P.

The pilot component is dataless, so its bit rate is 0 Hz.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_bit_rate(::Type{<:GPSL1CPReplica}) = 0Hz

@inline get_supported_prns(::Type{<:GPSL1CPReplica}) = 1:63
@inline get_secondary_code_length(::Type{<:GPSL1CPReplica}) = 1800
