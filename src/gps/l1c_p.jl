"""
    GPSL1C_P{C, M, R} <: AbstractGPSSignal{C}
    GPSL1C_P()
    GPSL1C_P(BOCsin)
    GPSL1C_P(TMBOC)

GPS L1C pilot signal (the dataless component of L1C, carrying 75% of
the L1C power per IS-GPS-800G; broadcast by Block III/IIIF satellites).

Broadcast as TMBOC(6,1,4/33) on the L1 band (1575.42 MHz): every 33
primary chips, four positions `{0, 4, 6, 29}` use BOC(6,1) and the
rest use BOC(1,1). 10230-chip primary code at 1.023 Mcps (same Weil
construction as L1C-D, different per-PRN parameters), modulo-2 added
with an 18 s, 1800-bit per-PRN overlay code (IS-GPS-800G §3.2.2.1.2)
exposed here as a [`PerPRNSecondaryCode`](@ref).

The replica's modulation is chosen at construction and carried in the type parameter `R`:

- `GPSL1C_P()`, the same as `GPSL1C_P(BOCsin)`, replicates a continuous BOC(1,1). The
  four BOC(6,1) chips in 33 correlate to almost nothing against it, so for an ideal
  wideband received signal the correlation loss is ≈ 1.12 dB (`20 log10(29/33)`); the
  loss depends on receiver filtering. Its code table uses two samples per primary chip
  instead of twelve.
- `GPSL1C_P(TMBOC)` replicates the full TMBOC(6,1,4/33) signal.

Both are replicas of one transmitted component: they share the primary and overlay codes,
rates, band, relative power, signal id (`:GPSL1C_P`) and name. [`get_modulation`](@ref)
tells them apart.

PRNs 1-63 supported.

# Example
```julia
gpsl1c_p = GPSL1C_P()
get_code_length(gpsl1c_p)            # 10230
get_secondary_code_length(gpsl1c_p)  # 1800
get_band(gpsl1c_p)                   # L1()
get_modulation(gpsl1c_p)             # BOCsin(1, 1)
get_modulation(GPSL1C_P(TMBOC))      # TMBOC(BOCsin(1, 1), BOCsin(6, 1), …)
```
"""
struct GPSL1C_P{C<:AbstractMatrix,M<:AbstractMatrix,R<:Modulation} <: AbstractGPSSignal{C}
    codes::C
    overlay_codes::M    # 1800 × 63 Int8 ±1 matrix, exposed via PerPRNSecondaryCode
    lut::SignalLUT      # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GPSL1C_P{<:AbstractMatrix,<:AbstractMatrix,BOCsin}}) = BOCsin(1, 1)
# The bare type stands for the default replica, as `GPSL1C_P()` does.
get_modulation(::Type{GPSL1C_P}) = BOCsin(1, 1)
get_modulation(::Type{<:GPSL1C_P{<:AbstractMatrix,<:AbstractMatrix,TMBOC}}) =
    TMBOC(BOCsin(1, 1), BOCsin(6, 1),
          ntuple(k -> (k - 1) ∈ L1C_TMBOC_BOC6_POSITIONS, Val(L1C_TMBOC_PERIOD)))

@inline get_band(::Type{<:GPSL1C_P}) = L1()

# 75 % of the L1C composite (IS-GPS-800J, Table 3.2-1), for either replica. See
# [`GPSL1C_D`](@ref)'s note and [`get_relative_power`](@ref).
@inline get_relative_power(::Type{<:GPSL1C_P}) = 0.75

@inline get_signal_name(::Type{<:GPSL1C_P}) = "GPS L1C-P"

function read_gpsl1c_p_codes()
    _l1c_build_primary_codes(L1C_P_WEIL_INDEX, L1C_P_INSERTION_INDEX)
end

GPSL1C_P() = GPSL1C_P(BOCsin)
GPSL1C_P(::Type{BOCsin}) = _gpsl1c_p(BOCsin)
GPSL1C_P(::Type{TMBOC}) = _gpsl1c_p(TMBOC)

function _gpsl1c_p(::Type{R}) where {R<:Modulation}
    codes = read_gpsl1c_p_codes()
    overlay = _l1c_build_overlay_codes()
    S = GPSL1C_P{typeof(codes),typeof(overlay),R}
    # The 1800-chip per-PRN overlay is far too long to bake (1800·10230·12 ≫ typemax(Int16)),
    # so it stays residual in the SignalLUT and is applied per primary period at gen time.
    S(codes, overlay, build_signal_lut(get_modulation(S), codes, PerPRNSecondaryCode(overlay)))
end

"""
$(SIGNATURES)

Get the primary code length for GPS L1C-P.

# Returns
- `Int`: 10230 chips
"""
@inline get_code_length(::Type{<:GPSL1C_P}) = L1C_PRIMARY_LENGTH

"""
$(SIGNATURES)

Get the secondary code (overlay) for GPS L1C-P.

Per IS-GPS-800G §3.2.2.1.2, every primary period (10 ms) is XOR'd with
one chip of a 1800-bit per-PRN LFSR-generated overlay code, giving an
18 s total cycle.

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 1800 × 63 overlay matrix
"""
@inline get_secondary_code(s::GPSL1C_P) = PerPRNSecondaryCode(s.overlay_codes)

"""
$(SIGNATURES)

Get the nominal chip rate of GPS L1C-P.

# Returns
- `Frequency`: 1.023 MHz
"""
@inline get_chip_rate(::Type{<:GPSL1C_P}) = 1_023_000Hz

"""
$(SIGNATURES)

Get the bit rate, as tracked, of GPS L1C-P.

The pilot component is dataless, so its bit rate is 0 Hz.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_bit_rate(::Type{<:GPSL1C_P}) = 0Hz

@inline get_supported_prns(::Type{<:GPSL1C_P}) = 1:63
@inline get_secondary_code_length(::Type{<:GPSL1C_P}) = 1800
