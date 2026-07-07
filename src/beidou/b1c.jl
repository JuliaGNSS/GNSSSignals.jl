include("b1c_constants.jl")

"""
    BeiDouB1C_D{C} <: AbstractBeiDouSignal{C}

BeiDou B1C data signal (the data-carrying component of B1C).

Sine-phased BOC(1,1)-modulated 10230-chip primary code at 1.023 Mcps on the L1
band (1575.42 MHz), carrying the B-CNAV1 navigation message at 100 symbols/s.
The data component has no secondary code (BDS-SIS-ICD-B1C-1.0 Table 5-1).

The primary code is a truncated Weil code (BDS-SIS-ICD-B1C-1.0 §5.2.1): a
length-10243 Weil code (modulo-2 sum of a Legendre sequence and a shifted copy)
cyclically truncated to 10230 chips, with a per-SVID phase difference `w` and
truncation point `p` (`B1C_DATA_WEIL_PARAMS`, ICD Table 5-2). PRNs 1-63.

# Example
```julia
b1c_d = BeiDouB1C_D()
get_code_length(b1c_d)   # 10230
get_modulation(b1c_d)    # BOCsin(1, 1)
get_band(b1c_d)          # L1()
```
"""
struct BeiDouB1C_D{C<:AbstractMatrix} <: AbstractBeiDouSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    BeiDouB1C_P{C, M} <: AbstractBeiDouSignal{C}

BeiDou B1C pilot signal (the dataless component of B1C, carrying 3/4 of the B1C
power).

The B1C pilot is specified as QMBOC(6,1,4/33) (BDS-SIS-ICD-B1C-1.0 §4.2): a
BOC(1,1) subcarrier on the in-phase arm and a BOC(6,1) subcarrier on the
quadrature arm, with 4/33 of the pilot power on the BOC(6,1) component. This
implementation generates a pure sine-phased **BOC(1,1)** replica — `get_modulation`
returns `BOCsin(1, 1)` and `gen_code!` emits it (needing only `fs ≥ 2·1.023 MHz`).

Note this is *not* the same kind of choice as [`GalileoE1B`](@ref) →
[`GalileoE1B_BOC11`](@ref). Galileo E1 is CBOC — a *real* weighted sum of
BOC(1,1) and BOC(6,1) — so a faithful real replica exists (the default type
bakes it) and `_BOC11` is a genuine lower-power approximation. QMBOC instead
puts the two components in phase *quadrature*: its faithful replica is complex,
and its real (in-phase) part simply **is** BOC(1,1). A real correlator tracking
the pilot therefore uses BOC(1,1), capturing the in-phase 29/33 of the power;
the quadrature BOC(6,1) arm cannot be captured by a single real replica. Hence
there is no separate "full" vs `_BOC11` pilot type here — BOC(1,1) is already
the in-phase replica, matching PocketSDR (`sdr_code.py` generates B1C pilot as
"BOC(1,1) instead of QMBOC(6,1,4/33)").

The 10230-chip primary code is a truncated Weil code (same construction as
[`BeiDouB1C_D`](@ref), with the pilot per-SVID parameters `B1C_PILOT_WEIL_PARAMS`,
ICD Table 5-3). It is overlaid with an 1800-chip per-SVID secondary code — a
truncated length-3607 Weil code (ICD Table 5-4) — for an 18 s tiered code,
exposed via [`get_secondary_code`](@ref) as a [`PerPRNSecondaryCode`](@ref).
PRNs 1-63.

# Example
```julia
b1c_p = BeiDouB1C_P()
get_code_length(b1c_p)            # 10230
get_secondary_code_length(b1c_p)  # 1800
get_data_frequency(b1c_p)         # 0 Hz
```
"""
struct BeiDouB1C_P{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractBeiDouSignal{C}
    codes::C
    overlay_codes::M    # 1800 × 63 Int8 ±1 matrix, exposed via PerPRNSecondaryCode
    lut::SignalLUT      # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

const B1C_PRIMARY_WEIL_LENGTH = 10243     # Legendre length for the primary codes
const B1C_PRIMARY_LENGTH = 10230
const B1C_SECONDARY_WEIL_LENGTH = 3607    # Legendre length for the pilot secondary
const B1C_SECONDARY_LENGTH = 1800
const B1C_NUM_PRNS = 63

# Build the 10230 × 63 truncated-Weil primary-code matrix for the given per-SVID (w, p) table.
function _read_beidou_b1c_codes(wp_table)
    L = _beidou_legendre(B1C_PRIMARY_WEIL_LENGTH)
    codes = Matrix{Int8}(undef, B1C_PRIMARY_LENGTH, B1C_NUM_PRNS)
    for prn = 1:B1C_NUM_PRNS
        w, p = wp_table[prn]
        codes[:, prn] = _beidou_weil_code(B1C_PRIMARY_WEIL_LENGTH, B1C_PRIMARY_LENGTH, w, p, L)
    end
    codes
end

read_beidou_b1c_d_codes() = _read_beidou_b1c_codes(B1C_DATA_WEIL_PARAMS)
read_beidou_b1c_p_codes() = _read_beidou_b1c_codes(B1C_PILOT_WEIL_PARAMS)

# Build the 1800 × 63 B1C pilot secondary matrix (truncated length-3607 Weil codes).
function _build_beidou_b1c_secondary()
    L = _beidou_legendre(B1C_SECONDARY_WEIL_LENGTH)
    codes = Matrix{Int8}(undef, B1C_SECONDARY_LENGTH, B1C_NUM_PRNS)
    for prn = 1:B1C_NUM_PRNS
        w, p = B1C_PILOT_SECONDARY_WP[prn]
        codes[:, prn] = _beidou_weil_code(B1C_SECONDARY_WEIL_LENGTH, B1C_SECONDARY_LENGTH, w, p, L)
    end
    codes
end

function BeiDouB1C_D()
    codes = widen_codes_to_storage(read_beidou_b1c_d_codes())
    lut = build_signal_lut(get_modulation(BeiDouB1C_D), codes, NoSecondaryCode())
    BeiDouB1C_D(codes, lut)
end

function BeiDouB1C_P()
    codes = widen_codes_to_storage(read_beidou_b1c_p_codes())
    overlay = _build_beidou_b1c_secondary()
    # The 1800-chip per-SVID overlay is far too long to bake, so it stays residual in the
    # SignalLUT and is applied per primary period at gen time.
    lut = build_signal_lut(get_modulation(BeiDouB1C_P), codes, PerPRNSecondaryCode(overlay))
    BeiDouB1C_P(codes, overlay, lut)
end

# Shared interface.

get_modulation(::Type{<:BeiDouB1C_D}) = BOCsin(1, 1)
@inline get_modulation(::BeiDouB1C_D) = BOCsin(1, 1)
get_modulation(::Type{<:BeiDouB1C_P}) = BOCsin(1, 1)   # in-phase (real) replica of QMBOC(6,1,4/33); see docstring
@inline get_modulation(::BeiDouB1C_P) = BOCsin(1, 1)

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

BeiDou B1C shares the L1 carrier frequency (1575.42 MHz), so this returns
[`L1`](@ref) — band identity is by RF, not by ICD label.
"""
@inline get_band(::Type{<:BeiDouB1C_D}) = L1()
@inline get_band(::Type{<:BeiDouB1C_P}) = L1()

"""
$(SIGNATURES)

Get the human-readable signal name.
"""
get_signal_name(::BeiDouB1C_D) = "BeiDou B1C data"
get_signal_name(::BeiDouB1C_P) = "BeiDou B1C pilot"

"""
$(SIGNATURES)

Get the code length for BeiDou B1C (10230 chips, both components).
"""
@inline get_code_length(::Type{<:BeiDouB1C_D}) = B1C_PRIMARY_LENGTH
@inline get_code_length(::Type{<:BeiDouB1C_P}) = B1C_PRIMARY_LENGTH

"""
$(SIGNATURES)

Get the code chipping rate for BeiDou B1C (1.023 MHz, both components).
"""
@inline get_code_frequency(::Type{<:BeiDouB1C_D}) = 1_023_000Hz
@inline get_code_frequency(::Type{<:BeiDouB1C_P}) = 1_023_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for BeiDou B1C data (B-CNAV1, 100 symbols/s).

# Returns
- `Frequency`: 100 Hz
"""
@inline get_data_frequency(::Type{<:BeiDouB1C_D}) = 100Hz

"""
$(SIGNATURES)

Get the data symbol rate for BeiDou B1C pilot (dataless).

# Returns
- `Frequency`: 0 Hz
"""
@inline get_data_frequency(::Type{<:BeiDouB1C_P}) = 0Hz

"""
$(SIGNATURES)

Get the secondary code for BeiDou B1C data.

The B1C data component has no secondary code (BDS-SIS-ICD-B1C-1.0 Table 5-1).

# Returns
- [`NoSecondaryCode`](@ref)
"""
@inline get_secondary_code(::BeiDouB1C_D) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the secondary code for BeiDou B1C pilot.

Every primary period (10 ms) is overlaid with one chip of an 1800-bit per-SVID
truncated Weil code (BDS-SIS-ICD-B1C-1.0 §5.2.2), giving an 18 s tiered code.

# Returns
- [`PerPRNSecondaryCode`](@ref) wrapping the 1800 × 63 overlay matrix
"""
@inline get_secondary_code(s::BeiDouB1C_P) = PerPRNSecondaryCode(s.overlay_codes)
