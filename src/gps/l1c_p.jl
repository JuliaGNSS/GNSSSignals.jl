"""
    GPSL1C_P{C, M} <: AbstractGNSSSignal{C}

GPS L1C pilot signal (the dataless component of L1C, carrying 75% of
the L1C power per IS-GPS-800G; broadcast by Block III/IIIF satellites).

TMBOC(6,1,4/33) modulated on the L1 band (1575.42 MHz): every 33
primary chips, four positions `{0, 4, 6, 29}` use BOC(6,1) and the
rest use BOC(1,1). 10230-chip primary code at 1.023 Mcps (same Weil
construction as L1C-D, different per-PRN parameters), modulo-2 added
with an 18 s, 1800-bit per-PRN overlay code (IS-GPS-800G §3.2.2.1.2)
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
struct GPSL1C_P{C<:AbstractMatrix, M<:AbstractMatrix} <: AbstractGNSSSignal{C}
    codes::C
    overlay_codes::M    # 1800 × 63 Int8 ±1 matrix, exposed via PerPRNSecondaryCode
    lut::SignalLUT      # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GPSL1C_P}) =
    TMBOC(BOCsin(1, 1), BOCsin(6, 1),
          ntuple(k -> (k - 1) ∈ L1C_TMBOC_BOC6_POSITIONS, Val(L1C_TMBOC_PERIOD)))
@inline get_modulation(s::GPSL1C_P) = get_modulation(typeof(s))

"""
$(SIGNATURES)

Get the band the signal is transmitted on.
"""
@inline get_band(::Type{<:GPSL1C_P}) = L1()

"""
$(SIGNATURES)

Get the human-readable signal name.
"""
get_signal_name(::GPSL1C_P) = "GPS L1C-P"

function read_gpsl1c_p_codes()
    _l1c_build_primary_codes(L1C_P_WEIL_INDEX, L1C_P_INSERTION_INDEX)
end

function GPSL1C_P()
    codes = widen_codes_to_storage(read_gpsl1c_p_codes())
    overlay = _l1c_build_overlay_codes()
    # The 1800-chip per-PRN overlay is far too long to bake (1800·10230·12 ≫ typemax(Int16)),
    # so it stays residual in the SignalLUT and is applied per primary period at gen time.
    lut = build_signal_lut(get_modulation(GPSL1C_P), codes, PerPRNSecondaryCode(overlay))
    GPSL1C_P(codes, overlay, lut)
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

Get the code chipping rate for GPS L1C-P.

# Returns
- `Frequency`: 1.023 MHz
"""
@inline get_code_frequency(::Type{<:GPSL1C_P}) = 1_023_000Hz

"""
$(SIGNATURES)

Get the data bit rate for GPS L1C-P.

The pilot component is dataless.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_data_frequency(::Type{<:GPSL1C_P}) = 0Hz
