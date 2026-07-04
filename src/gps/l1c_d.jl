"""
    GPSL1C_D{C} <: AbstractGNSSSignal{C}

GPS L1C data signal (the data-carrying component of L1C, broadcast by
Block III/IIIF satellites alongside the legacy L1 C/A).

BOC(1,1) sine-phased on the L1 band (1575.42 MHz). 10230-chip primary
code at 1.023 Mcps, derived from a Weil-code construction with a 7-chip
expansion inserted at a PRN-specific point (IS-GPS-800G §3.2.2.1.1).
No secondary code; carries the CNAV-2 message at 100 sps (50 bps after
rate-½ LDPC decoding, per IS-GPS-800G §3.2.3 — `get_data_frequency`
exposes the broadcast symbol rate, matching the convention used for all
other GPS / Galileo signals in this package).

# Example
```julia
gpsl1c_d = GPSL1C_D()
get_code_length(gpsl1c_d)   # 10230
get_band(gpsl1c_d)          # L1()
```

PRNs 1-63 are supported; PRNs 64-210 (IS-GPS-800G Table 6.3-1) are not
implemented in this package.
"""
struct GPSL1C_D{C<:AbstractMatrix} <: AbstractGNSSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GPSL1C_D}) = BOCsin(1, 1)
@inline get_modulation(::GPSL1C_D) = BOCsin(1, 1)

"""
$(SIGNATURES)

Get the band the signal is transmitted on.
"""
@inline get_band(::Type{<:GPSL1C_D}) = L1()

"""
$(SIGNATURES)

Get the human-readable signal name.
"""
get_signal_name(::GPSL1C_D) = "GPS L1C-D"

function read_gpsl1c_d_codes()
    _l1c_build_primary_codes(L1C_D_WEIL_INDEX, L1C_D_INSERTION_INDEX)
end

function GPSL1C_D()
    codes = widen_codes_to_storage(read_gpsl1c_d_codes())
    lut = build_signal_lut(get_modulation(GPSL1C_D), codes, NoSecondaryCode())
    GPSL1C_D(codes, lut)
end

"""
$(SIGNATURES)

Get the primary code length for GPS L1C-D.

# Returns
- `Int`: 10230 chips
"""
@inline get_code_length(::Type{<:GPSL1C_D}) = L1C_PRIMARY_LENGTH

"""
$(SIGNATURES)

Get the secondary code for GPS L1C-D.

The data component has no secondary code (the overlay applies only to
the pilot, L1C-P).
"""
@inline get_secondary_code(::GPSL1C_D) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the code chipping rate for GPS L1C-D.

# Returns
- `Frequency`: 1.023 MHz
"""
@inline get_code_frequency(::Type{<:GPSL1C_D}) = 1_023_000Hz

"""
$(SIGNATURES)

Get the symbol rate of the CNAV-2 message broadcast on GPS L1C-D.

Per IS-GPS-800G §3.2.3, each frame of 9 + 1200 + 548 + 24 + 24 - 5 = 1800
LDPC-encoded symbols is broadcast at 100 sps. The post-decode
information bit rate is 50 bps (rate-½ LDPC), but `get_data_frequency`
returns the broadcast symbol rate to stay consistent with the other
signals in this package — e.g. GPS L5I and Galileo E1B both report the
post-FEC channel symbol rate, not the pre-FEC information rate.

# Returns
- `Frequency`: 100 Hz
"""
@inline get_data_frequency(::Type{<:GPSL1C_D}) = 100Hz
