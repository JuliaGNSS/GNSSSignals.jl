"""
    GPSL2CM{C} <: AbstractGPSSignal{C}

GPS L2 CM signal (the moderate-length, data-carrying component of the GPS L2
civil signal L2C).

BPSK-modulated 10230-chip code at 511.5 kcps on the L2 band (1227.6 MHz),
giving a 20 ms code period (IS-GPS-200N §3.2.1.4). It carries the CNAV
message, tracked at 50 bit/s: the symbol rate of its rate-½ convolutional code,
the message carrying 25 bit/s (IS-GPS-200N §3.2.2; see [`get_bit_rate`](@ref)).

In the broadcast L2C signal the L2 CM-code is time-division multiplexed
chip-by-chip with the dataless [`GPSL2CL`](@ref) pilot at 1.023 MHz; like
GNSS-SDR and PocketSDR this implementation models the CM component on its own
at its native 511.5 kcps chip rate, which is what a receiver replicates when
acquiring and tracking L2 CM. There is no secondary/overlay code.

The CM and CL codes share one degree-27 modular shift-register generator
(IS-GPS-200N §3.3.2.4, `gen_l2c_code`) and differ only in the per-PRN initial
register state (`GPS_L2CM_INITIAL_STATES`) and the short-cycle period. PRNs
1-63 are supported.

# Example
```julia
gpsl2cm = GPSL2CM()
get_code_length(gpsl2cm)      # 10230
get_chip_rate(gpsl2cm)        # 511500 Hz
get_band(gpsl2cm)             # L2()
```
"""
struct GPSL2CM{C<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GPSL2CL{C} <: AbstractGPSSignal{C}

GPS L2 CL signal (the long, dataless pilot component of the GPS L2 civil
signal L2C).

BPSK-modulated 767250-chip code at 511.5 kcps on the L2 band (1227.6 MHz),
giving a 1.5 s code period (IS-GPS-200N §3.2.1.5). As the pilot component it
carries no navigation data, so [`get_bit_rate`](@ref) returns 0 Hz, and
it has no secondary/overlay code (the 1.5 s code is itself the long code).

In the broadcast L2C signal the dataless CL-code is time-division multiplexed
chip-by-chip with the [`GPSL2CM`](@ref) data component at 1.023 MHz; like
GNSS-SDR and PocketSDR this implementation models the CL component on its own
at its native 511.5 kcps chip rate.

The CL code uses the same degree-27 modular shift-register generator as
[`GPSL2CM`](@ref) (IS-GPS-200N §3.3.2.4, `gen_l2c_code`) and differs only in
the per-PRN initial register state (`GPS_L2CL_INITIAL_STATES`) and the longer
767250-chip short-cycle period. PRNs 1-63 are supported.

# Example
```julia
gpsl2cl = GPSL2CL()
get_code_length(gpsl2cl)      # 767250
get_bit_rate(gpsl2cl)         # 0 Hz
get_band(gpsl2cl)             # L2()
```
"""
struct GPSL2CL{C<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
$(SIGNATURES)

Generate one GPS L2 civil code (CM or CL) from its initial shift-register
state `initial_state`, returning a length-`code_length` `Vector{Int8}` of ±1
chips.

The L2 CM- and L2 CL-codes share a single degree-27 modular (Galois)
shift-register generator (IS-GPS-200N §3.3.2.4); they differ only in the
per-PRN initial state and the number of chips before the register is
short-cycled (10230 for CM, 767250 for CL). At each step the register outputs
its least-significant bit (0 → +1, 1 → −1, the IS-GPS-200N / PocketSDR chip
convention), then shifts right and XORs the feedback mask
`GPS_L2C_FEEDBACK_MASK` back in when the output bit was 1.

```julia-repl
julia> code = gen_l2c_code(GPS_L2CM_INITIAL_STATES[1], GPS_L2CM_CODE_LENGTH);

julia> length(code)
10230
```
"""
function gen_l2c_code(initial_state, code_length)
    register = UInt32(initial_state)
    code = Vector{Int8}(undef, code_length)
    @inbounds for i = 1:code_length
        output = register & one(register)
        code[i] = iszero(output) ? Int8(1) : Int8(-1)
        register = (register >> 1) ⊻ (GPS_L2C_FEEDBACK_MASK * output)
    end
    return code
end

# Build the `(code_length, GPS_L2C_NUM_PRNS)` primary-code matrix from a
# per-PRN initial-state table, preallocated and filled column-by-column (the
# pattern used by `_l1c_build_primary_codes`).
function _build_l2c_codes(initial_states::AbstractVector, code_length::Integer)
    codes = Matrix{Int8}(undef, code_length, GPS_L2C_NUM_PRNS)
    for prn = 1:GPS_L2C_NUM_PRNS
        codes[:, prn] = gen_l2c_code(initial_states[prn], code_length)
    end
    return codes
end

read_gpsl2cm_codes() = _build_l2c_codes(GPS_L2CM_INITIAL_STATES, GPS_L2CM_CODE_LENGTH)
read_gpsl2cl_codes() = _build_l2c_codes(GPS_L2CL_INITIAL_STATES, GPS_L2CL_CODE_LENGTH)

function GPSL2CM()
    codes = widen_codes_to_storage(read_gpsl2cm_codes())
    lut = build_signal_lut(get_modulation(GPSL2CM), codes, NoSecondaryCode())
    GPSL2CM(codes, lut)
end

function GPSL2CL()
    codes = widen_codes_to_storage(read_gpsl2cl_codes())
    lut = build_signal_lut(get_modulation(GPSL2CL), codes, NoSecondaryCode())
    GPSL2CL(codes, lut)
end

# Shared interface (modulation, band, frequencies).

get_modulation(::Type{<:GPSL2CM}) = BPSK()
get_modulation(::Type{<:GPSL2CL}) = BPSK()

# IS-GPS-200N Table 3-III designates the L2 P(Y) carrier as the band's in-phase (I)
# reference; the civil L2C component (the CM/CL chip-by-chip multiplex) rides the
# quadrature carrier. §3.3.1.5.1 makes this configurable: CNAV Type 10 message bit 273
# signals either phase quadrature (L2C lagging L2P(Y) by 90° → −1/4 cycle) or in-phase (0);
# with no CNAV it is "fixed in phase quadrature". We return that nominal/default −1/4 cycle.
@inline get_carrier_phase_offset(::Type{<:GPSL2CM}) = -0.25
@inline get_carrier_phase_offset(::Type{<:GPSL2CL}) = -0.25

@inline get_band(::Type{<:GPSL2CM}) = L2()
@inline get_band(::Type{<:GPSL2CL}) = L2()

# CM and CL are chip-by-chip time-multiplexed, so each carries half the average
# L2C power. This 50/50 split follows from the multiplexing rather than a tabulated
# figure — the ICD gives the combined L2C only (−160.0 dBW, IS-GPS-200N, Table 3-Va,
# IIR-M/IIF worst case), and that total is the unit here. See
# [`get_relative_power`](@ref).
@inline get_relative_power(::Type{<:GPSL2CM}) = 0.5
@inline get_relative_power(::Type{<:GPSL2CL}) = 0.5

@inline get_signal_name(::Type{<:GPSL2CM}) = "GPS L2CM"
@inline get_signal_name(::Type{<:GPSL2CL}) = "GPS L2CL"

"""
$(SIGNATURES)

Get the code length for GPS L2 CM (10230 chips, 20 ms at 511.5 kcps).

# Examples
```julia-repl
julia> get_code_length(GPSL2CM())
10230
```
"""
@inline get_code_length(::Type{<:GPSL2CM}) = GPS_L2CM_CODE_LENGTH

"""
$(SIGNATURES)

Get the code length for GPS L2 CL (767250 chips, 1.5 s at 511.5 kcps).

# Examples
```julia-repl
julia> get_code_length(GPSL2CL())
767250
```
"""
@inline get_code_length(::Type{<:GPSL2CL}) = GPS_L2CL_CODE_LENGTH

"""
$(SIGNATURES)

Get the secondary code for GPS L2 CM.

The L2 CM-code has no secondary/overlay code.

# Returns
- [`NoSecondaryCode`](@ref)
"""
@inline get_secondary_code(::GPSL2CM) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the secondary code for GPS L2 CL.

The L2 CL-code has no secondary/overlay code; the 1.5 s code is itself the
long code.

# Returns
- [`NoSecondaryCode`](@ref)
"""
@inline get_secondary_code(::GPSL2CL) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the nominal chip rate of the GPS L2 civil signals (511.5 kHz, both
components).

# Examples
```julia-repl
julia> get_chip_rate(GPSL2CM())
511500 Hz
```
"""
@inline get_chip_rate(::Type{<:GPSL2CM}) = GPS_L2C_CHIP_RATE
@inline get_chip_rate(::Type{<:GPSL2CL}) = GPS_L2C_CHIP_RATE

"""
$(SIGNATURES)

Get the bit rate, as tracked, of GPS L2 CM.

The L2 CM-code carries the CNAV message, tracked at 50 bit/s: the symbol rate of its
rate-½ convolutional code, the message carrying 25 bit/s (IS-GPS-200N §3.2.2).

# Returns
- `Frequency`: 50 Hz
"""
@inline get_bit_rate(::Type{<:GPSL2CM}) = 50Hz

"""
$(SIGNATURES)

Get the bit rate, as tracked, of GPS L2 CL.

The L2 CL-code is a dataless pilot, so its bit rate is 0 Hz.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_bit_rate(::Type{<:GPSL2CL}) = 0Hz
