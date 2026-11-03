"""
    GPSL2CM{C} <: AbstractGPSSignal{C}

GPS L2 CM signal (the moderate-length, data-carrying component of the GPS L2
civil signal L2C).

BPSK-modulated 10230-chip code at 511.5 kcps on the L2 band (1227.6 MHz),
giving a 20 ms code period (IS-GPS-200N §3.2.1.4). It carries the CNAV
message, transmitted at 50 symbols/s: the symbol rate of its rate-½ convolutional code,
the message carrying 25 bit/s (IS-GPS-200N §3.2.2; see [`get_bit_rate`](@ref)).

In the broadcast L2C signal the CM-code is time-division multiplexed chip-by-chip with
the dataless [`GPSL2CL`](@ref) pilot at 1.023 MHz; like GNSS-SDR and PocketSDR this
implementation models CM on its own at its native 511.5 kcps. There is no
secondary/overlay code.

CM and CL share one degree-27 modular shift-register generator (IS-GPS-200N §3.3.2.4,
`gen_l2c_code`) and differ only in the per-PRN initial register state
(`GPS_L2CM_INITIAL_STATES`) and the short-cycle period. PRNs 1-63 are supported.

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
giving a 1.5 s code period (IS-GPS-200N §3.2.1.5). As the pilot it carries no
navigation data, so [`get_bit_rate`](@ref) returns 0 Hz, and it has no
secondary/overlay code (the 1.5 s code is itself the long code).

In the broadcast L2C signal CL is time-division multiplexed chip-by-chip with
[`GPSL2CM`](@ref) at 1.023 MHz; like GNSS-SDR and PocketSDR this implementation models
CL on its own at its native 511.5 kcps.

CL uses the same generator as [`GPSL2CM`](@ref) (IS-GPS-200N §3.3.2.4, `gen_l2c_code`),
with its own per-PRN initial register states (`GPS_L2CL_INITIAL_STATES`) and the longer
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

The CM- and CL-codes share a degree-27 modular (Galois) shift-register generator
(IS-GPS-200N §3.3.2.4), differing in the per-PRN initial state and the short-cycle length
(10230 for CM, 767250 for CL). At each step the register outputs its least-significant
bit (0 → +1, 1 → −1, the IS-GPS-200N / PocketSDR chip convention), then shifts right and
XORs `GPS_L2C_FEEDBACK_MASK` back in when the output bit was 1.

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

# Build the `(code_length, GPS_L2C_NUM_PRNS)` primary-code matrix from a per-PRN
# initial-state table.
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
    codes = read_gpsl2cm_codes()
    lut = build_signal_lut(get_modulation(GPSL2CM), codes, NoSecondaryCode())
    GPSL2CM(codes, lut)
end

function GPSL2CL()
    codes = read_gpsl2cl_codes()
    lut = build_signal_lut(get_modulation(GPSL2CL), codes, NoSecondaryCode())
    GPSL2CL(codes, lut)
end

# Shared interface (modulation, band, frequencies).

get_modulation(::Type{<:GPSL2CM}) = BPSK()
get_modulation(::Type{<:GPSL2CL}) = BPSK()

# IS-GPS-200N Table 3-III: the L2 P(Y) carrier is the in-phase (I) reference; L2C (the
# CM/CL multiplex) rides the quadrature carrier. §3.3.1.5.1: CNAV Type 10 message bit 273
# signals quadrature (L2C lagging L2P(Y) by 90° → −1/4 cycle) or in-phase (0); with no CNAV
# it is "fixed in phase quadrature". We return the nominal −1/4 cycle.
@inline get_carrier_phase_offset(::Type{<:GPSL2CM}) = -0.25
@inline get_carrier_phase_offset(::Type{<:GPSL2CL}) = -0.25

@inline get_band(::Type{<:GPSL2CM}) = L2()
@inline get_band(::Type{<:GPSL2CL}) = L2()

# CM and CL are chip-by-chip time-multiplexed, so each carries half the average L2C power.
# The split follows from the multiplexing, not a tabulated figure: the ICD gives only the
# combined L2C (−160.0 dBW, IS-GPS-200N, Table 3-Va, IIR-M/IIF worst case), the unit here.
# See [`get_relative_power`](@ref).
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

The L2 CM-code has no secondary/overlay code ([`NoSecondaryCode`](@ref)).
"""
@inline get_secondary_code(::GPSL2CM) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the secondary code for GPS L2 CL.

The L2 CL-code has no secondary/overlay code ([`NoSecondaryCode`](@ref)); the 1.5 s code
is itself the long code.
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
"""
@inline get_bit_rate(::Type{<:GPSL2CM}) = 50Hz

"""
$(SIGNATURES)

Get the bit rate, as tracked, of GPS L2 CL.

The L2 CL-code is a dataless pilot, so its bit rate is 0 Hz.
"""
@inline get_bit_rate(::Type{<:GPSL2CL}) = 0Hz

@inline get_supported_prns(::Type{<:Union{GPSL2CM,GPSL2CL}}) = 1:63
@inline get_secondary_code_length(::Type{<:Union{GPSL2CM,GPSL2CL}}) = 1
