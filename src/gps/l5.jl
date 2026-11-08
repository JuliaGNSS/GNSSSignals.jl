"""
    GPSL5I{C} <: AbstractGPSSignal{C}

GPS L5-I signal (the in-phase, data-carrying component of GPS L5).

BPSK-modulated 10230-chip primary code at 10.23 Mcps on the L5 band
(1176.45 MHz), with a 10-bit Neuman-Hofman secondary code (NH10) overlaying
the data channel.

# Example
```julia
gpsl5i = GPSL5I()
get_code_length(gpsl5i)            # 10230
get_secondary_code_length(gpsl5i)  # 10
get_band(gpsl5i)                   # L5()
```
"""
struct GPSL5I{C<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GPSL5Q{C} <: AbstractGPSSignal{C}

GPS L5-Q signal (the quadrature, dataless pilot component of GPS L5).

BPSK-modulated 10230-chip primary code at 10.23 Mcps on the L5 band
(1176.45 MHz), with a 20-bit Neuman-Hofman secondary code (NH20) overlaying
the pilot channel for a 20 ms tiered code. As the pilot component it carries
no navigation data, so [`get_bit_rate`](@ref) returns 0 Hz.

The Q5 primary code uses the same generator as [`GPSL5I`](@ref) (IS-GPS-705J §3.2.1.1,
`gen_l5_code`) with its own per-PRN initial XB states (`INITIAL_XB_CODE_STATES_Q5`).
PRNs 1-37 are supported.

# Example
```julia
gpsl5q = GPSL5Q()
get_code_length(gpsl5q)            # 10230
get_secondary_code_length(gpsl5q)  # 20
get_bit_rate(gpsl5q)               # 0 Hz
get_band(gpsl5q)                   # L5()
```
"""
struct GPSL5Q{C<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

#= GPS L5 primary code generation (IS-GPS-705J §3.3.2.2).

Each L5 primary code is the modulo-2 sum of two 13-stage shift registers, XA
and XB, truncated to 10230 chips. XA is common to every PRN (all-ones start,
short-cycled: reset to all-ones one chip before its natural period); XB
carries a per-PRN initial state. I5 and Q5 differ only in that initial XB state,
so `gen_l5_code` serves both, selected by the table passed in. =#

"""
$(SIGNATURES)

Takes the status of the registers as an Int `register`, and an array of register
`indices` to calculate and return the register output and the new register value.

```julia-repl
julia> output, register = shift_register(8191, [9, 10, 12, 13])
(1, 4095)
```
"""
function shift_register(register, indices)
    update = 0
    for i in indices
        update = update ⊻ ((register >> (13 - i)) & 1)
    end
    (register & 1), (register >> 1) + 2^12 * update
end

"""
$(SIGNATURES)

Calculate a GPS L5 PRN `satellite_code` for the initial XB register states
`initial_xb_code_states`.

Used by both [`GPSL5I`](@ref) and [`GPSL5Q`](@ref) (IS-GPS-705J §3.2.1.1), selected by the
initial-state table passed in (`INITIAL_XB_CODE_STATES_I5` or `INITIAL_XB_CODE_STATES_Q5`).
```julia-repl
julia> initial_states_PRN_num_1_I = [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0]
julia> prn_code_sat_1_I_signal = gen_l5_code(initial_states_PRN_num_1_I)
```
"""
function gen_l5_code(initial_xb_code_states)
    XA = 8191 # all-ones start state: thirteen 1-bits filling the 13-stage XA register
    XB = initial_xb_code_states' * [4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1]
    satellite_code = zeros(Int8, 10230)
    XA_indices = [9, 10, 12, 13]
    XB_indices = [1, 3, 4, 6, 7, 8, 12, 13]
    for i = 1:10230
        output_xa, XA = shift_register(XA, XA_indices)
        output_xb, XB = shift_register(XB, XB_indices)
        satellite_code[i] = 2 * (output_xa ⊻ output_xb) - 1
        if (i == 8190)
            XA = 8191
        end
    end
    return satellite_code
end

function read_gpsl5i_codes()
    mapreduce(prn -> gen_l5_code(INITIAL_XB_CODE_STATES_I5[prn]), hcat, 1:37)
end

function read_gpsl5q_codes()
    mapreduce(prn -> gen_l5_code(INITIAL_XB_CODE_STATES_Q5[prn]), hcat, 1:37)
end

function GPSL5I()
    codes = read_gpsl5i_codes()
    lut = build_signal_lut(get_modulation(GPSL5I), codes, _gpsl5i_secondary_code())
    GPSL5I(codes, lut)
end

function GPSL5Q()
    codes = read_gpsl5q_codes()
    lut = build_signal_lut(get_modulation(GPSL5Q), codes, _gpsl5q_secondary_code())
    GPSL5Q(codes, lut)
end

# Shared interface (modulation, band, frequencies).

get_modulation(::Type{<:GPSL5I}) = BPSK()
get_modulation(::Type{<:GPSL5Q}) = BPSK()

# L5I data rides the in-phase carrier (default 0); L5Q pilot is in quadrature.
# IS-GPS-705J §3.3.1.5: the Q5 carrier lags I5 by 90° (cos(ωt−π/2) = sin(ωt)), so
# −1/4 cycle. Opposite to Galileo E5aQ, whose ICD (Eq. 1, I·cos − Q·sin) puts the pilot
# at +1/4 cycle.
@inline get_carrier_phase_offset(::Type{<:GPSL5Q}) = -0.25

@inline get_band(::Type{<:GPSL5I}) = L5()
@inline get_band(::Type{<:GPSL5Q}) = L5()

# The ICD tabulates I5 and Q5 at the same −157.9 dBW (IS-GPS-705J, Table 3-III, Block IIF
# worst case), a 50/50 split. The block dependence drops out: GPS III/IIIF is 0.9 dB
# stronger at −157.0 but lifts both components together. See [`get_relative_power`](@ref).
@inline get_relative_power(::Type{<:GPSL5I}) = 0.5
@inline get_relative_power(::Type{<:GPSL5Q}) = 0.5

@inline get_signal_name(::Type{<:GPSL5I}) = "GPS L5-I"
@inline get_signal_name(::Type{<:GPSL5Q}) = "GPS L5-Q"

"""
$(SIGNATURES)

Get the code length for GPS L5 (10230 chips, both components).

# Examples
```julia-repl
julia> get_code_length(GPSL5I())
10230
```
"""
@inline get_code_length(::Type{<:GPSL5I}) = 10230
@inline get_code_length(::Type{<:GPSL5Q}) = 10230

"""
$(SIGNATURES)

Get the nominal chip rate of GPS L5 (10.23 MHz, both components).

# Examples
```julia-repl
julia> get_chip_rate(GPSL5I())
10230000 Hz
```
"""
@inline get_chip_rate(::Type{<:GPSL5I}) = 10_230_000Hz
@inline get_chip_rate(::Type{<:GPSL5Q}) = 10_230_000Hz

"""
$(SIGNATURES)

Get the bit rate, as tracked, of GPS L5-I.

The L5-I channel carries the CNAV navigation message, tracked at 100 bit/s: the symbol
rate of its rate-½ convolutional code, the message carrying 50 bit/s.

# Examples
```julia-repl
julia> get_bit_rate(GPSL5I())
100 Hz
```
"""
@inline get_bit_rate(::Type{<:GPSL5I}) = 100Hz

"""
$(SIGNATURES)

Get the bit rate, as tracked, of GPS L5-Q.

The Q5 component is a dataless pilot, so its bit rate is 0 Hz.

# Examples
```julia-repl
julia> get_bit_rate(GPSL5Q())
0 Hz
```
"""
@inline get_bit_rate(::Type{<:GPSL5Q}) = 0Hz

"""
$(SIGNATURES)

Get the secondary (Neuman-Hofman NH10) code for GPS L5-I.

NH10 is shared across all PRNs: every primary code period (1 ms) is XOR'd with one chip
of the 10-bit sequence `0000110101`, mapped to `±1`. Returns a
[`SharedSecondaryCode`](@ref) of length 10.

# Examples
```julia-repl
julia> get_secondary_code(GPSL5I())
SharedSecondaryCode{10, Int8}((1, 1, 1, 1, -1, -1, 1, -1, 1, -1))
```
"""
@inline function get_secondary_code(::GPSL5I)
    _gpsl5i_secondary_code()
end

# NH10 secondary, factored out so the `GPSL5I` constructor can build the `SignalLUT` before
# an instance exists.
@inline function _gpsl5i_secondary_code()
    SharedSecondaryCode(
        Int8(1), Int8(1), Int8(1), Int8(1), Int8(-1),
        Int8(-1), Int8(1), Int8(-1), Int8(1), Int8(-1),
    )
end

"""
$(SIGNATURES)

Get the secondary (Neuman-Hofman NH20) code for GPS L5-Q.

NH20 is shared across all PRNs: every primary code period (1 ms) is XOR'd with
one chip of the 20-bit sequence `00000100110101001110` (IS-GPS-705J §3.3.2.3),
mapped to `±1`, giving a 20 ms tiered code. Returns a [`SharedSecondaryCode`](@ref) of
length 20.

# Examples
```julia-repl
julia> get_secondary_code(GPSL5Q())
SharedSecondaryCode{20, Int8}((1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1))
```
"""
@inline function get_secondary_code(::GPSL5Q)
    _gpsl5q_secondary_code()
end

# NH20 secondary, factored out so the `GPSL5Q` constructor can build the `SignalLUT` before
# an instance exists.
@inline function _gpsl5q_secondary_code()
    SharedSecondaryCode(
        Int8(1), Int8(1), Int8(1), Int8(1), Int8(1),
        Int8(-1), Int8(1), Int8(1), Int8(-1), Int8(-1),
        Int8(1), Int8(-1), Int8(1), Int8(-1), Int8(1),
        Int8(1), Int8(-1), Int8(-1), Int8(-1), Int8(1),
    )
end

@inline get_supported_prns(::Type{<:GPSL5I}) = 1:37
@inline get_secondary_code_length(::Type{<:GPSL5I}) = 10

@inline get_supported_prns(::Type{<:GPSL5Q}) = 1:37
@inline get_secondary_code_length(::Type{<:GPSL5Q}) = 20
