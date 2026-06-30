"""
    GPSL5I{C} <: AbstractGNSSSignal{C}

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
struct GPSL5I{C<:AbstractMatrix} <: AbstractGNSSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GPSL5Q{C} <: AbstractGNSSSignal{C}

GPS L5-Q signal (the quadrature, dataless pilot component of GPS L5).

BPSK-modulated 10230-chip primary code at 10.23 Mcps on the L5 band
(1176.45 MHz), with a 20-bit Neuman-Hofman secondary code (NH20) overlaying
the pilot channel for a 20 ms tiered code. As the pilot component it carries
no navigation data, so [`get_data_frequency`](@ref) returns 0 Hz.

The Q5 primary code uses the same XA/XB shift-register generator as
[`GPSL5I`](@ref) (IS-GPS-705 §3.2.1.1, `gen_l5_code`) and differs only in the
per-PRN initial XB code state (`INITIAL_XB_CODE_STATES_Q5`). PRNs 1-37 are
supported, matching [`GPSL5I`](@ref).

# Example
```julia
gpsl5q = GPSL5Q()
get_code_length(gpsl5q)            # 10230
get_secondary_code_length(gpsl5q)  # 20
get_data_frequency(gpsl5q)         # 0 Hz
get_band(gpsl5q)                   # L5()
```
"""
struct GPSL5Q{C<:AbstractMatrix} <: AbstractGNSSSignal{C}
    codes::C
    # Cached element-wise negation of `codes`, exactly as in `GPSL5I`: the
    # NH20 secondary chips are ±1, so for the half of primary periods where
    # the secondary chip is -1 we read from `negated_codes` instead of
    # multiplying every chip by `sec_val`, restoring the fused
    # load-broadcast pattern in the hot inner loop.
    negated_codes::C
end

#= GPS L5 primary code generation (IS-GPS-705 §3.2.1.1).

Each L5 primary code is the modulo-2 sum of two 13-stage shift registers, XA
and XB, truncated to 10230 chips. XA is common to every PRN (all-ones start,
and short-cycled: reset to all-ones one chip before its natural period); XB
carries a per-PRN initial state. The I5 and Q5 codes share these identical XA
and XB generators and differ only in the per-PRN initial XB code state, so
`gen_l5_code` serves both — the component is selected purely by which
initial-state table is passed in (`INITIAL_XB_CODE_STATES_I5` for I5,
`INITIAL_XB_CODE_STATES_Q5` for Q5). =#

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

The I5 and Q5 codes share the same XA and XB shift-register generators
(IS-GPS-705 §3.2.1.1) and differ only in the per-PRN initial XB code state, so
this generator is used by both [`GPSL5I`](@ref) and [`GPSL5Q`](@ref) — the
component is selected purely by which initial-state table is passed in
(`INITIAL_XB_CODE_STATES_I5` for I5, `INITIAL_XB_CODE_STATES_Q5` for Q5).
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

#= These are the initial XB code states for the I5 code.
INITIAL_XB_CODE_STATES_I5[1] is a 13-chip array holding the shift-register
values for PRN 1; INITIAL_XB_CODE_STATES_I5[3][4] is the 4th shift-register
stage for the GPS signal with PRN number 3. =#
const INITIAL_XB_CODE_STATES_I5 = [             #sat PRN number
    [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0],    #01
    [1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1],    #02
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],    #03
    [1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0],    #04
    [1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1],    #05
    [0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0],    #06
    [1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1],    #07
    [1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0],    #08
    [1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1],    #09
    [0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0],    #10
    [0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0],    #11
    [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1],    #12
    [0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],    #13
    [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1],    #14
    [0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0],    #15
    [0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1],    #16
    [0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1],    #17
    [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0],    #18
    [1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1],    #19
    [0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1],    #20
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],    #21
    [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1],    #22
    [1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],    #23
    [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0],    #24
    [1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1],    #25
    [1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0],    #26
    [0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],    #27
    [0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],    #28
    [0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1],    #29
    [1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],    #30
    [0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0],    #31
    [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1],    #32
    [1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1],    #33
    [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1],    #34
    [1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],    #35
    [1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0],    #36
    [0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0],    #37
]

#= These are the initial XB code register states for the Q5 code, in the same
format as `INITIAL_XB_CODE_STATES_I5`: `INITIAL_XB_CODE_STATES_Q5[3][4]`
is the 4th shift-register stage for the satellite with PRN number 3.

The Q5 code shares the I5 XA/XB generators and differs only in these initial
XB states (IS-GPS-705 §3.2.1.2). The values equal the all-ones XB register
advanced by the IS-GPS-705 Q5 "XB code advance" for each PRN (the same advances
published by GNSS-SDR's `GPS_L5q_INIT_REG` and PocketSDR's `L5Q_XB_adv`); the
analogous derivation reproduces the I5 table in `INITIAL_XB_CODE_STATES_I5`
exactly. =#
const INITIAL_XB_CODE_STATES_Q5 = [             #sat PRN number
    [1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0],    #01
    [0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0],    #02
    [1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1],    #03
    [0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0],    #04
    [0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0],    #05
    [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1],    #06
    [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1],    #07
    [0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0],    #08
    [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1],    #09
    [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0],    #10
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1],    #11
    [0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1],    #12
    [0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1],    #13
    [1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],    #14
    [1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1],    #15
    [1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1],    #16
    [1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0],    #17
    [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0],    #18
    [0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1],    #19
    [1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],    #20
    [0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],    #21
    [0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0],    #22
    [1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1],    #23
    [0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1],    #24
    [0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1],    #25
    [0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0],    #26
    [1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0],    #27
    [1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0],    #28
    [0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0],    #29
    [1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1],    #30
    [0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1],    #31
    [1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0],    #32
    [1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0],    #33
    [1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0],    #34
    [0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1],    #35
    [0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1],    #36
    [0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1],    #37
]

function read_gpsl5i_codes()
    mapreduce(sat -> gen_l5_code(INITIAL_XB_CODE_STATES_I5[sat]), hcat, 1:37)
end

function read_gpsl5q_codes()
    mapreduce(sat -> gen_l5_code(INITIAL_XB_CODE_STATES_Q5[sat]), hcat, 1:37)
end

function GPSL5I()
    codes = widen_codes_to_storage(read_gpsl5i_codes())
    lut = build_signal_lut(get_modulation(GPSL5I), codes, _gpsl5i_secondary_code())
    GPSL5I(codes, lut)
end

function GPSL5Q()
    codes = widen_codes_to_storage(read_gpsl5q_codes())
    GPSL5Q(codes, .-codes)
end

# Shared interface (modulation, band, frequencies).

get_modulation(::Type{<:GPSL5I}) = LOC()
@inline get_modulation(::GPSL5I) = LOC()
get_modulation(::Type{<:GPSL5Q}) = LOC()
@inline get_modulation(::GPSL5Q) = LOC()

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

# Examples
```julia-repl
julia> get_band(GPSL5I())
L5()
```
"""
@inline get_band(::GPSL5I) = L5()
@inline get_band(::GPSL5Q) = L5()

"""
$(SIGNATURES)

Get the human-readable signal name.

# Examples
```julia-repl
julia> get_signal_name(GPSL5I())
"GPS L5-I"
```
"""
get_signal_name(::GPSL5I) = "GPS L5-I"
get_signal_name(::GPSL5Q) = "GPS L5-Q"

"""
$(SIGNATURES)

Get the code length for GPS L5 (10230 chips, both components).

# Examples
```julia-repl
julia> get_code_length(GPSL5I())
10230
```
"""
@inline get_code_length(::GPSL5I) = 10230
@inline get_code_length(::GPSL5Q) = 10230

"""
$(SIGNATURES)

Get the code chipping rate for GPS L5 (10.23 MHz, both components).

# Examples
```julia-repl
julia> get_code_frequency(GPSL5I())
10230000 Hz
```
"""
@inline get_code_frequency(::GPSL5I) = 10_230_000Hz
@inline get_code_frequency(::GPSL5Q) = 10_230_000Hz

"""
$(SIGNATURES)

Get the data symbol rate for GPS L5-I.

The L5-I channel carries the CNAV navigation message at 100 sps (50 bps with
rate-1/2 convolutional coding); `get_data_frequency` returns the broadcast
symbol rate, matching the convention used across this package (see
[`GPSL1C_D`](@ref)).

# Returns
- `Frequency`: 100 Hz

# Examples
```julia-repl
julia> get_data_frequency(GPSL5I())
100 Hz
```
"""
@inline get_data_frequency(::GPSL5I) = 100Hz

"""
$(SIGNATURES)

Get the data symbol rate for GPS L5-Q.

The Q5 component is a dataless pilot, so its data frequency is 0 Hz.

# Returns
- `Frequency`: 0 Hz

# Examples
```julia-repl
julia> get_data_frequency(GPSL5Q())
0 Hz
```
"""
@inline get_data_frequency(::GPSL5Q) = 0Hz

"""
$(SIGNATURES)

Get the secondary (Neuman-Hofman NH10) code for GPS L5-I.

NH10 is shared across all PRNs: every primary code period (1 ms) is XOR'd
with one chip of the 10-bit sequence `0000110101`, mapped to `±1`.

# Returns
- [`SharedSecondaryCode`](@ref) of length 10

# Examples
```julia-repl
julia> get_secondary_code(GPSL5I())
SharedSecondaryCode{10, Int8}((1, 1, 1, 1, -1, -1, 1, -1, 1, -1))
```
"""
@inline function get_secondary_code(::GPSL5I)
    _gpsl5i_secondary_code()
end

# NH10 secondary, shared across PRNs. Factored out so the `GPSL5I` constructor can build the
# embedded `SignalLUT` (which needs the secondary) before an instance exists.
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
one chip of the 20-bit sequence `00000100110101001110` (IS-GPS-705 §3.2.1.2),
mapped to `±1`, giving a 20 ms tiered code.

# Returns
- [`SharedSecondaryCode`](@ref) of length 20

# Examples
```julia-repl
julia> get_secondary_code(GPSL5Q())
SharedSecondaryCode{20, Int8}((1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1))
```
"""
@inline function get_secondary_code(::GPSL5Q)
    SharedSecondaryCode(
        Int8(1), Int8(1), Int8(1), Int8(1), Int8(1),
        Int8(-1), Int8(1), Int8(1), Int8(-1), Int8(-1),
        Int8(1), Int8(-1), Int8(1), Int8(-1), Int8(1),
        Int8(1), Int8(-1), Int8(-1), Int8(-1), Int8(1),
    )
end

# GPSL5I and GPSL5Q pre-negate their primary code matrix (NH10/NH20 secondary
# chips are ±1); the shared `_select_codes_for` fast path for such signals lives
# on `NegatedPrimaryCacheSignal` in `common.jl`.
