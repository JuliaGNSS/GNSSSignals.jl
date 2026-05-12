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
end

get_modulation(::Type{<:GPSL5I}) = LOC()
@inline get_modulation(::GPSL5I) = LOC()

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

#=These are the initial XB Code States for the I5 code,
initial_xb_code_states[1] is a 1 3 chip array which represent the shift
register values initial_xb_code_states[3][4] represents the 4th shift register
of the GPS Signal with PRN numver 3  =#
const INITIAL_XB_CODE_STATES = [                #sat PRN number
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
    [0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0],     #37
]

"""
$(SIGNATURES)

Takes the status of the registers as an `integer` and returns them as an array.
# Examples
```julia-repl
julia> reshape(8190)
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
```
"""
function reshape(integer)
    a = zeros(13)
    for i = 1:13
        b = (integer >> (i - 1)) & 1
        a[14-i] = b
    end
    a
end

"""
$(SIGNATURES)

Takes the status of the registers as an Int `registers`, and an array of register `indices`
to calculate and return the new register values and the register output.

```julia-repl
julia> reshape(8190)
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
julia> output, registers = shift_register(8910,  [9, 10, 12, 13])
julia> reshape(registers)
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
julia> output == 1
true
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

Calculate the GPS L5-I PRN `satellite_code` for the initial XB register states
`initial_xb_code_states`.
```julia-repl
julia> initial_states_PRN_num_1_I = [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0]
julia> prn_code_sat_1_I_signal = gen_l5i_code(initial_states_PRN_num_1_I)
```
"""
function gen_l5i_code(initial_xb_code_states)
    XA = 8191 # int with 3 leading zeros and then 13*1
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
    mapreduce(sat -> gen_l5i_code(INITIAL_XB_CODE_STATES[sat]), hcat, 1:37)
end

function GPSL5I()
    GPSL5I(widen_codes_to_storage(read_gpsl5i_codes()))
end

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
    SharedSecondaryCode((Int8(1), Int8(1), Int8(1), Int8(1), Int8(-1), Int8(-1), Int8(1), Int8(-1), Int8(1), Int8(-1)))
end

"""
$(SIGNATURES)

Get the code length for GPS L5-I.

# Returns
- `Int`: 10230 chips

# Examples
```julia-repl
julia> get_code_length(GPSL5I())
10230
```
"""
@inline function get_code_length(::GPSL5I)
    10230
end

"""
$(SIGNATURES)

Get the code chipping rate for GPS L5-I.

# Returns
- `Frequency`: 10.23 MHz

# Examples
```julia-repl
julia> get_code_frequency(GPSL5I())
10230000 Hz
```
"""
@inline function get_code_frequency(::GPSL5I)
    10_230_000Hz
end

"""
$(SIGNATURES)

Get the data bit rate for GPS L5-I.

# Returns
- `Frequency`: 100 Hz

# Examples
```julia-repl
julia> get_data_frequency(GPSL5I())
100 Hz
```
"""
@inline function get_data_frequency(::GPSL5I)
    100Hz
end
