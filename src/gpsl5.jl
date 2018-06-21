#=These are the initial XB Code States for the I5 code, initial_xb_code_states[1] is a 1 3 chip array which represent the shift register values
initial_xb_code_states[3][4] represents the 4th shift register of the GPS Signal with PRN numver 3  =#
INITIAL_XB_CODE_STATES = [                      #sat PRN number
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
        [0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0]     #37
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
        b = (integer >> (i-1)) & 1
        a[14-i] = b
    end
    a
end

"""
$(SIGNATURES)

Takes the status of the registers as an Int `registers`, and an array of register `indices` to calculate and return the new register values and the register_output.

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

Calculate the gps L5 PRN `satellite_code` for the initial XB register states 'initial_xb_code_states'.
```julia-repl
julia> initial_states_PRN_num_1_I = [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0]
julia> prn_code_sat_1_I_signal = gen_l5_code(initial_states_PRN_num_1_I)
```
"""
function gen_l5_code(initial_xb_code_states)
    XA = 8191 # int with 3 leading zeros and then 13*1
    XB = initial_xb_code_states' * [4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1]
    satellite_code = zeros(Int8, 10230)
    XA_indices = [9, 10, 12, 13]
    XB_indices = [1, 3, 4, 6, 7, 8, 12, 13]
    for i = 1:10230
        output_xa, XA = shift_register(XA, XA_indices)
        output_xb, XB = shift_register(XB, XB_indices)
        satellite_code[i] =  2 * (output_xa ⊻ output_xb) - 1
        if (i == 8190)
            XA = 8191
        end
    end
    return satellite_code
end

"""
$(SIGNATURES)

Returns functions to generate sampled code and code phase for the GPS L5 (I5) signal.
# Examples
```julia-repl
julia> gen_gpsi5_code, get_i5_code_phase = init_gpsl5_i5_codes()
julia> gen_gpsi5_code(samples, f, φ₀, f_s, sat)
julia> get_i5_code_phase(sample, f, φ₀, f_s)
```
"""

function gen_neuman_hofman_sequence(initial_xb_code_states)
    satellite_code = gen_l5_code(initial_xb_code_states) # = satellite_code .⊻ 0
    ones_sat_code = 1 .⊻ satellite_code
    # 10 digit NH code 0000110101
    nh_satellite_code = vcat(satellite_code, satellite_code, satellite_code, satellite_code, #0000
     ones_sat_code, ones_sat_code, satellite_code,  #110
      ones_sat_code, satellite_code, ones_sat_code) #101

end

function init_gpsl5_code()
    code_length = 102300
    codes = zeros(Int8, 102300*37)
    codes = mapreduce(sat -> gen_neuman_hofman_sequence(INITIAL_XB_CODE_STATES[sat]), hcat, 1:37)
    gen_sampled_code(samples, f, φ₀, f_s, sat) = gen_sat_code(samples, f, φ₀, f_s, codes[:,sat])
    get_sampled_code_phase(sample, f, φ₀, f_s) = get_sat_code_phase(sample, f, φ₀, f_s, code_length)
    gen_sampled_code, get_sampled_code_phase
end
