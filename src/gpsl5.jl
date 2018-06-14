
function init_shift_register(init_registers, update_func)
    registers = CircularBuffer{Int}(13)
    append!(registers, init_registers)
    () -> begin
        output = registers[13]
        update = update_func(registers)
        unshift!(registers, update)
        output, registers 
    end
end

function gen_L5_I5_code_v2(initial_xb_code_states)
    reg_xa = init_shift_register(
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        XA -> (XA[9] + XA[10] + XA[12] + XA[13]) % 2
        )
    reg_xb = init_shift_register(
        initial_xb_code_states,
        XB -> (XB[1] + XB[3] + XB[4] + XB[6] + XB[7] + XB[8] + XB[12] + XB[13]) % 2
        )
    satellite_code = zeros(10230)
    for i = 1:10230
        #reg_xa, output_xa = reg_xa()
        #reg_xb, output_xb = reg_xb()
        output_xa, = reg_xa()
        output_xb, = reg_xb()
        satellite_code[i] =  2* ((output_xa + output_xb) % 2) - 1
        if (i == 8190)
            reg_xa = init_shift_register(
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            XA -> (XA[9] + XA[10] + XA[12] + XA[13]) % 2
            )
        end
        if (i == 8191)
            reg_xb = init_shift_register(initial_xb_code_states,
                                 XB -> (XB[1] + XB[3] + XB[4] + XB[6] + XB[7] + XB[8] + XB[12] + XB[13]) % 2
                                 )
        end
        
    end
    return satellite_code
end
"""
Generate L5 PRN satellite code withe the `initial_xb_code_states`.
# Examples
```julia-repl
julia> gen_sat_code(1:4000, 1023e3, 2, 4e6, [1, -1, 1, 1, 1])
```
"""

function gen_L5_I5_code(initial_xb_code_states)
    XA = CircularBuffer{Int}(13)
    XB = CircularBuffer{Int}(13)
    satellite_code = zeros(10230)
    append!(XA, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    append!(XB, initial_xb_code_states)
    for i = 1:10230
        xa = XA[9] + XA[10] + XA[12] + XA[13]
        xb = XB[1] + XB[3] + XB[4] + XB[6] + XB[7] + XB[8] + XB[12] + XB[13]
        output = (XA[13] + XB[13]) % 2
        unshift!(XA, xa % 2)
        unshift!(XB, xb % 2)
        satellite_code[i] = 2 * output - 1
        if (i == 8190)
            append!(XA, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        end
        if (i == 8191)
            append!(XB, initial_xb_code_states)
        end
    end
    return satellite_code
end


"""
$(SIGNATURES)

Returns functions to generate sampled code and code phase for the GPS L1, signal.
# Examples
```julia-repl
julia> gen_gpsl1_code, get_gpsl1_code_phase = init_gpsl1_codes()
julia> gen_gpsl1_code(samples, f, φ₀, f_s, sat)
julia> get_code_phase(sample, f, φ₀, f_s)
```
"""
function init_gpsl5_code()
    code_length = 10230
    #=These are the initial XB Code States for the I5 code, initial_xb_code_states[1] is a 1 3 chip array which represent the shift register values
    initial_xb_code_states[3][4] represents the 4th shift register of the GPS Signal with PRN numver 3  =#
    initial_xb_code_states = [                      #sat PRN number
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
    codes = zeros(37*code_length)
    for i = 1:37
        codes[(code_length * (i-1) + 1):(code_length * i)] = gen_L5_I5_code(initial_xb_code_states[i])
    end
    gen_sampled_code(samples, f, φ₀, f_s, sat) = gen_sat_code(samples, f, φ₀, f_s, codes[(code_length * (sat-1) + 1):(code_length * sat)])
    get_sampled_code_phase(sample, f, φ₀, f_s) = get_sat_code_phase(sample, f, φ₀, f_s, code_length)
    gen_sampled_code, get_sampled_code_phase
end