"""
Generate L5 satellite code withe the  `initial_xb_code_states` with the code phase `φ₀`.
# Examples
```julia-repl
julia> gen_sat_code(1:4000, 1023e3, 2, 4e6, [1, -1, 1, 1, 1])
```
"""
    
function gen_L5_I5_code(initial_xb_code_states, φ₀)
    # ToDo include phase offset into calculations
    XA = CircularBuffer{Int}(13)
    XB = CircularBuffer{Int}(13)
    satellite_code = zeros(10230)
    append!(XA, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    append!(XB, initial_xb_code_states)
    for i = 1:10230
        if (i == 8190)
            append!(XA, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
            continue
        end
        if (i == 8191)
            append!(XB, initial_xb_code_states)
        end
        xa = 1 + XA[9] + XA[10] + XA[12] + XA[13]
        xb = 1 + XB[1] + XB[3] + XB[4] + XB[6] + XB[7] + XB[8] + XB[12] + XB[13]
        new_value = (xa+xb) % 2
        unshift!(XB, new_value)
        unshift!(XA, new_value)
        satellite_code[i] = new_value
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
function init_gpsl5_code(sat)
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
    satellite_code = gen_L5_I5_code(initial_xb_code_states[sat], phase)
    gen_sampled_code(samples, f, φ₀, f_s) = gen_sat_code(samples, f, φ₀, f_s, satellite_code)
    get_sampled_code_phase(sample, f, φ₀, f_s) = get_sat_code_phase(sample, f, φ₀, f_s, code_length)
    gen_sampled_code, get_sampled_code_phase
end