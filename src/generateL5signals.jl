"""
Generate sampled code at sample points `samples` with the code frequency `f`, code phase `φ₀` and sampling 
frequency `f_s`. The code is provided by `code`.
# Examples
```julia-repl
julia> gen_sat_code(1:4000, 1023e3, 2, 4e6, [1, -1, 1, 1, 1])
```
"""
    
function gen_L5_I5_code(prn_code, phase)
    XA = CircularBuffer{Int}(13)
    XB = CircularBuffer{Int}(13)
    satellite_code = zeros(10230)
    append!(XA, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    append!(XB, prn_code)
    for i = 1:10230
        if (i == 8190)
            append!(XA, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
            continue
        end
        if (i == 8191)
            append!(XB, prn_code)
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

