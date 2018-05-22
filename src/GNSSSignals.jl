module GNSSSignals

    using Yeppp

    export gen_carrier, get_carrier_phase, gen_sat_code, get_sat_code_phase, init_gpsl1_codes

    include("gpsl1.jl")

    function gen_carrier(t, f, φ₀, f_s)
        arg = (2 * π * f / f_s) .* t .+ φ₀
        sin_sig, cos_sig = Yeppp.sin(arg), Yeppp.cos(arg) # use Yeppp for better performance
        return complex.(cos_sig, sin_sig)
        #return cis.(arg) # cis(...) = exp(1im * ...)
    end

    function get_carrier_phase(t, f, φ₀, f_s)
        return mod2pi((2 * π * f / f_s) * t + φ₀)
    end

    function gen_sat_code(t, f, φ₀, f_s, code)
        code_indices = floor.(Int16, f / f_s .* t .+ φ₀)
        code_indices .= 1 .+ mod.(code_indices - 1, length(code))
        @inbounds code_sampled = code[code_indices]
        return code_sampled
    end

    function get_sat_code_phase(t, f, φ₀, f_s, code_length)
        return mod(f / f_s * t + φ₀ + code_length / 2, code_length) - code_length / 2
    end

    function read_in_codes(filename, code_length)
        file_stats = stat(filename)
        num_prn_codes = floor(Int, file_stats.size / code_length)
        codes = open(filename) do file_stream
            float(read(file_stream, Int8, code_length, num_prn_codes))
        end
    end 

end