module GNSSSignals

    using Yeppp, DocStringExtensions, DataStructures

    export 
        gen_carrier,
        get_carrier_phase,
        gen_code,
        calc_code_phase,
        GPSL1,
        GPSL5

    abstract type AbstractGNSSSystem end

    struct GPSL1 <: AbstractGNSSSystem
        codes::Array{Int8, 2}
        code_length::Int
        code_freq::Float64
        center_freq::Float64
    end

    struct GPSL5 <: AbstractGNSSSystem
        codes::Array{Int8, 2}
        code_length::Int
        code_freq::Float64
        center_freq::Float64
        code_length_wo_neuman::Int
    end

    """
    $(SIGNATURES)

    Reads codes from a file with filename `filename` (including the path). The code length is provided 
    by `code_length`.
    # Examples
    ```julia-repl
    julia> read_in_codes("/data/gpsl1codes.bin", 1023)
    ```
    """
    function read_in_codes(filename, code_length)
        file_stats = stat(filename)
        num_prn_codes = floor(Int, file_stats.size / code_length)
        codes = open(filename) do file_stream
            read(file_stream, Int8, code_length, num_prn_codes)
        end
    end 

    include("gpsl1.jl")
    include("gpsl5.jl")
    include("sampling.jl")

end