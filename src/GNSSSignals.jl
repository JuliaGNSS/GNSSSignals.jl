module GNSSSignals

    using DocStringExtensions, DataStructures, StaticArrays
    using Unitful: Hz, ms

    export
        gen_carrier,
        gen_carrier_fast,
        calc_carrier_phase,
        calc_code_phase_unsafe,
        get_code_length,
        gen_code,
        get_code_unsafe,
        calc_code_phase,
        GPSL1,
        GPSL5,
        AbstractGNSSSystem

    abstract type AbstractGNSSSystem end

    struct GPSL1 <: AbstractGNSSSystem
        codes::Array{Int8, 2}
        code_length::Int
        code_freq::typeof(1.0Hz)
        center_freq::typeof(1.0Hz)
        num_prns_per_bit::Int
    end

    struct GPSL5 <: AbstractGNSSSystem
        codes::Array{Int8, 2}
        code_length::Int
        code_freq::typeof(1.0Hz)
        center_freq::typeof(1.0Hz)
        neuman_hofman_code::Vector{Int8}
        code_length_wo_neuman_hofman_code::Int
        num_prns_per_bit::Int
    end

    """
    $(SIGNATURES)

    Reads codes from a file with filename `filename` (including the path).
    The code length is provided by `code_length`.
    # Examples
    ```julia-repl
    julia> read_in_codes("/data/gpsl1codes.bin", 1023)
    ```
    """
    function read_in_codes(filename, code_length)
        file_stats = stat(filename)
        num_prn_codes = floor(Int, file_stats.size / code_length)
        codes = open(filename) do file_stream
            read!(file_stream, Array{Int8}(undef, code_length, num_prn_codes))
        end
    end

    function get_code_length(system::AbstractGNSSSystem)
        system.code_length
    end

    include("gpsl1.jl")
    include("gpsl5.jl")
    include("sampling.jl")

end
