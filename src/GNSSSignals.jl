module GNSSSignals

    using DocStringExtensions, StaticArrays
    using Unitful: Hz

    export
        AbstractGNSSSystem,
        GPSL1,
        GPSL5,
        get_codes,
        get_code_length,
        get_shortest_code_length,
        get_center_frequency,
        get_code_frequency,
        get_code_unsafe,
        get_code,
        get_data_frequency,
        get_code_center_frequency_ratio,
        get_carrier_fast_unsafe,
        get_carrier_vfast_unsafe

    abstract type AbstractGNSSSystem end

    struct GPSL1 <: AbstractGNSSSystem end

    struct GPSL5 <: AbstractGNSSSystem end

    """
    $(SIGNATURES)

    Reads Int8 encoded codes from a file with filename `filename` (including
    the path). The code length must be provided by `code_length` and the
    number of PRNs by `num_prns`.
    # Examples
    ```julia-repl
    julia> read_in_codes("/data/gpsl1codes.bin", 32, 1023)
    ```
    """
    function read_in_codes(filename, num_prns, code_length)
        open(filename) do file_stream
            read!(file_stream, Array{Int8}(undef, code_length, num_prns))
        end
    end

    include("gpsl1.jl")
    include("gpsl5.jl")
    include("carrier.jl")
    include("common.jl")

end
