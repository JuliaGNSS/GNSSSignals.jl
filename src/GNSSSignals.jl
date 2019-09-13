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
        get_data_frequency

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

    """
    $(SIGNATURES)

    Get code of type <: `AbstractGNSSSystem` at phase `phase` of prn `prn`.
    ```julia-repl
    julia> get_code(GPSL1, 1200.3, 1)
    ```
    """
    function get_code(::Type{T}, phase, prn::Int) where T <: AbstractGNSSSystem
        get_code_unsafe(T, mod(phase, get_code_length(T)), prn)
    end

    """
    $(SIGNATURES)

    Get code of type <: `AbstractGNSSSystem` at phase `phase` of prn `prn`.
    The phase will not be wrapped by the code length. The phase has to smaller
    than the code length.
    ```julia-repl
    julia> get_code_unsafe(GPSL1, 10.3, 1)
    ```
    """
    function get_code_unsafe(::Type{T}, phase, prn::Int) where T <: AbstractGNSSSystem
        get_code_unsafe(T, floor(Int, phase), prn::Int)
    end

    include("gpsl1.jl")
    include("gpsl5.jl")
    include("carrier.jl")

end
