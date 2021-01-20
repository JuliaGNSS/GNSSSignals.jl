module GNSSSignals

    using
        DocStringExtensions,
        LoopVectorization,
        StructArrays,
        Statistics,
        FixedPointSinCosApproximations

    using Unitful: Frequency, Hz

    export
        AbstractGNSSSystem,
        GPSL1,
        GPSL5,
        GalileoE1B,
        GenericBOC,
        get_codes,
        get_code_length,
        get_secondary_code_length,
        get_center_frequency,
        get_code_frequency,
        get_subcarrier_frequency,
        get_code_spectrum,
        get_code_unsafe,
        get_code,
        get_data_frequency,
        get_code_center_frequency_ratio,
        get_carrier_fast_unsafe,
        get_carrier_vfast_unsafe,
        get_quadrant_size_power,
        get_carrier_amplitude_power,
        fpcarrier_phases!,
        fpcarrier!,
        min_bits_for_code_length

    abstract type AbstractGNSSSystem end

    struct GPSL1 <: AbstractGNSSSystem end

    struct GPSL5 <: AbstractGNSSSystem end

    struct GalileoE1B <: AbstractGNSSSystem end

    struct GenericBOC{T<:AbstractGNSSSystem, m, n} <: AbstractGNSSSystem end


    """
    $(SIGNATURES)

    Reads Int8 encoded codes from a file with filename `filename` (including the path). The
    code length must be provided by `code_length` and the number of PRNs by `num_prns`.
    # Examples
    ```julia-repl
    julia> read_in_codes("/data/gpsl1codes.bin", 32, 1023)
    ```
    """
    function read_in_codes(filename, num_prns, code_length)
        code_int8 = open(filename) do file_stream
            read!(file_stream, Array{Int8}(undef, code_length, num_prns))
        end
        Int16.(code_int8)
    end

    function extend_front_and_back(codes)
        [codes[end, :]'; codes; codes[1,:]'; codes[2,:]']
    end

    include("gps_l1.jl")
    include("gps_l5.jl")
    include("galileo_e1b.jl")
    include("generic_boc.jl")
    include("carrier.jl")
    include("common.jl")
end
