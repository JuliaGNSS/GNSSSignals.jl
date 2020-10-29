module GNSSSignals

    using
        DocStringExtensions,
        LoopVectorization,
        StructArrays,
        Statistics,
        FixedPointSinCosApproximations,
        CUDA

    using Unitful: Hz

    export
        AbstractGNSSSystem,
        GPSL1,
        GPSL5,
        GalileoE1B,
        get_codes,
        get_code_length,
        get_secondary_code_length,
        get_center_frequency,
        get_code_frequency,
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
        min_bits_for_code_length,
        length

    const use_gpu = CUDA.functional()

    abstract type AbstractGNSSSystem{T} end

    struct GPSL1{T} <: AbstractGNSSSystem{T} 
        codes::T
    end

    function GPSL1()
        GPSL1(
            extend_front_and_back(read_in_codes(
            joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_gps_l1.bin"),
            37,
            1023
            ))
        )
    end

    struct GPSL5{T} <: AbstractGNSSSystem{T} 
        codes::T
    end

    # function GPSL5()
    #     codes = [                #sat PRN number
    #     [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0],    #01
    #     [1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1],    #02
    #     [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],    #03
    #     [1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0],    #04
    #     [1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1],    #05
    #     [0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0],    #06
    #     [1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1],    #07
    #     [1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0],    #08
    #     [1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1],    #09
    #     [0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0],    #10
    #     [0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0],    #11
    #     [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1],    #12
    #     [0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],    #13
    #     [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1],    #14
    #     [0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0],    #15
    #     [0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1],    #16
    #     [0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1],    #17
    #     [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0],    #18
    #     [1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1],    #19
    #     [0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1],    #20
    #     [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],    #21
    #     [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1],    #22
    #     [1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],    #23
    #     [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0],    #24
    #     [1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1],    #25
    #     [1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0],    #26
    #     [0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],    #27
    #     [0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],    #28
    #     [0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1],    #29
    #     [1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],    #30
    #     [0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0],    #31
    #     [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1],    #32
    #     [1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1],    #33
    #     [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1],    #34
    #     [1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],    #35
    #     [1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0],    #36
    #     [0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0]     #37
    # ]
    #     GPSL5(
    #         use_gpu ? CuArray{ComplexF32}(code_int8) : Int16.(code_int8)
    #     )
    # end

    struct GalileoE1B{T} <: AbstractGNSSSystem{T} 
        codes::T
    end


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
        use_gpu ? CuArray{ComplexF32}(code_int8) : Int16.(code_int8)
    end

    function extend_front_and_back(codes)
        [codes[end, :]'; codes; codes[1,:]'; codes[2,:]']
    end

    include("gps_l1.jl")
    include("gps_l5.jl")
    include("galileo_e1b.jl")
    include("carrier.jl")
    include("common.jl")
end
