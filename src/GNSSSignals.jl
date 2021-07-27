module GNSSSignals

    using
        DocStringExtensions,
        Statistics

    using Unitful: Hz

    using CUDA
    const use_gpu = Ref(false)

    export
        AbstractGNSS,
        GPSL1,
        GPSL5,
        GalileoE1B,
        BOCcos,
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


    abstract type AbstractGNSS end
    abstract type AbstractGNSSBOCcos{M, N} <: AbstractGNSS end

    Base.Broadcast.broadcastable(system::AbstractGNSS) = Ref(system)

    function __init__()
        # use_gpu[] = CUDA.functional()
        if use_gpu[]
            @info "Found CUDA, activating GPU signal processing. Call GNSSSignals.use_gpu[] = false to override this. Beware of any created objects, you may need to reconstruct them for the override to take place."
        else
            @info "CUDA not found. Using solely CPU signal processing."
        end
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
        use_gpu[] ? CuArray{Float32}(code_int8) : extend_front_and_back(Int16.(code_int8), code_length)
    end

    function extend_front_and_back(codes, code_length)
        [codes[end - code_length + 1:end,:]; codes; codes[1:code_length,:]]
    end

    include("gps_l1.jl")
    include("gps_l5.jl")
    include("galileo_e1b.jl")
    include("boc.jl")
    include("common.jl")
end
