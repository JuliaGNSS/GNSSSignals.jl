module GNSSSignals

    using Core: toInt16
    using DocStringExtensions
    using FixedPointNumbers
    using Statistics
    using Unitful: Frequency, Hz

    import Base.show

    export
        AbstractGNSS,
        GPSL1,
        GPSL5,
        GalileoE1B,
        gen_code!,
        gen_code,
        get_codes,
        get_code_type,
        get_code_length,
        get_secondary_code_length,
        get_center_frequency,
        get_code_frequency,
        get_code,
        get_code_unsafe,
        get_data_frequency,
        get_code_center_frequency_ratio,
        get_subcarrier_frequency,
        get_code_spectrum,
        get_system_string,
        min_bits_for_code_length,
        get_modulation


    abstract type AbstractGNSS{C} end

    Base.Broadcast.broadcastable(system::AbstractGNSS) = Ref(system)

    Base.show(io::IO, x::AbstractGNSS) = print("$(typeof(x))()")

    """
    $(SIGNATURES)

    Reads Int8 encoded codes from a file with filename `filename` (including the path). The
    code length must be provided by `code_length` and the number of PRNs by `num_prns`.
    # Examples
    ```julia-repl
    julia> read_in_codes("/data/gpsl1codes.bin", 32, 1023)
    ```
    """
    function read_in_codes(type, filename, num_prns, code_length)
        open(filename) do file_stream
            read!(file_stream, Array{type}(undef, code_length, num_prns))
        end
    end


    include("modulation.jl")
    include("gps_l1.jl")
    include("gps_l5.jl")
    include("galileo_e1b.jl")
    include("common.jl")
end
