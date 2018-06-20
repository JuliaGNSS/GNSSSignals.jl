module GNSSSignals

    using Yeppp, DocStringExtensions, DataStructures

    export gen_carrier, get_carrier_phase, gen_sat_code, get_sat_code_phase, init_gpsl1_codes, init_gpsl5_i5_codes

    include("gpsl1.jl")
    include("gpsl5.jl")
    include("sampling.jl")

end