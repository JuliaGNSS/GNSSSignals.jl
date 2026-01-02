using Test, GNSSSignals, Statistics, Aqua
import Unitful: Hz, MHz, Frequency
import GNSSSignals: BOCsin, BOCcos, CBOC

@testset "Aqua.jl" begin
    Aqua.test_all(GNSSSignals; ambiguities=false, deps_compat=(check_extras=false,))
end

include("test_codes.jl")
include("modulation.jl")
include("gps_l1.jl")
include("gps_l5.jl")
include("galileo_e1b.jl")
include("common.jl")
