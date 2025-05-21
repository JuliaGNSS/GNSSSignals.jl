using Test, GNSSSignals, Statistics
import Unitful: Hz, MHz, Frequency
import GNSSSignals: BOCsin, BOCcos, CBOC

include("test_codes.jl")
include("modulation.jl")
include("gps_l1.jl")
include("gps_l5.jl")
include("galileo_e1b.jl")
include("common.jl")
