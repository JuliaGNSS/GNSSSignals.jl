using Test, GNSSSignals, Statistics, Aqua
import Unitful: Hz, MHz, Frequency
import GNSSSignals: BOCsin, BOCcos, CBOC
using CodecZlib: GzipDecompressorStream

@testset "Aqua.jl" begin
    Aqua.test_all(GNSSSignals; deps_compat=(check_extras=false,))
end

include("test_codes.jl")
include("bands.jl")
include("modulation.jl")
include("gps/l1ca.jl")
include("gps/l5i.jl")
include("gps/l1c_codes.jl")
include("gps/l1c_d.jl")
include("gps/l1c_p.jl")
include("galileo/e1b.jl")
include("common.jl")
