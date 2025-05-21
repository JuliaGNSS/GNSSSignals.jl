using BenchmarkTools
using GNSSSignals
using Unitful: Hz

const SUITE = BenchmarkGroup()

num_samples = 2000
sampled_code = zeros(Int16, num_samples)

SUITE["code"]["code generation"]["GPSL1"] =
    @benchmarkable gen_code!($sampled_code, $(GPSL1()), $1, $(2e6Hz)) evals = 10 samples =
        10000

SUITE["code"]["code generation"]["GPSL5"] =
    @benchmarkable gen_code!($sampled_code, $(GPSL5()), $1, $(20e6Hz)) evals = 10 samples =
        10000

sampled_code_f32 = zeros(Float32, num_samples)
SUITE["code"]["code generation"]["GalileoE1B"] =
    @benchmarkable gen_code!($sampled_code_f32, $(GalileoE1B()), $1, $(15e6Hz)) evals = 10 samples =
        10000