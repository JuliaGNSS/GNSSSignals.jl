using BenchmarkTools
using GNSSSignals
using Unitful: Hz

const SUITE = BenchmarkGroup()
SUITE["code"] = BenchmarkGroup(["Code"])

num_samples = 2000
gpsl1 = GPSL1()
sampling_frequency = 4e6Hz
sampled_code = zeros(Int32, num_samples)

SUITE["code"]["sampling"] =
    @benchmarkable gen_code!($sampled_code, $gpsl1, $1, $sampling_frequency) evals = 10 samples =
        10000