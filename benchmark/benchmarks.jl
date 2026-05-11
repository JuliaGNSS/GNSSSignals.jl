using BenchmarkTools
using GNSSSignals
using Unitful: Hz

const SUITE = BenchmarkGroup()

num_samples = 2000
sampled_code = zeros(Int16, num_samples)

SUITE["code"]["code generation"]["GPSL1CA"] = @benchmarkable gen_code!(
    $sampled_code,
    $(GPSL1CA()),
    $1,
    $(2e6Hz),
    $(1023e3Hz),
    $0.0,
    $0,
    $(Val(2e6Hz)),
) evals = 10 samples = 10000

SUITE["code"]["code generation"]["GPSL5I"] = @benchmarkable gen_code!(
    $sampled_code,
    $(GPSL5I()),
    $1,
    $(20e6Hz),
    $(10230e3Hz),
    $0.0,
    $0,
    $(Val(20e6Hz)),
) evals = 10 samples = 10000

sampled_code_f32 = zeros(Float32, num_samples)
SUITE["code"]["code generation"]["GalileoE1B"] = @benchmarkable gen_code!(
    $sampled_code_f32,
    $(GalileoE1B()),
    $1,
    $(15e6Hz),
    $(1023e3Hz),
    $0.0,
    $0,
    $(Val(15e6Hz)),
) evals = 10 samples = 10000
