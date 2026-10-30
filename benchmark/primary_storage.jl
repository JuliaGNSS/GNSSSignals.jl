# Compare primary-code matrix storage while sharing the exact same baked LUT.
# Run from the repository root with BenchmarkTools available in the Julia environment:
#   julia --project=. benchmark/primary_storage.jl
using BenchmarkTools
using GNSSSignals
using Unitful: Hz

const N = 65_536

# Fractional phases also exercise the subcarrier for the CBOC case.
function scalar_sum(signal)
    result = 0.0
    for i in 0:N-1
        result += get_code(signal, i + 0.25, 1)
    end
    result
end

function compare_storage(S, fs)
    original = S()
    small = S(Int8.(original.codes), original.lut)
    large = S(Int16.(original.codes), original.lut)
    @assert small.lut === large.lut
    @assert all(get_code(small, i + 0.25, 1) == get_code(large, i + 0.25, 1) for i in 0:N-1)

    small_out, large_out = zeros(Int8, N), zeros(Int8, N)
    gen_code!(small_out, small, 1, fs)
    gen_code!(large_out, large, 1, fs)
    @assert small_out == large_out
    @assert scalar_sum(small) == scalar_sum(large)

    println("\n", nameof(S), " (", N, " samples, ", fs, ")")
    for (label, signal, out) in (("Int8", small, small_out), ("Int16", large, large_out))
        scalar = @benchmarkable scalar_sum($signal) evals = 1 samples = 100
        generated = @benchmarkable gen_code!($out, $signal, $1, $fs) evals = 1 samples = 100
        scalar_trial = run(scalar; seconds = 1)
        generated_trial = run(generated; seconds = 1)
        println(label, ": primary matrix = ", Base.summarysize(signal.codes), " bytes; ",
            "scalar minimum = ", minimum(scalar_trial).time / 1e6, " ms; ",
            "gen_code! minimum = ", minimum(generated_trial).time / 1e6, " ms")
    end
end

for (S, fs) in ((GPSL1CA, 5e6Hz), (GPSL2CL, 2.5e6Hz), (GalileoE1B_CBOC, 15e6Hz))
    compare_storage(S, fs)
end
