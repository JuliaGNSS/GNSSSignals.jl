using BenchmarkTools
using GNSSSignals
using Unitful: Hz, ustrip, @u_str

# Strip a Unitful frequency to a plain Float64 Hz value (for computing N).
ustrip_hz(x) = Float64(ustrip(u"Hz", x))

# Use the v2 names when available, fall back to the pre-v2 names so the same
# benchmark script can run against master and this branch. Remove the fallback
# once master is on v2.
const _GPSL1 = isdefined(GNSSSignals, :GPSL1CA) ? GNSSSignals.GPSL1CA : GNSSSignals.GPSL1
const _GPSL5 = isdefined(GNSSSignals, :GPSL5I) ? GNSSSignals.GPSL5I : GNSSSignals.GPSL5

const SUITE = BenchmarkGroup()

# ── original gen_code! rows (legacy fixed-size buffers, kept for continuity) ──
num_samples = 2000
sampled_code = zeros(Int16, num_samples)

SUITE["code"]["code generation"]["GPSL1"] = @benchmarkable gen_code!(
    $sampled_code,
    $(_GPSL1()),
    $1,
    $(2e6Hz),
    $(1023e3Hz),
    $0.0,
    $0,
) evals = 10 samples = 10000

SUITE["code"]["code generation"]["GPSL5"] = @benchmarkable gen_code!(
    $sampled_code,
    $(_GPSL5()),
    $1,
    $(20e6Hz),
    $(10230e3Hz),
    $0.0,
    $0,
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
) evals = 10 samples = 10000

# ─────────────────────────────────────────────────────────────────────────────
# Realistic 1 ms integration: per signal, compare the original `gen_code!`
# (plain signal) against a prebuilt `CodeReplicaLUT` plan + `gen_code!(plan)`.
# N = round(Int, fs·1e-3); fs chosen ≥ fc·subchip_factor per signal. The plan is
# built ONCE outside the timed region (a receiver reuses it across integrations),
# so only the per-call resample is timed.
# ─────────────────────────────────────────────────────────────────────────────

# (name, signal, prn, fs, fc) — fc is the primary chip rate (1.023 MHz family
# except L5I at 10.23 MHz). fs picked above fc·subchip_factor:
#   L1CA   LOC   P=1  → 5 MHz
#   L1C_D  BOC11 P=2  → 5 MHz   (2.046 MHz)
#   E1B_B  BOC11 P=2  → 8 MHz
#   L5I    LOC   P=1  → 40 MHz
#   L1C_P  TMBOC P=12 → 25 MHz  (12.276 MHz)
const _LUT_CASES = let cases = Any[]
    push!(cases, ("GPSL1CA",  _GPSL1(),                  1, 5e6Hz,  1023e3Hz))
    push!(cases, ("GPSL5I",   _GPSL5(),                  1, 40e6Hz, 10230e3Hz))
    isdefined(GNSSSignals, :GPSL1C_D) &&
        push!(cases, ("GPSL1C_D", GNSSSignals.GPSL1C_D(), 1, 5e6Hz,  1023e3Hz))
    isdefined(GNSSSignals, :GalileoE1B_BOC11) &&
        push!(cases, ("GalileoE1B_BOC11", GNSSSignals.GalileoE1B_BOC11(), 1, 8e6Hz, 1023e3Hz))
    isdefined(GNSSSignals, :GPSL1C_P) &&
        push!(cases, ("GPSL1C_P", GNSSSignals.GPSL1C_P(), 1, 25e6Hz, 1023e3Hz))
    cases
end

for (name, signal, prn, fs, fc) in _LUT_CASES
    N = round(Int, ustrip_hz(fs) * 1e-3)
    out16 = zeros(Int16, N)
    SUITE["code"]["1ms gen_code! original"][name] = @benchmarkable gen_code!(
        $out16, $signal, $prn, $fs, $fc, $0.0, $0,
    ) evals = 10 samples = 1000

    # Plan built ONCE, outside the timed region (receiver reuses it).
    if isdefined(GNSSSignals, :CodeReplicaLUT)
        plan = GNSSSignals.CodeReplicaLUT(signal, prn)
        out8 = zeros(Int8, N)
        SUITE["code"]["1ms gen_code! plan"][name] = @benchmarkable gen_code!(
            $out8, $plan, $fs, $fc, $0.0, $0,
        ) evals = 10 samples = 1000
    end

    # Continuing generator: DDA init paid ONCE in the CodeGeneratorLUT constructor (outside
    # the timed region), then `gen_code!(out, gen)` per integration continues the state with
    # no rate setup / re-init and a single-stream+scalar tail — the fast repeated path.
    if isdefined(GNSSSignals, :CodeGeneratorLUT)
        gen = GNSSSignals.CodeGeneratorLUT(GNSSSignals.CodeReplicaLUT(signal, prn), fs, fc)
        out8g = zeros(Int8, N)
        SUITE["code"]["1ms gen_code! generator"][name] = @benchmarkable gen_code!(
            $out8g, $gen,
        ) evals = 10 samples = 1000
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Threaded multi-channel: build a Vector{CodeReplicaLUT} for 8 PRNs once, then
# time filling 8 per-channel buffers in parallel with `Threads.@threads`. The
# plans are immutable/read-only, so this is thread-safe.
# ─────────────────────────────────────────────────────────────────────────────
if isdefined(GNSSSignals, :CodeReplicaLUT)
    let nch = 8, fs = 5e6Hz, fc = 1023e3Hz
        N = round(Int, ustrip_hz(fs) * 1e-3)
        plans = [GNSSSignals.CodeReplicaLUT(_GPSL1(), prn) for prn in 1:nch]
        outs = [zeros(Int8, N) for _ in 1:nch]
        SUITE["code"]["threaded multi-channel"]["GPSL1CA x8"] = @benchmarkable begin
            Threads.@threads for ch in 1:$nch
                gen_code!($outs[ch], $plans[ch], $fs, $fc, 0.0, 0)
            end
        end evals = 1 samples = 1000
    end
end
