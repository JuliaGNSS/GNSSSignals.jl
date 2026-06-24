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
# 1 ms integration — old gen_code! vs new LUT, head to head (see the loop comment).
# fs picked ≥ fc·subchip_factor per signal; N = round(Int, fs·1e-3).
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

# Per signal, fill N = round(Int, fs·1e-3) samples (one 1 ms integration) at the SAME
# fs/fc three ways, all under one parent `code/1 ms integration/<signal>/…` so the rows
# sit adjacent in the sorted table:
#   "original"      — classic gen_code!(out, signal, prn, fs, fc)              (Int16 out)
#   "LUT generator" — warm CodeGeneratorLUT reused per call (steady state)     (Int8 out, 0 alloc)
#   "LUT one-shot"  — gen_code!(out, plan, fs, fc): rebuilds the DDA each call (Int8 out)
# Plan/generator are built ONCE outside the timed region (a receiver reuses them). Output
# eltype differs by design: Int16 is the original's correlator type, Int8 the LUT's native ±1.
for (name, signal, prn, fs, fc) in _LUT_CASES
    N = round(Int, ustrip_hz(fs) * 1e-3)
    g = SUITE["code"]["1 ms integration"][name]
    out16 = zeros(Int16, N)
    g["original"] = @benchmarkable gen_code!(
        $out16, $signal, $prn, $fs, $fc, $0.0, $0,
    ) evals = 10 samples = 1000

    # Warm continuing generator: DDA init paid once in the constructor (outside the timed
    # region); gen_code!(out, gen) continues the state with no re-init — the fast repeat path.
    if isdefined(GNSSSignals, :CodeGeneratorLUT)
        gen = GNSSSignals.CodeGeneratorLUT(GNSSSignals.CodeReplicaLUT(signal, prn), fs, fc)
        out8g = zeros(Int8, N)
        g["LUT generator"] = @benchmarkable gen_code!($out8g, $gen) evals = 10 samples = 1000
    end

    # One-shot: builds the generator + DDA init on every call (the drop-in gen_code! form).
    if isdefined(GNSSSignals, :CodeReplicaLUT)
        plan = GNSSSignals.CodeReplicaLUT(signal, prn)
        out8 = zeros(Int8, N)
        g["LUT one-shot"] = @benchmarkable gen_code!(
            $out8, $plan, $fs, $fc, $0.0, $0,
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

if isdefined(GNSSSignals, :GalileoE1C)
    SUITE["code"]["code generation"]["GalileoE1C"] = @benchmarkable gen_code!(
        $sampled_code_f32,
        $(GNSSSignals.GalileoE1C()),
        $1,
        $(15e6Hz),
        $(1023e3Hz),
        $0.0,
        $0,
    ) evals = 10 samples = 10000
end

if isdefined(GNSSSignals, :GalileoE1C_BOC11)
    SUITE["code"]["code generation"]["GalileoE1C_BOC11"] = @benchmarkable gen_code!(
        $sampled_code,
        $(GNSSSignals.GalileoE1C_BOC11()),
        $1,
        $(15e6Hz),
        $(1023e3Hz),
        $0.0,
        $0,
    ) evals = 10 samples = 10000
end

# ─────────────────────────────────────────────────────────────────────────────
# OSR sweep — old gen_code! vs new LUT across oversampling ratios, at a small and a
# steady-state buffer. Grouped by signal / osr / size under `code/osr sweep/…` so
# "original" and "LUT" sit adjacent. Same N + fs/fc for both, and the LUT side uses the
# warm generator (0-alloc), matching the original's 0-alloc fill. The LUT is ~flat in OSR
# (one permute/sample); the original's run-fill speeds up with OSR — so the LUT wins most
# at low OSR and the original catches up high (crossover ~OSR 8-16 BPSK, later BOC). Two
# representative signals (BPSK + BOC(1,1)); both at fc = 1.023 MHz.
const _OSR_SIGS = let s = Any[("GPSL1CA", _GPSL1(), 1, 1)]   # (name, signal, prn, subchip_factor P)
    isdefined(GNSSSignals, :GalileoE1B_BOC11) &&
        push!(s, ("GalileoE1B_BOC11", GNSSSignals.GalileoE1B_BOC11(), 1, 2))
    s
end
let fc = 1023e3Hz
    for (name, signal, prn, P) in _OSR_SIGS
        for osr in (2, 8, 32), (slabel, n) in (("4k", 4096), ("64k", 65536))
            osr < P && continue                      # LUT needs fs ≥ fc·P
            fs = osr * fc
            g = SUITE["code"]["osr sweep"][name]["osr$(lpad(osr, 2, '0'))"][slabel]
            o16 = zeros(Int16, n)
            g["original"] =
                @benchmarkable gen_code!($o16, $signal, $prn, $fs, $fc, $0.0, $0) evals = 1 samples = 300
            if isdefined(GNSSSignals, :CodeGeneratorLUT)
                gen = GNSSSignals.CodeGeneratorLUT(GNSSSignals.CodeReplicaLUT(signal, prn), fs, fc)
                o8 = zeros(Int8, n)
                g["LUT"] =
                    @benchmarkable gen_code!($o8, $gen) evals = 1 samples = 300
            end
        end
    end
end
