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

# Detect the LUT era so the buffer eltype below can be era-aware across the two revisions
# benchpkg compares: `_EMBEDDED_LUT` marks this branch's embedded per-signal LUT (no plan),
# whose `gen_code!` is Int8-only, vs the PR #69 base's `CodeReplicaLUT` plan.
const _HAS_PLAN     = isdefined(GNSSSignals, :CodeReplicaLUT)
const _EMBEDDED_LUT = isdefined(GNSSSignals, :code_engine) && !_HAS_PLAN
const SUITE = BenchmarkGroup()

# ── gen_code! rows (legacy fixed-size buffers, kept for continuity) ──
# Buffer eltype is era-aware: on this branch the embedded-LUT `gen_code!` is Int8-only, so the
# baseline's Int16/Float32 buffers would MethodError. Using Int8 here keeps the SAME SUITE keys
# so benchpkg can diff the old generator (Int16/Float32 on the base) against the new embedded
# generator (Int8 here) head-to-head across the two revisions.
num_samples = 2000
sampled_code = zeros(_EMBEDDED_LUT ? Int8 : Int16, num_samples)

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

sampled_code_f32 = zeros(_EMBEDDED_LUT ? Int8 : Float32, num_samples)
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

# Per signal, one 1 ms integration (N = round(Int, fs·1e-3) samples) via a single one-shot
# `gen_code!(out, signal, prn, fs, fc)` row under `code/1 ms integration/<signal>`. Same SUITE
# key in every era, so benchpkg diffs the old generator (Int16 buffer on the base) against the
# new embedded-LUT generator (Int8 here) head-to-head. Buffer eltype is era-aware (the embedded
# `gen_code!` is Int8-only; the base's Int16 buffer would MethodError here).
for (name, signal, prn, fs, fc) in _LUT_CASES
    N = round(Int, ustrip_hz(fs) * 1e-3)
    out = zeros(_EMBEDDED_LUT ? Int8 : Int16, N)
    SUITE["code"]["1 ms integration"][name] = @benchmarkable gen_code!(
        $out, $signal, $prn, $fs, $fc, $0.0, $0,
    ) evals = 10 samples = 1000
end

# ── Galileo E1B full CBOC — same `code/1 ms integration/GalileoE1B` group. CBOC needs a Float32
# buffer on the base (the subcarrier amplitudes are irrational); on this branch the embedded LUT
# bakes an Int8 integer approximation (default (19,6)). One one-shot gen_code! row, same SUITE
# key for the cross-revision diff. fs ≥ fc·subchip_factor = 12·1.023 MHz; use 15 MHz.
let signal = GalileoE1B(), prn = 1, fs = 15e6Hz, fc = 1023e3Hz
    N = round(Int, ustrip_hz(fs) * 1e-3)
    out = zeros(_EMBEDDED_LUT ? Int8 : Float32, N)
    SUITE["code"]["1 ms integration"]["GalileoE1B"] = @benchmarkable gen_code!(
        $out, $signal, $prn, $fs, $fc, $0.0, $0,
    ) evals = 10 samples = 1000
end

# ─────────────────────────────────────────────────────────────────────────────
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
# Oversampling sweep — gen_code! across oversampling ratios, at a small and a steady-state
# buffer. "Oversampling ratio" = samples per code chip = fs / fc (2x = sample twice per chip).
# One one-shot gen_code! row per (oversampling, size) under `code/oversampling sweep/…`, same key
# in every era so benchpkg diffs old vs new. The LUT uses the windowed permute (flat in the
# oversampling ratio) at low oversampling and switches to a broadcast run-fill once the baked
# table is oversampled ≳7× (AVX-512) / ≳4× (AVX2). One representative signal (GPS L1 C/A, plain
# BPSK): its oversampling ratio equals the table-oversampling, so 2× exercises the permute path
# and 8×/32× the run-fill it switches to across the threshold.
let name = "GPSL1CA", signal = _GPSL1(), prn = 1, fc = 1023e3Hz
    for oversampling in (2, 8, 32), (slabel, n) in (("4k", 4096), ("64k", 65536))
        fs = oversampling * fc
        out = zeros(_EMBEDDED_LUT ? Int8 : Int16, n)   # Int8 embedded gen_code! on this branch
        SUITE["code"]["oversampling sweep"][name]["$(lpad(oversampling, 2, '0'))x"][slabel] =
            @benchmarkable gen_code!($out, $signal, $prn, $fs, $fc, $0.0, $0) evals = 1 samples = 300
    end
end
