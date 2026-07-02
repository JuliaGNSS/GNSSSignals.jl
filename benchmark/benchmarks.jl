using BenchmarkTools
using GNSSSignals
using Unitful: Hz, ustrip, @u_str

# Strip a Unitful frequency to a plain Float64 Hz value (for computing N).
ustrip_hz(x) = Float64(ustrip(u"Hz", x))
# Sampling frequency as an MHz string for a benchmark label ("5 MHz", "2.5 MHz").
_mhz(fs) = (m = ustrip_hz(fs) / 1e6; (m == round(m) ? string(Int(m)) : string(m)) * " MHz")

# Use the v2 names when available, fall back to the pre-v2 names so the same
# benchmark script can run against master and this branch. Remove the fallback
# once master is on v2.
const _GPSL1 = isdefined(GNSSSignals, :GPSL1CA) ? GNSSSignals.GPSL1CA : GNSSSignals.GPSL1
const _GPSL5 = isdefined(GNSSSignals, :GPSL5I) ? GNSSSignals.GPSL5I : GNSSSignals.GPSL5

# Detect the LUT era so the buffer eltype can be era-aware across the two revisions benchpkg
# compares: `_EMBEDDED_LUT` marks this branch's embedded per-signal LUT (no plan), whose
# `gen_code!` is Int8-only, vs the PR #69 base's `CodeReplicaLUT` plan.
const _HAS_PLAN     = isdefined(GNSSSignals, :CodeReplicaLUT)
const _EMBEDDED_LUT = isdefined(GNSSSignals, :code_engine) && !_HAS_PLAN
const SUITE = BenchmarkGroup()

# Era-aware output buffer: Int8 on this branch (the embedded-LUT `gen_code!` is Int8-only, so a
# baseline Int16/Float32 buffer would MethodError), else the signal's natural `get_code_type`
# (Int16, or Float32 for the CBOC E1B/E1C) so a pre-embedded base matches. Same SUITE key either
# way, so benchpkg diffs the old generator against the new embedded one head-to-head.
_buf(signal, n) = zeros(_EMBEDDED_LUT ? Int8 : get_code_type(signal), n)

# ─────────────────────────────────────────────────────────────────────────────
# 2000 sample code generation — fixed 2000-sample spot-check at a couple of representative rates. A
# small, stable sanity row for the plain (GPSL1/GPSL5) and CBOC (GalileoE1B, Float32 on a
# pre-embedded base) generators; full per-signal coverage lives in "1 ms code generation" below.
# ─────────────────────────────────────────────────────────────────────────────
for (name, signal, fs) in (("GPSL1", _GPSL1(), 2e6Hz), ("GPSL5", _GPSL5(), 20e6Hz),
                           ("GalileoE1B", GalileoE1B(), 15e6Hz))
    fc = get_code_frequency(signal)
    out = _buf(signal, 2000)
    SUITE["2000 sample code generation"]["$name @ $(_mhz(fs))"] = @benchmarkable gen_code!(
        $out, $signal, $1, $fs, $fc, $0.0, $0,
    ) evals = 10 samples = 10000
end

# ─────────────────────────────────────────────────────────────────────────────
# 1 ms code generation — one gen_code! fill of an N = round(Int, fs·1e-3)-sample epoch per signal, for
# EVERY signal the package defines. `fs` is chosen per (code_frequency, subchip_factor P) family,
# so signals sharing a chip rate + modulation order share a rate and timing differences reflect the
# signal/modulation rather than a different fs; each row carries its fs in the label. fs clears the
# LUT's fs ≥ code_frequency·P requirement with headroom:
#   1.023 MHz  P=1  LOC         →  5 MHz   GPSL1CA
#   1.023 MHz  P=2  BOC(1,1)    →  5 MHz   GPSL1C_D, GalileoE1B_BOC11, GalileoE1C_BOC11
#   1.023 MHz  P=12 CBOC/TMBOC  → 15 MHz   GPSL1C_P, GalileoE1B, GalileoE1C
#   0.5115 MHz P=1  LOC (L2C)   →  2.5 MHz GPSL2CM, GPSL2CL
#   10.23 MHz  P=1  LOC         → 40 MHz   GPSL5I, GPSL5Q, GalileoE5aI, GalileoE5aQ
# Signals absent on the compared revision (e.g. an older base) are skipped via isdefined.
# ─────────────────────────────────────────────────────────────────────────────
const _MS_CASES = let cases = Any[]
    push!(cases, ("GPSL1CA", _GPSL1(), 5e6Hz))
    push!(cases, ("GPSL5I",  _GPSL5(), 40e6Hz))
    for (sym, fs) in (
            (:GPSL1C_D, 5e6Hz), (:GPSL1C_P, 15e6Hz),
            (:GPSL2CM, 2.5e6Hz), (:GPSL2CL, 2.5e6Hz), (:GPSL5Q, 40e6Hz),
            (:GalileoE1B, 15e6Hz), (:GalileoE1B_BOC11, 5e6Hz),
            (:GalileoE1C, 15e6Hz), (:GalileoE1C_BOC11, 5e6Hz),
            (:GalileoE5aI, 40e6Hz), (:GalileoE5aQ, 40e6Hz))
        isdefined(GNSSSignals, sym) &&
            push!(cases, (String(sym), getfield(GNSSSignals, sym)(), fs))
    end
    cases
end
for (name, signal, fs) in _MS_CASES
    fc = get_code_frequency(signal)
    N = round(Int, ustrip_hz(fs) * 1e-3)
    out = _buf(signal, N)
    SUITE["1 ms code generation"]["$name @ $(_mhz(fs))"] = @benchmarkable gen_code!(
        $out, $signal, $1, $fs, $fc, $0.0, $0,
    ) evals = 10 samples = 1000
end

# ─────────────────────────────────────────────────────────────────────────────
# oversampling sweep — gen_code! across oversampling ratios (fs/fc = samples per chip) at a small
# (4k) and steady-state (64k) buffer. Guards the LUT's ratio-dependent path selection: the windowed
# permute (flat in the ratio) at low oversampling, switching to a broadcast run-fill once the baked
# table is oversampled ≳7× (AVX-512) / ≳4× (AVX2). One representative signal (GPS L1 C/A, plain
# BPSK): its oversampling ratio equals the table-oversampling, so 2× exercises the permute path and
# 8×/32× the run-fill across the threshold. 17×/24× are non-power-of-two run-fill ratios: their
# per-chip run length (17/24) is padded to a non-power-of-two store width by `_runfill_pad`, so
# these rows track the run-fill broadcast-store width tuning (the power-of-two ratios never do).
# ─────────────────────────────────────────────────────────────────────────────
let signal = _GPSL1(), fc = 1023e3Hz
    for oversampling in (2, 8, 17, 24, 32), (slabel, n) in (("4k", 4096), ("64k", 65536))
        fs = oversampling * fc
        out = _buf(signal, n)
        # evals = 10 (average out per-call scheduling jitter — the 4k fills run in a few
        # hundred ns, where evals = 1 made the reported minimum swing ±25% on shared runners)
        # and samples = 10000, matching the "2000 sample" group, for a stable minimum.
        SUITE["oversampling sweep"]["$(lpad(oversampling, 2, '0'))x"][slabel] =
            @benchmarkable gen_code!($out, $signal, $1, $fs, $fc, $0.0, $0) evals = 10 samples = 10000
    end
end
