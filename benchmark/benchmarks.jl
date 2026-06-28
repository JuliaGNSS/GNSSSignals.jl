using BenchmarkTools
using GNSSSignals
using SIMD: Vec, vload, vstore, shufflevector
import SinCosLUT
using SinCosLUT: SinCosTable, generate_carrier!, cycles_per_sample,
                 carrier_engine, carrier_state, carrier_lookup, carrier_advance, carrier_width
using Unitful: Hz, ustrip, @u_str

# Polyester is optional: it only backs the parallel multi-channel correlation rows below.
# Guard the import so the rest of the suite still loads (and diffs against a baseline) on an
# environment without it. Added to the benchpkg env via the workflow's `--add` list.
const _HAS_POLYESTER = try
    @eval import Polyester
    true
catch
    false
end

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

# ── Galileo E1B full CBOC, head to head — same `code/1 ms integration/GalileoE1B/…` group
# as the rows above, but special-cased because CBOC differs from the ±1 cases:
#   • the ORIGINAL `gen_code!` needs a Float32 buffer (the CBOC subcarrier amplitudes are
#     irrational), whereas the LUT bakes an Int8 integer approximation (default (19,6));
#   • the LUT plan only supports CBOC on this branch, so the LUT rows are probed and skipped
#     on a baseline (the PR base #69) that errors on CBOC — the `original` row still compares.
# fs ≥ fc·subchip_factor = 12·1.023 MHz; use 15 MHz (matches the legacy E1B row).
const _LUT_CBOC_OK = isdefined(GNSSSignals, :CodeReplicaLUT) && try
    GNSSSignals.CodeReplicaLUT(GalileoE1B(), 1)   # CBOC supported here, errors on the baseline
    true
catch
    false
end
let signal = GalileoE1B(), prn = 1, fs = 15e6Hz, fc = 1023e3Hz
    N = round(Int, ustrip_hz(fs) * 1e-3)
    g = SUITE["code"]["1 ms integration"]["GalileoE1B"]
    out_f32 = zeros(Float32, N)
    g["original"] = @benchmarkable gen_code!(
        $out_f32, $signal, $prn, $fs, $fc, $0.0, $0,
    ) evals = 10 samples = 1000
    if _LUT_CBOC_OK
        gen = GNSSSignals.CodeGeneratorLUT(GNSSSignals.CodeReplicaLUT(signal, prn), fs, fc)
        out8g = zeros(Int8, N)
        g["LUT generator"] = @benchmarkable gen_code!($out8g, $gen) evals = 10 samples = 1000

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

# ─────────────────────────────────────────────────────────────────────────────
# Oversampling sweep — old gen_code! vs new LUT across oversampling ratios, at a small
# and a steady-state buffer. "Oversampling ratio" = samples per code chip = fs / fc;
# the level is labelled as a multiplier (2x = sample twice per chip). Grouped by signal /
# oversampling / size under `code/oversampling sweep/…` so "original" and "LUT" sit
# adjacent. Same N + fs/fc for both, and the LUT side uses the warm generator (0-alloc),
# matching the original's 0-alloc fill. The LUT uses the windowed permute (flat in the
# oversampling ratio) at low oversampling and switches to a broadcast run-fill once the
# baked table is oversampled ≳7× (AVX-512) / ≳4× (AVX2), so it tracks the original's
# run-fill speed-up at high oversampling instead of staying flat. It still wins outright at
# low oversampling and for modulated signals (baked subcarrier).
# Two representative signals (BPSK + BOC(1,1)); both at fc = 1.023 MHz.
const _SWEEP_SIGS = let s = Any[("GPSL1CA", _GPSL1(), 1, 1)]   # (name, signal, prn, subchip_factor P)
    isdefined(GNSSSignals, :GalileoE1B_BOC11) &&
        push!(s, ("GalileoE1B_BOC11", GNSSSignals.GalileoE1B_BOC11(), 1, 2))
    s
end
let fc = 1023e3Hz
    for (name, signal, prn, P) in _SWEEP_SIGS
        for oversampling in (2, 8, 32), (slabel, n) in (("4k", 4096), ("64k", 65536))
            oversampling < P && continue             # LUT needs fs ≥ fc·P
            fs = oversampling * fc
            g = SUITE["code"]["oversampling sweep"][name]["$(lpad(oversampling, 2, '0'))x"][slabel]
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

# ── Run-fill threshold crossover sweep (validates the permute↔run-fill selection) ──────
# The coarse sweep above (2/8/32×) straddles but never lands in the crossover region, so a
# threshold retune produces no visible diff there. This dense sweep hits the exact (m, N)
# cells where the per-backend run-fill threshold flips which kernel runs — m around the
# measured crossovers (3 for AVX2/NEON, 5–7 for AVX-512) at short and long fills (the
# short ones exercise the AVX-512 N-aware step). GPS L1 C/A has subchip_factor P = 1, so
# m = oversampling exactly. LUT only — the threshold is internal to the LUT path.
if isdefined(GNSSSignals, :CodeGeneratorLUT)
    let fc = 1023e3Hz, signal = _GPSL1(), prn = 1
        for m in (3, 5, 6, 7), (slabel, n) in (("512", 512), ("2k", 2048), ("64k", 65536))
            fs = m * fc
            gen = GNSSSignals.CodeGeneratorLUT(GNSSSignals.CodeReplicaLUT(signal, prn), fs, fc)
            o8 = zeros(Int8, n)
            SUITE["code"]["runfill crossover"]["$(lpad(m, 2, '0'))x"][slabel] =
                @benchmarkable gen_code!($o8, $gen) evals = 1 samples = 500
        end
    end
end

# ── NEON run-fill crossover probe (Apple Silicon only) ────────────────────────────────
# The NEON crossover can't be measured on the x86 CI box (x86 can't emit NEON), so probe the
# two kernels directly on the macos-14 runner: time the windowed permute vs the broadcast
# run-fill at several m, so the NEON crossover can be read off the table and the threshold
# (currently 3) confirmed/retuned. Gated to NEON hosts — forcing the Neon kernel on x86 would
# fail to compile — so these keys appear only in the Apple-Silicon results. Base and head call
# identical kernels here, so both columns report the same NEON times (we want the absolutes).
if isdefined(GNSSSignals, :CodeGeneratorLUT) &&
   GNSSSignals.CodeLUT.default_backend() isa GNSSSignals.CodeLUT.Neon
    let SD = GNSSSignals.CodeLUT._STEP_DEN, be = GNSSSignals.CodeLUT.Neon()
        tbl = GNSSSignals.CodeReplicaLUT(_GPSL1(), 1).mc.table
        for m in (2, 3, 4, 5), (slabel, n) in (("512", 512), ("8k", 8192))
            sn = SD ÷ m
            o = zeros(Int8, n)
            grp = SUITE["code"]["neon crossover probe"]["$(lpad(m, 2, '0'))x"][slabel]
            grp["permute"] =
                @benchmarkable GNSSSignals.CodeLUT._generate!($o, $tbl, $sn, $SD, 0, $be) evals = 1 samples = 500
            grp["runfill"] =
                @benchmarkable GNSSSignals.CodeLUT._generate_runfill!($o, $tbl, $sn, $SD, 0) evals = 1 samples = 500
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Correlation signal path — FUSED vs UNFUSED, full Early / Prompt / Late (E/P/L).
#
# Models one tracking update with the three standard correlators. Front-end words:
#   measurement — COMPLEX Int16, re/im interleaved [reₙ, imₙ, …] (12-bit ADC, |m| ≲ 2048),
#                 a single Vector{Int16} of length 2·N.
#   carrier     — Int8, small amplitude (local carrier cos + j·sin); SHARED by E/P/L.
#   code        — Int8 ±1 replica (CBOC/BOC can reach ±2); Early/Late are the Prompt code
#                 shifted by ∓D samples (½-chip spacing ⇒ D = round(0.5·fs/fc)).
# Per-correlator output is the complex sum of measurement × codeₓ × conj(cos + j·sin):
#   Iₓ = Σ codeₓ·(mᵣ·cos + mᵢ·sin),   Qₓ = Σ codeₓ·(mᵢ·cos − mᵣ·sin),   x ∈ {E, P, L}.
#
# SHARED CARRIER WIPE-OFF. The carrier is identical for all three correlators, so the
# per-sample wipe DI = mᵣ·cos + mᵢ·sin, DQ = mᵢ·cos − mᵣ·sin is computed ONCE and each
# correlator reduces to a cheap Σ codeₓ·DI / Σ codeₓ·DQ — ~1.8× faster than three
# independent correlations. The wipe is a per-sample complex multiply: deinterleave mᵣ/mᵢ and
# multiply by the carrier. The arithmetic width is chosen at compile time (see `choose_carrier`):
#   • Int16 fast path — when 2·max|meas|·amplitude ≤ typemax(Int16) (true for the 12-bit ADC at the
#     auto-chosen amplitude), DI/DQ stay Int16 (`vpmullw`, half the SIMD register width). The code
#     step Σ codeₓ·DI is then a fused widening multiply-add: `vpmaddwd` (Int16×Int16→Int32 pairwise)
#     into half-width accumulators on AVX2/AVX-512, or `SMLAL`/`vmlal` (per-lane widening MAC) into
#     full-width accumulators on NEON. NOTE: this fused MAC pays off for the *reduced-sum* code step,
#     but NOT for the per-sample wipe (`vpmaddwd`'s tiling shuffles + 4× carrier bandwidth lose there
#     — measured); so the wipe stays a straight deinterleave-multiply. Net ~2× on AVX-512, ~1.2× on
#     AVX2; NEON expected to gain from the narrower wipe + SMLAL (verified on the macos-14 runner).
#   • Int32 fallback — exact for any amplitude (per-sample products can reach 2048·127·2 > Int16), and
#     the only path on the scalar `Portable` backend. Deinterleave, widen to Int32, `vpmulld`.
# Both paths are bit-identical (the Int16 path provably cannot overflow when selected). The
# amplitude is auto-maximised within the Int16-safe bound, so a coarser ADC keeps a finer carrier.
#
# EARLY/LATE = PROMPT SHIFTED BY D SAMPLES. One Prompt code stream is built; Early/Late index
# it ∓D samples away. This generates the code ONCE (N+2D samples) instead of three times and
# — crucially — it FUSES: the fused kernel looks the Prompt code up once per W-wide chunk and
# derives Early/Late by register lane-shifts across the carried neighbouring chunks (a
# `shufflevector`). With width W the shift spans ⌈D/W⌉ neighbour chunks, so the kernel carries
# a small sliding pipeline of P = ⌈D/W⌉+⌊D/W⌋+2 code chunks (registers only; one lookup + one
# DDA advance + two shuffles per chunk, independent of the correlator count). The alternative
# (generate Early/Late independently — three lookups + advances) is ~20–25 % slower and, when
# fused, can't even express a sub-chip spacing (the value-based code engine takes an
# integer-chip phase offset), so the sample-shift is the natural fused form.
#
#   UNFUSED — materialise ONE extended Prompt code buffer (samples −D … N−1+D) + the shared
#             carrier into Int8 buffers, then stream them through the shared-wipe loop with
#             Early/Prompt/Late as offset reads (ext[n] / ext[n+D] / ext[n+2D]).
#   FUSED   — never materialise: drive the value-based Prompt code engine + carrier engine
#             straight into the loop; Early/Late are lane-shifts of the looked-up Prompt code.
#
# Code before sample 0 is taken as zero (a D-sample edge on Early only; applied identically to
# fused / unfused / a scalar reference, so all three stay bit-for-bit equal). The < W
# measurement tail rides on the zero-padding of `meas` to a W multiple (padding · anything = 0).
# Benchmarked at GPS L1 C/A 1 ms for two rates: 5 MHz (N=5000, D=2) and 40 MHz (N=40000, D=20;
# there D > W on NEON, exercising the multi-chunk pipeline). Plans / tables / engines / buffers
# are built once in setup (a receiver reuses them); the timed region runs the correlation loop.
# ─────────────────────────────────────────────────────────────────────────────
# SIMD width of the Int8 carrier/code permute backend on this host (compile-time const).
# Derived from SinCosLUT (always available); the code LUT uses the same Int8 width per
# backend, so a single `W` drives both fused iterators. `default_backend` and the backend
# singletons are taken from SinCosLUT so this file still loads on a baseline without the
# code-LUT API (the correlation benchmarks below are then simply skipped).
const _CORR_W = let be = SinCosLUT.default_backend(Int8, 64)
    be isa SinCosLUT.AVX512 ? 64 :
    be isa SinCosLUT.AVX2 ? 32 :
    be isa SinCosLUT.Neon ? 16 : 1
end
# Strip-mining block length (samples) for the "hybrid-blocked" kernel below. A multiple of every
# backend `W`; sized so its three L1-resident scratch buffers (≈ 3·_CORR_BLK bytes per channel)
# stay in L1 and are reused across blocks rather than spilling to the LLC.
const _CORR_BLK = 8192

# ── Int16 carrier-wipe fast path ────────────────────────────────────────────
# The shared wipe DI = mᵣ·cos + mᵢ·sin / DQ = mᵢ·cos − mᵣ·sin is the dominant per-sample cost.
# It can run entirely in Int16 (half the SIMD register width of the Int32 form) *iff* it cannot
# overflow: |DI| ≤ 2·max|meas|·amplitude ≤ typemax(Int16). The code accumulate Σ codeₓ·DI then uses
# a fused widening multiply-add: `vpmaddwd` on AVX2/AVX-512 (pairwise, half-width accumulators) or
# `SMLAL`/`vmlal` on AArch64 NEON (per-lane widening MAC, full-width accumulators — emitted by LLVM
# from the widening pattern). So the fast path runs on any SIMD backend (W>1); only the truly scalar
# `Portable` (W=1) backend forgoes it.
const _HAS_INT16 = _CORR_W in (16, 32, 64)     # NEON (W=16) / AVX2 (W=32) / AVX-512 (W=64)

# Pick (carrier amplitude, wipe intermediate type) for a declared maximum |measurement component|.
# REQUIRED input, NO default: under-declaring `max_meas` would silently overflow the Int16 wipe and
# corrupt results, while over-declaring only forgoes the speedup. We take the LARGEST Int16-safe
# amplitude (carrier rounding error is ±0.5 regardless of amplitude, so relative error ≈ 0.5/amp —
# bigger amp = finer carrier), capped at the Int8 storage limit. If no amplitude ≥1 is safe (very
# large `max_meas`) or the backend is scalar, fall back to the exact Int32 wipe at full Int8
# amplitude (Int32 has the headroom, so best fidelity there).
function choose_carrier(max_meas::Integer)
    _HAS_INT16 || return (Int(typemax(Int8)), Int32)
    a = Int(typemax(Int16)) ÷ (2 * Int(max_meas))           # ⌊32767 / (2·max_meas)⌋, may be 0
    a >= 1 ? (min(a, Int(typemax(Int8))), Int16) : (Int(typemax(Int8)), Int32)
end

# The benchmark measurement is a 12-bit ADC sample, so |meas| ≤ 2^11. One source of truth for the
# carrier amplitude + wipe type, shared by every kernel so all four stay mutually bit-identical.
const _ADC_BITS = 12
const _MAX_MEAS = 1 << (_ADC_BITS - 1)                       # 2048
const _CARRIER_AMP, _WIPE_TI = choose_carrier(_MAX_MEAS)     # (7, Int16) on AVX2/AVX-512

@inline _wide32(v::Vec{W,Int8}) where {W} = convert(Vec{W,Int32}, v)
@inline _wide32(v::Vec{W,Int16}) where {W} = convert(Vec{W,Int32}, v)
@inline _wide16(v::Vec{W,Int8}) where {W} = convert(Vec{W,Int16}, v)

# Deinterleave [re, im, …] → real (even lanes) / imag (odd lanes), Vec{2W,Int16} → Vec{W,Int16}.
@inline _re(v::Vec{M,Int16}) where {M} = shufflevector(v, Val(ntuple(i -> 2 * (i - 1), Val(M ÷ 2))))
@inline _im(v::Vec{M,Int16}) where {M} = shufflevector(v, Val(ntuple(i -> 2 * (i - 1) + 1, Val(M ÷ 2))))

# ── vpmaddwd: Int16×Int16 → Int32 pairwise-add — the dot-product primitive for the code accumulate
# Σ codeₓ·DI. Native 512-/256-bit intrinsics tiled to width W; x86-only (the Int16 path is gated on
# `_HAS_MADDWD`, so these are never instantiated on other backends). Output is Vec{W÷2,Int32}
# (adjacent samples pre-summed); the final `sum` over the accumulator is identical to a per-lane
# Int32 reduction (integer addition is associative), so the result is bit-exact with the Int32 path.
@static if Sys.ARCH in (:x86_64, :i686)
    @inline _madd_tile(a::Vec{M,Int16}, ::Val{o}, ::Val{t}) where {M,o,t} =
        shufflevector(a, Val(ntuple(i -> i - 1 + o, Val(t))))
    @inline function _madd512(a::Vec{32,Int16}, b::Vec{32,Int16})
        Vec(Base.llvmcall(("""
            declare <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16>, <32 x i16>)
            define <16 x i32> @entry(<32 x i16> %a, <32 x i16> %b) #0 {
              %r = call <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16> %a, <32 x i16> %b)
              ret <16 x i32> %r }
            attributes #0 = { alwaysinline }""", "entry"), NTuple{16,Base.VecElement{Int32}},
            Tuple{NTuple{32,Base.VecElement{Int16}}, NTuple{32,Base.VecElement{Int16}}}, a.data, b.data))
    end
    @inline function _madd256(a::Vec{16,Int16}, b::Vec{16,Int16})
        Vec(Base.llvmcall(("""
            declare <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16>, <16 x i16>)
            define <8 x i32> @entry(<16 x i16> %a, <16 x i16> %b) #0 {
              %r = call <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16> %a, <16 x i16> %b)
              ret <8 x i32> %r }
            attributes #0 = { alwaysinline }""", "entry"), NTuple{8,Base.VecElement{Int32}},
            Tuple{NTuple{16,Base.VecElement{Int16}}, NTuple{16,Base.VecElement{Int16}}}, a.data, b.data))
    end
    @inline function _maddacc(a::Vec{64,Int16}, b::Vec{64,Int16})        # AVX-512: 2×512-bit tiles
        lo = _madd512(_madd_tile(a, Val(0), Val(32)),  _madd_tile(b, Val(0), Val(32)))
        hi = _madd512(_madd_tile(a, Val(32), Val(32)), _madd_tile(b, Val(32), Val(32)))
        shufflevector(lo, hi, Val(ntuple(i -> i - 1, Val(32))))
    end
    @inline function _maddacc(a::Vec{32,Int16}, b::Vec{32,Int16})        # AVX2: 2×256-bit tiles
        lo = _madd256(_madd_tile(a, Val(0), Val(16)),  _madd_tile(b, Val(0), Val(16)))
        hi = _madd256(_madd_tile(a, Val(16), Val(16)), _madd_tile(b, Val(16), Val(16)))
        shufflevector(lo, hi, Val(ntuple(i -> i - 1, Val(16))))
    end
end

# Per-sample shared carrier wipe → (DI, DQ). Two paths, selected at compile time by the wipe type:
#   Int32 — deinterleave mᵣ/mᵢ, widen to Int32, multiply. Exact for any amplitude; the fallback.
#   Int16 — same but stays in Int16 (vpmullw), valid only when 2·max|meas|·amplitude ≤ typemax(Int16)
#           (see `choose_carrier`). Half the register width on AVX, and feeds the vpmaddwd accumulate.
@inline function _wipe(::Type{Int32}, v::Vec{M,Int16}, sinv::Vec{W,Int8}, cosv::Vec{W,Int8}) where {M,W}
    mr = _wide32(_re(v)); mi = _wide32(_im(v)); cs = _wide32(cosv); sn = _wide32(sinv)
    (mr * cs + mi * sn, mi * cs - mr * sn)
end
@inline function _wipe(::Type{Int16}, v::Vec{M,Int16}, sinv::Vec{W,Int8}, cosv::Vec{W,Int8}) where {M,W}
    mr = _re(v); mi = _im(v); cs = _wide16(cosv); sn = _wide16(sinv)
    (mr * cs + mi * sn, mi * cs - mr * sn)
end

# Six E/P/L accumulators as a tuple (IE,QE,IP,QP,IL,QL); `_acc_zeros`/`_acc_sum` bracket the loop
# and `_acc` folds one chunk. The accumulate dispatches on the DI/DQ element type the wipe produced,
# and (for the Int16 path) on the backend width:
#   Int32 DI                → widen code to Int32, vpmulld into full-width (Vec{W,Int32}) accumulators.
#   Int16 DI, W∈{32,64} AVX  → vpmaddwd (pairwise) into HALF-width (Vec{W÷2,Int32}) accumulators.
#   Int16 DI, W=16   NEON    → per-lane widening MAC → SMLAL, into FULL-width (Vec{16,Int32}) accs.
@inline _acc_zeros(::Type{Int32}, ::Val{W}) where {W} = ntuple(_ -> zero(Vec{W,Int32}), Val(6))
@inline _acc_zeros(::Type{Int16}, ::Val{16}) = ntuple(_ -> zero(Vec{16,Int32}), Val(6))      # NEON: per-lane
@inline _acc_zeros(::Type{Int16}, ::Val{W}) where {W} = ntuple(_ -> zero(Vec{W ÷ 2,Int32}), Val(6))  # AVX: pairwise
@inline _acc_sum(a) = ((sum(a[1]), sum(a[2])), (sum(a[3]), sum(a[4])), (sum(a[5]), sum(a[6])))
@inline function _acc(a, cE::Vec{W,Int8}, cP::Vec{W,Int8}, cL::Vec{W,Int8}, DI::Vec{W,Int32}, DQ::Vec{W,Int32}) where {W}
    e = _wide32(cE); p = _wide32(cP); l = _wide32(cL)
    (a[1] + e * DI, a[2] + e * DQ, a[3] + p * DI, a[4] + p * DQ, a[5] + l * DI, a[6] + l * DQ)
end
# NEON Int16 path (W=16): widen code & DI to Int32 from Int16 so LLVM emits SMLAL (sext-i16 operands).
@inline function _acc(a, cE::Vec{16,Int8}, cP::Vec{16,Int8}, cL::Vec{16,Int8}, DI::Vec{16,Int16}, DQ::Vec{16,Int16})
    e = _wide32(_wide16(cE)); p = _wide32(_wide16(cP)); l = _wide32(_wide16(cL))
    di = _wide32(DI); dq = _wide32(DQ)
    (a[1] + e * di, a[2] + e * dq, a[3] + p * di, a[4] + p * dq, a[5] + l * di, a[6] + l * dq)
end
# AVX Int16 path (W∈{32,64}): vpmaddwd pairwise into half-width accumulators.
@inline function _acc(a, cE::Vec{W,Int8}, cP::Vec{W,Int8}, cL::Vec{W,Int8}, DI::Vec{W,Int16}, DQ::Vec{W,Int16}) where {W}
    e = _wide16(cE); p = _wide16(cP); l = _wide16(cL)
    (a[1] + _maddacc(e, DI), a[2] + _maddacc(e, DQ), a[3] + _maddacc(p, DI),
     a[4] + _maddacc(p, DQ), a[5] + _maddacc(l, DI), a[6] + _maddacc(l, DQ))
end

# ── FUSED E/P/L: one Prompt code lookup per chunk; Early/Late are lane-shifts across a sliding
# pipeline of P = ⌈D/W⌉+⌊D/W⌋+2 carried code chunks. `@generated` so the pipeline length, the
# per-correlator chunk indices and the shuffle masks fold to literals (fully unrolled, 0-alloc)
# for any (W, D). The Prompt code is a single value-based stream; nothing is materialised.
@generated function correlate_epl_fused(meas::Vector{Int16}, ceng1, reng, ::Val{W}, ::Val{D}, ::Val{TI}) where {W,D,TI}
    g = cld(D, W); fl = fld(D, W); P = g + fl + 2
    re = g * W - D; rl = D - fl * W
    csym = Vector{Symbol}(undef, P); for i in 1:P; csym[i] = Symbol(:c_, i); end
    em = Vector{Int}(undef, W); for i in 1:W; em[i] = (i - 1) + re; end; emask = Tuple(em)
    lm = Vector{Int}(undef, W); for i in 1:W; lm[i] = (i - 1) + rl; end; lmask = Tuple(lm)
    prime = Expr[]                                           # c_1..c_g = 0 (samples < 0); c_{g+1}.. = lookups
    for i in 1:g; push!(prime, :($(csym[i]) = z)); end
    for i in g+1:P
        push!(prime, :($(csym[i]) = code_lookup(ceng1, cs))); push!(prime, :(cs = code_advance(ceng1, cs)))
    end
    shift = Expr[]                                           # slide the pipeline left, push the next lookup
    for i in 1:P-1; push!(shift, :($(csym[i]) = $(csym[i+1]))); end
    push!(shift, :($(csym[P]) = code_lookup(ceng1, cs))); push!(shift, :(cs = code_advance(ceng1, cs)))
    quote
        acc = _acc_zeros($TI, Val($W))
        cs = code_state(ceng1, 0); rs = carrier_state(reng, 0); z = zero(Vec{$W,Int8})
        $(prime...)
        n = 0; lim = length(meas) ÷ 2
        @inbounds while n < lim
            sinv, cosv = carrier_lookup(reng, rs); rs = carrier_advance(reng, rs, 1)
            DI, DQ = _wipe($TI, vload(Vec{$(2W),Int16}, meas, 2n + 1), sinv, cosv)
            early = shufflevector($(csym[1]), $(csym[2]), Val($emask))
            late = shufflevector($(csym[g+fl+1]), $(csym[g+fl+2]), Val($lmask))
            acc = _acc(acc, early, $(csym[g+1]), late, DI, DQ)
            $(shift...)
            n += $W
        end
        _acc_sum(acc)
    end
end

# ── UNFUSED E/P/L: materialise the extended Prompt code (samples −D … N−1+D, leading D zero)
# + the shared carrier, then run the shared-wipe loop with Early/Prompt/Late as offset reads.
# `ext[k]` holds the code at sample k−1−D, so at output sample n Early/Prompt/Late read
# ext[n+1] / ext[n+D+1] / ext[n+2D+1]. Generator is warm (continuing); timed region is 0-alloc.
function correlate_epl_unfused!(ext::Vector{Int8}, csin::Vector{Int8}, ccos::Vector{Int8},
                                meas::Vector{Int16}, code_gen, tbl, fs, freq, N::Int,
                                ::Val{W}, ::Val{D}, ::Val{TI}) where {W,D,TI}
    gen_code!(view(ext, (D + 1):(D + N + D)), code_gen)     # samples 0 … N+D−1 (ext[1:D] stay 0)
    generate_carrier!(view(csin, 1:N), view(ccos, 1:N), tbl; frequency = freq, sampling_frequency = fs)
    acc = _acc_zeros(TI, Val(W))
    n = 0
    @inbounds while n < length(csin)
        DI, DQ = _wipe(TI, vload(Vec{2W,Int16}, meas, 2n + 1),
                       vload(Vec{W,Int8}, csin, n + 1), vload(Vec{W,Int8}, ccos, n + 1))
        acc = _acc(acc, vload(Vec{W,Int8}, ext, n + 1), vload(Vec{W,Int8}, ext, n + D + 1),
                   vload(Vec{W,Int8}, ext, n + 2D + 1), DI, DQ)
        n += W
    end
    _acc_sum(acc)
end

# ── HYBRID E/P/L: materialise ONLY the code (fast run-fill `gen_code!`, same `ext` layout as
# unfused), but keep the CARRIER in-register via the value-based engine and FUSE the wipe +
# E/P/L accumulate. Code is read at Early/Prompt/Late offsets from `ext` (cheap loads); the
# carrier never touches memory. So this writes/reads ONE scratch buffer per channel (≈ N Int8)
# instead of the unfused's three (code + sin + cos), while still getting run-fill's fast code
# generation that the fully-fused path lacks. Arithmetic is bit-identical to both kernels.
function correlate_epl_hybrid!(ext::Vector{Int8}, meas::Vector{Int16}, code_gen, reng, N::Int,
                               ::Val{W}, ::Val{D}, ::Val{TI}) where {W,D,TI}
    gen_code!(view(ext, (D + 1):(D + N + D)), code_gen)     # samples 0 … N+D−1 (ext[1:D] stay 0)
    acc = _acc_zeros(TI, Val(W))
    rs = carrier_state(reng, 0)
    n = 0; lim = length(meas) ÷ 2
    @inbounds while n < lim
        sinv, cosv = carrier_lookup(reng, rs); rs = carrier_advance(reng, rs, 1)
        DI, DQ = _wipe(TI, vload(Vec{2W,Int16}, meas, 2n + 1), sinv, cosv)
        acc = _acc(acc, vload(Vec{W,Int8}, ext, n + 1), vload(Vec{W,Int8}, ext, n + D + 1),
                   vload(Vec{W,Int8}, ext, n + 2D + 1), DI, DQ)
        n += W
    end
    _acc_sum(acc)
end

# ── HYBRID-BLOCKED E/P/L: the hybrid's strategy taken to its conclusion — strip-mine the whole
# correlation into BLK-sample blocks and regenerate BOTH the code and the carrier into SMALL,
# L1-resident scratch that is REUSED across blocks. So, unlike the plain hybrid (which keeps a
# full ≈N-sample code buffer per channel), this materialises NOTHING of length N: its working set
# is just `3·BLK` scratch bytes per channel (code block + sin + cos), independent of N.
#
# Per block it runs two tight loops — fill (code via run-fill `gen_code!`, carrier via the value
# engine in a 4-way-unrolled fill, both shuffle-port bound) then correlate (cheap `ext`/carrier
# loads + the multiply-port-bound wipe + E/P/L accumulate) — so each loop hits peak port
# utilisation the way the unfused's two full-length passes do, but without ever touching the LLC.
# That combination is the point: it is fused-class in LLC footprint (so it keeps scaling with
# cores once the working set spills, where even the plain hybrid's code buffer eventually loses to
# the fused's zero-buffer kernel) yet hybrid-class in per-sample cost (cheap loads, not per-chunk
# value-engine permutes). Arithmetic is bit-identical to all three other kernels.
#
# `extb` holds samples (b−D) … (b+len−1+D) for the current block at absolute start `b` in entries
# 1 … len+2D (Early/Prompt/Late at output n read entries n+1 / n+D+1 / n+2D+1, as in the unfused).
# The leading D entries of block 0 are the sample-<0 zero edge; every later block carries the
# 2D-sample overlap from the previous block's tail (the run-fill generator is continuing, so the
# code stays exact across the block seam). `BLK` is a multiple of `W`, so every block length is too.

# Fill `len` carrier samples starting at absolute sample `start` into csb/ccb, 4-way unrolled so
# the permute lookups pipeline (the same scheme `generate_carrier!` uses internally, but driven by
# the value engine at an exact absolute start, so it is bit-identical to the fused carrier).
@inline function _epl_fill_carrier!(csb::Vector{Int8}, ccb::Vector{Int8}, reng,
                                    start::Int, len::Int, ::Val{W}) where {W}
    s0 = carrier_state(reng, start);      s1 = carrier_state(reng, start + W)
    s2 = carrier_state(reng, start + 2W); s3 = carrier_state(reng, start + 3W)
    n = 0
    @inbounds while n + 4W <= len
        a0, b0 = carrier_lookup(reng, s0); a1, b1 = carrier_lookup(reng, s1)
        a2, b2 = carrier_lookup(reng, s2); a3, b3 = carrier_lookup(reng, s3)
        vstore(a0, csb, n + 1);      vstore(b0, ccb, n + 1)
        vstore(a1, csb, n + W + 1);  vstore(b1, ccb, n + W + 1)
        vstore(a2, csb, n + 2W + 1); vstore(b2, ccb, n + 2W + 1)
        vstore(a3, csb, n + 3W + 1); vstore(b3, ccb, n + 3W + 1)
        s0 = carrier_advance(reng, s0, 4); s1 = carrier_advance(reng, s1, 4)
        s2 = carrier_advance(reng, s2, 4); s3 = carrier_advance(reng, s3, 4)
        n += 4W
    end
    @inbounds while n < len                                  # < 4W tail, one chunk at a time
        a0, b0 = carrier_lookup(reng, s0); vstore(a0, csb, n + 1); vstore(b0, ccb, n + 1)
        s0 = carrier_advance(reng, s0, 1)
        n += W
    end
    nothing
end

function correlate_epl_hybrid_blocked!(extb::Vector{Int8}, csb::Vector{Int8}, ccb::Vector{Int8},
                                       meas::Vector{Int16}, code_gen, reng,
                                       ::Val{W}, ::Val{D}, ::Val{BLK}, ::Val{TI}) where {W,D,BLK,TI}
    acc = _acc_zeros(TI, Val(W))
    lim = length(meas) ÷ 2
    b = 0; prevlen = 0
    @inbounds while b < lim
        len = min(BLK, lim - b)
        if prevlen == 0                                      # block 0: leading D zeros, fill 0…len+D−1
            for i in 1:D; extb[i] = Int8(0); end
            gen_code!(view(extb, (D + 1):(D + len + D)), code_gen)
        else                                                 # carry 2D overlap from the prev tail, fill len
            for i in 1:2D; extb[i] = extb[prevlen + i]; end
            gen_code!(view(extb, (2D + 1):(2D + len)), code_gen)
        end
        _epl_fill_carrier!(csb, ccb, reng, b, len, Val(W))
        n = 0
        @inbounds while n < len
            DI, DQ = _wipe(TI, vload(Vec{2W,Int16}, meas, 2(b + n) + 1),
                           vload(Vec{W,Int8}, csb, n + 1), vload(Vec{W,Int8}, ccb, n + 1))
            acc = _acc(acc, vload(Vec{W,Int8}, extb, n + 1), vload(Vec{W,Int8}, extb, n + D + 1),
                       vload(Vec{W,Int8}, extb, n + 2D + 1), DI, DQ)
            n += W
        end
        b += len; prevlen = len
    end
    _acc_sum(acc)
end

# Prompt code engine (a single stream — Early/Late are derived by lane-shift) + carrier engine.
# Built once and reused. `code_engine`'s parameterless `default_backend()` const-folds to the
# host backend, so the engine types are concrete and the fused kernel stays 0-alloc.
@inline _epl_engines(plan, tbl, fs, fc, freq) =
    (GNSSSignals.code_engine(plan, fs, fc, Val(1)),
     carrier_engine(tbl; frequency = freq, sampling_frequency = fs))

# Register the E/P/L benchmarks only where the value-based code engine exists (skipped on a
# baseline without the code-LUT API, so AirspeedVelocity can still diff against it). Two rates:
# 5 MHz (N=5000, D=2; D<W on every backend) and 40 MHz (N=40000, D=20; D>W on NEON).
if isdefined(GNSSSignals, :code_engine)
    for (label, fs, fc) in (
        ("GPSL1CA E/P/L 5MHz 1ms", 5e6, 1.023e6),
        ("GPSL1CA E/P/L 40MHz 1ms", 40e6, 1.023e6),
    )
        let
            freq = 1234.0                              # residual carrier / Doppler to wipe off (Hz)
            N = round(Int, fs * 1e-3)                  # 1 ms
            D = round(Int, 0.5 * fs / fc)              # ½-chip Early/Late spacing, in samples
            # Carrier amplitude + wipe type auto-chosen from the 12-bit ADC range (Int16 fast path).
            plan = GNSSSignals.CodeReplicaLUT(_GPSL1(), 1)
            code_gen = GNSSSignals.CodeGeneratorLUT(plan, fs, fc)   # warm generator for the unfused fill
            tbl = SinCosTable(Int8; steps = 64, amplitude = _CARRIER_AMP)
            ceng1, reng = _epl_engines(plan, tbl, fs, fc, freq)

            # Zero-pad to a whole number of W-wide chunks (Npad samples). meas is complex with
            # re/im interleaved → 2·Npad Int16 words. ext spans samples −D … (Npad−1)+D.
            Npad = cld(N, _CORR_W) * _CORR_W
            lim = Int16(_MAX_MEAS)                                  # 12-bit signed range ±2048
            meas = zeros(Int16, 2 * Npad)                          # [re₀, im₀, re₁, im₁, …]
            meas[1:2N] .= rand(-lim:(lim - Int16(1)), 2N)
            ext = zeros(Int8, Npad + 2D); csin = zeros(Int8, Npad); ccos = zeros(Int8, Npad)
            # hybrid-blocked: small L1-resident scratch reused across blocks (size independent of N).
            extb = zeros(Int8, _CORR_BLK + 2D); csb = zeros(Int8, _CORR_BLK); ccb = zeros(Int8, _CORR_BLK)

            g = SUITE["correlation"][label]
            g["fused"] = @benchmarkable correlate_epl_fused(
                $meas, $ceng1, $reng, $(Val(_CORR_W)), $(Val(D)), $(Val(_WIPE_TI)),
            ) evals = 1 samples = 1000
            g["unfused"] = @benchmarkable correlate_epl_unfused!(
                $ext, $csin, $ccos, $meas, $code_gen, $tbl, $fs, $freq, $N,
                $(Val(_CORR_W)), $(Val(D)), $(Val(_WIPE_TI)),
            ) evals = 1 samples = 1000
            g["hybrid"] = @benchmarkable correlate_epl_hybrid!(
                $ext, $meas, $code_gen, $reng, $N, $(Val(_CORR_W)), $(Val(D)), $(Val(_WIPE_TI)),
            ) evals = 1 samples = 1000
            g["hybrid-blocked"] = @benchmarkable correlate_epl_hybrid_blocked!(
                $extb, $csb, $ccb, $meas, $code_gen, $reng,
                $(Val(_CORR_W)), $(Val(D)), $(Val(_CORR_BLK)), $(Val(_WIPE_TI)),
            ) evals = 1 samples = 1000
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# PARALLEL multi-channel correlation — FUSED vs UNFUSED vs HYBRID vs HYBRID-BLOCKED across many
# satellites (Polyester).
#
# This is where the fused/unfused trade-off INVERTS. Per channel (single thread) the fused is
# execution-port bound and ~10% slower than the unfused (see the single-channel rows above).
# But the fused materialises NOTHING — it drives the code+carrier engines straight into the
# correlation, so its only memory traffic is reading the shared measurement. The unfused
# instead writes a per-channel extended-code + carrier scratch (≈ 3·N Int8 per channel) and
# reads it back. Run one channel per core and that per-channel scratch — `n_channels · 3 · N`
# bytes in flight at once — eventually overflows the shared last-level cache and saturates
# memory bandwidth: the unfused stops scaling with cores while the fused keeps scaling, because
# it never touches the bus. So the fused WINS once the working set exceeds the LLC.
#
# The HYBRID is the best of those two: it materialises only the CODE (the fast run-fill
# `gen_code!` the fully-fused path can't use) and keeps the carrier in-register, so it touches ONE
# buffer per channel (≈ N Int8) instead of three. It beats fused at moderate channel counts (cheap
# `ext` loads replace the per-chunk value-engine permute + DDA) and, once the working set spills
# the LLC, also beats unfused (a third of the bandwidth). But its per-channel ≈ N code buffer is
# still O(N): push to many channels / long integration (24 ch, 30 ms) and it too spills hard and
# falls BEHIND the zero-buffer fused (measured ~1.25× fused there).
#
# The HYBRID-BLOCKED removes that last buffer: it strip-mines into BLK-sample blocks and
# regenerates BOTH code and carrier into small L1-resident scratch reused across blocks, so its
# working set is O(BLK), not O(N) — fused-class LLC footprint with hybrid-class per-sample cost
# (see the kernel comment above). Measured (40 MHz, 24 threads, Zen5 / AVX-512), it is the fastest
# kernel in every spilling regime: 10 ch 20 ms ~0.85× fused / ~0.69× unfused; 24 ch 20 ms ~0.83×
# fused / ~0.81× hybrid; 24 ch 30 ms ~0.82× fused / ~0.65× hybrid. It also beats fused and hybrid
# when cache-resident (10 ch 5 ms ~0.87× fused / ~0.97× hybrid); the ONLY regime it does not win
# is pure cache-resident vs the unfused (~1.04×), where both do identical work and the unfused's
# two full-length passes are the ceiling. So: hybrid-blocked is the all-round pick, with the
# unfused still marginally ahead only when everything fits in cache.
#
# The crossover is governed by working-set vs LLC size, and the effect only appears with
# `Threads.nthreads() > 1` (run julia with `-t auto` / `JULIA_NUM_THREADS`). Two representative
# points at 40 MHz, 10 channels: 5 ms (≈6 MB, fits L3) and 20 ms (≈24 MB, spills L3). At lower
# sample rates / fewer channels the working set is smaller, so the crossover moves to longer
# integration.
#
# Each channel gets its own PRN, Doppler, code/carrier engines and (unfused / hybrid /
# hybrid-blocked) scratch buffers; the measurement is shared. State is bundled in `_EPLChannels`
# so Polyester sees only a struct and the loop index — it does NOT rewrite the per-channel
# vectors into stride-arrays (which would not match the kernels' `Vector` signatures). One
# channel per `@batch` task.
# ─────────────────────────────────────────────────────────────────────────────
# Defined unconditionally (a struct must be top level); the drivers/registration that use
# Polyester's `@batch` are guarded below.
struct _EPLChannels{W,D,CE,RE,CG,TB}
    meas::Vector{Int16}
    cengs::Vector{CE}
    rengs::Vector{RE}
    exts::Vector{Vector{Int8}}
    extsH::Vector{Vector{Int8}}   # hybrid's own code buffer (so it never aliases the unfused's)
    csins::Vector{Vector{Int8}}
    ccoss::Vector{Vector{Int8}}
    # hybrid-blocked's per-channel L1 scratch (block code + sin + cos), size independent of N
    extbs::Vector{Vector{Int8}}
    csbs::Vector{Vector{Int8}}
    ccbs::Vector{Vector{Int8}}
    code_gens::Vector{CG}
    tbls::Vector{TB}
    fs::Float64
    freqs::Vector{Float64}
    N::Int
    out::Vector{Int}     # per-channel sink (prompt I) — prevents dead-code elimination
end

if _HAS_POLYESTER && isdefined(GNSSSignals, :code_engine)
    # Build `nch` independent channels (distinct PRN + Doppler) sharing one measurement.
    function _build_channels(nch::Int, fs, fc, ms)
        W = _CORR_W
        N = round(Int, fs * 1e-3 * ms)
        D = round(Int, 0.5 * fs / fc)
        Npad = cld(N, W) * W
        lim = Int16(_MAX_MEAS)                                 # 12-bit ADC ±2048
        meas = zeros(Int16, 2 * Npad)
        meas[1:2N] .= rand(-lim:(lim - Int16(1)), 2N)
        plans = [GNSSSignals.CodeReplicaLUT(_GPSL1(), prn) for prn in 1:nch]
        freqs = [1234.0 + 137.0 * (s - 1) for s in 1:nch]      # spread the residual Dopplers
        tbls = [SinCosTable(Int8; steps = 64, amplitude = _CARRIER_AMP) for _ in 1:nch]
        engs = [_epl_engines(plans[s], tbls[s], fs, fc, freqs[s]) for s in 1:nch]
        cengs = [e[1] for e in engs]
        rengs = [e[2] for e in engs]
        code_gens = [GNSSSignals.CodeGeneratorLUT(plans[s], fs, fc) for s in 1:nch]
        exts = [zeros(Int8, Npad + 2D) for _ in 1:nch]
        extsH = [zeros(Int8, Npad + 2D) for _ in 1:nch]
        csins = [zeros(Int8, Npad) for _ in 1:nch]
        ccoss = [zeros(Int8, Npad) for _ in 1:nch]
        extbs = [zeros(Int8, _CORR_BLK + 2D) for _ in 1:nch]
        csbs = [zeros(Int8, _CORR_BLK) for _ in 1:nch]
        ccbs = [zeros(Int8, _CORR_BLK) for _ in 1:nch]
        ch = _EPLChannels{W,D,eltype(cengs),eltype(rengs),eltype(code_gens),eltype(tbls)}(
            meas, cengs, rengs, exts, extsH, csins, ccoss, extbs, csbs, ccbs,
            code_gens, tbls, Float64(fs), freqs, N, zeros(Int, nch))
        ch
    end

    # One channel of work (indexing lives here, not in the `@batch` body, so Polyester leaves
    # the per-channel vectors as plain `Vector`s).
    @inline function _epl_one_fused!(ch::_EPLChannels{W,D}, s) where {W,D}
        @inbounds r = correlate_epl_fused(ch.meas, ch.cengs[s], ch.rengs[s], Val(W), Val(D), Val(_WIPE_TI))
        @inbounds ch.out[s] = r[2][1]
        nothing
    end
    @inline function _epl_one_unfused!(ch::_EPLChannels{W,D}, s) where {W,D}
        @inbounds r = correlate_epl_unfused!(ch.exts[s], ch.csins[s], ch.ccoss[s], ch.meas,
            ch.code_gens[s], ch.tbls[s], ch.fs, ch.freqs[s], ch.N, Val(W), Val(D), Val(_WIPE_TI))
        @inbounds ch.out[s] = r[2][1]
        nothing
    end
    @inline function _epl_one_hybrid!(ch::_EPLChannels{W,D}, s) where {W,D}
        @inbounds r = correlate_epl_hybrid!(ch.extsH[s], ch.meas, ch.code_gens[s], ch.rengs[s],
            ch.N, Val(W), Val(D), Val(_WIPE_TI))
        @inbounds ch.out[s] = r[2][1]
        nothing
    end
    @inline function _epl_one_hybrid_blocked!(ch::_EPLChannels{W,D}, s) where {W,D}
        @inbounds r = correlate_epl_hybrid_blocked!(ch.extbs[s], ch.csbs[s], ch.ccbs[s], ch.meas,
            ch.code_gens[s], ch.rengs[s], Val(W), Val(D), Val(_CORR_BLK), Val(_WIPE_TI))
        @inbounds ch.out[s] = r[2][1]
        nothing
    end
    _epl_fused_par!(ch::_EPLChannels) =
        (Polyester.@batch for s in eachindex(ch.cengs); _epl_one_fused!(ch, s); end; ch.out)
    _epl_unfused_par!(ch::_EPLChannels) =
        (Polyester.@batch for s in eachindex(ch.cengs); _epl_one_unfused!(ch, s); end; ch.out)
    _epl_hybrid_par!(ch::_EPLChannels) =
        (Polyester.@batch for s in eachindex(ch.cengs); _epl_one_hybrid!(ch, s); end; ch.out)
    _epl_hybrid_blocked_par!(ch::_EPLChannels) =
        (Polyester.@batch for s in eachindex(ch.cengs); _epl_one_hybrid_blocked!(ch, s); end; ch.out)

    let nch = 10, fs = 40e6, fc = 1.023e6
        for ms in (5, 20)                              # 5 ms: fits L3; 20 ms: spills L3 (see note)
            ch = _build_channels(nch, fs, fc, ms)
            g = SUITE["correlation"]["GPSL1CA E/P/L 40MHz $(nch)ch $(ms)ms ‖"]
            g["fused"]   = @benchmarkable _epl_fused_par!($ch) evals = 1 samples = 200
            g["unfused"] = @benchmarkable _epl_unfused_par!($ch) evals = 1 samples = 200
            g["hybrid"]  = @benchmarkable _epl_hybrid_par!($ch) evals = 1 samples = 200
            g["hybrid-blocked"] = @benchmarkable _epl_hybrid_blocked_par!($ch) evals = 1 samples = 200
        end
    end
end
