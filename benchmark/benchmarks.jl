using BenchmarkTools
using GNSSSignals
using SIMD: Vec, vload, shufflevector
import SinCosLUT
using SinCosLUT: SinCosTable, generate_carrier!, cycles_per_sample,
                 carrier_engine, carrier_state, carrier_lookup, carrier_advance, carrier_width
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
# Oversampling sweep — old gen_code! vs new LUT across oversampling ratios, at a small
# and a steady-state buffer. "Oversampling ratio" = samples per code chip = fs / fc;
# the level is labelled as a multiplier (2x = sample twice per chip). Grouped by signal /
# oversampling / size under `code/oversampling sweep/…` so "original" and "LUT" sit
# adjacent. Same N + fs/fc for both, and the LUT side uses the warm generator (0-alloc),
# matching the original's 0-alloc fill. The LUT is ~flat in the oversampling ratio (one
# permute/sample); the original's run-fill speeds up as it grows — so the LUT wins most at
# low oversampling and the original catches up high (crossover ~8-16x BPSK, later BOC).
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
# independent correlations. The wipe is a per-sample complex multiply: deinterleave mᵣ/mᵢ,
# widen everything to Int32 and multiply (per-sample products reach 2048·8·2 = 32768 > Int16,
# so DI/DQ are Int32). The code step is then an Int32 multiply (code ±1, CBOC/BOC ±2) into the
# accumulators. NOTE: `vpmaddwd` is deliberately NOT used here. It packs an Int16×Int16→Int32
# *pairwise sum*, which only pays off when the output is a reduced sum; the wipe output is
# per-sample, so the interleave + shuffles `vpmaddwd` needs are pure overhead — the straight
# deinterleave-and-multiply is ~17–25 % faster (measured on AVX-512) and identical on NEON,
# which has no `vpmaddwd` anyway. (Keeping DI/DQ in Int16 to re-enable `vpmaddwd` for the code
# step was tried and is a wash: the kernel is bound by the lookup/shuffle work, not the
# multiply width.)
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

@inline _wide32(v::Vec{W,Int8}) where {W} = convert(Vec{W,Int32}, v)
@inline _wide32(v::Vec{W,Int16}) where {W} = convert(Vec{W,Int32}, v)

# Deinterleave [re, im, …] → real (even lanes) / imag (odd lanes), Vec{2W,Int16} → Vec{W,Int16}.
@inline _re(v::Vec{M,Int16}) where {M} = shufflevector(v, Val(ntuple(i -> 2 * (i - 1), Val(M ÷ 2))))
@inline _im(v::Vec{M,Int16}) where {M} = shufflevector(v, Val(ntuple(i -> 2 * (i - 1) + 1, Val(M ÷ 2))))

# Per-sample shared carrier wipe (interleaved meas + sin/cos Int8) → (DI, DQ)::Vec{W,Int32}:
# deinterleave mᵣ/mᵢ, widen to Int32, multiply. Same expressions on every backend (no intrinsic).
@inline function _wipe(v::Vec{M,Int16}, sinv, cosv) where {M}
    mr = _wide32(_re(v)); mi = _wide32(_im(v)); cs = _wide32(cosv); sn = _wide32(sinv)
    (mr * cs + mi * sn, mi * cs - mr * sn)
end

# Accumulate one chunk into the six Int32 E/P/L accumulators (code ±1/±2 × Int32 DI/DQ).
@inline _accum_epl(IE, QE, IP, QP, IL, QL, cE, cP, cL, DI, DQ) =
    (IE + cE * DI, QE + cE * DQ, IP + cP * DI, QP + cP * DQ, IL + cL * DI, QL + cL * DQ)

# ── FUSED E/P/L: one Prompt code lookup per chunk; Early/Late are lane-shifts across a sliding
# pipeline of P = ⌈D/W⌉+⌊D/W⌋+2 carried code chunks. `@generated` so the pipeline length, the
# per-correlator chunk indices and the shuffle masks fold to literals (fully unrolled, 0-alloc)
# for any (W, D). The Prompt code is a single value-based stream; nothing is materialised.
@generated function correlate_epl_fused(meas::Vector{Int16}, ceng1, reng, ::Val{W}, ::Val{D}) where {W,D}
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
        IE = zero(Vec{$W,Int32}); QE = zero(Vec{$W,Int32}); IP = zero(Vec{$W,Int32})
        QP = zero(Vec{$W,Int32}); IL = zero(Vec{$W,Int32}); QL = zero(Vec{$W,Int32})
        cs = code_state(ceng1, 0); rs = carrier_state(reng, 0); z = zero(Vec{$W,Int8})
        $(prime...)
        n = 0; lim = length(meas) ÷ 2
        @inbounds while n < lim
            sinv, cosv = carrier_lookup(reng, rs); rs = carrier_advance(reng, rs, 1)
            DI, DQ = _wipe(vload(Vec{$(2W),Int16}, meas, 2n + 1), sinv, cosv)
            early = shufflevector($(csym[1]), $(csym[2]), Val($emask))
            late = shufflevector($(csym[g+fl+1]), $(csym[g+fl+2]), Val($lmask))
            cE = _wide32(early); cP = _wide32($(csym[g+1])); cL = _wide32(late)
            IE, QE, IP, QP, IL, QL = _accum_epl(IE, QE, IP, QP, IL, QL, cE, cP, cL, DI, DQ)
            $(shift...)
            n += $W
        end
        ((sum(IE), sum(QE)), (sum(IP), sum(QP)), (sum(IL), sum(QL)))
    end
end

# ── UNFUSED E/P/L: materialise the extended Prompt code (samples −D … N−1+D, leading D zero)
# + the shared carrier, then run the shared-wipe loop with Early/Prompt/Late as offset reads.
# `ext[k]` holds the code at sample k−1−D, so at output sample n Early/Prompt/Late read
# ext[n+1] / ext[n+D+1] / ext[n+2D+1]. Generator is warm (continuing); timed region is 0-alloc.
function correlate_epl_unfused!(ext::Vector{Int8}, csin::Vector{Int8}, ccos::Vector{Int8},
                                meas::Vector{Int16}, code_gen, tbl, fs, freq, N::Int,
                                ::Val{W}, ::Val{D}) where {W,D}
    gen_code!(view(ext, (D + 1):(D + N + D)), code_gen)     # samples 0 … N+D−1 (ext[1:D] stay 0)
    generate_carrier!(view(csin, 1:N), view(ccos, 1:N), tbl; frequency = freq, sampling_frequency = fs)
    IE = zero(Vec{W,Int32}); QE = zero(Vec{W,Int32}); IP = zero(Vec{W,Int32})
    QP = zero(Vec{W,Int32}); IL = zero(Vec{W,Int32}); QL = zero(Vec{W,Int32})
    n = 0
    @inbounds while n < length(csin)
        DI, DQ = _wipe(vload(Vec{2W,Int16}, meas, 2n + 1),
                       vload(Vec{W,Int8}, csin, n + 1), vload(Vec{W,Int8}, ccos, n + 1))
        cE = _wide32(vload(Vec{W,Int8}, ext, n + 1))
        cP = _wide32(vload(Vec{W,Int8}, ext, n + D + 1))
        cL = _wide32(vload(Vec{W,Int8}, ext, n + 2D + 1))
        IE, QE, IP, QP, IL, QL = _accum_epl(IE, QE, IP, QP, IL, QL, cE, cP, cL, DI, DQ)
        n += W
    end
    ((sum(IE), sum(QE)), (sum(IP), sum(QP)), (sum(IL), sum(QL)))
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
            carrier_amplitude = 8                      # carrier ∈ ±amplitude (stored Int8)
            adc_bits = 12                              # measurement is a 12-bit ADC sample (Int16)
            plan = GNSSSignals.CodeReplicaLUT(_GPSL1(), 1)
            code_gen = GNSSSignals.CodeGeneratorLUT(plan, fs, fc)   # warm generator for the unfused fill
            tbl = SinCosTable(Int8; steps = 64, amplitude = carrier_amplitude)
            ceng1, reng = _epl_engines(plan, tbl, fs, fc, freq)

            # Zero-pad to a whole number of W-wide chunks (Npad samples). meas is complex with
            # re/im interleaved → 2·Npad Int16 words. ext spans samples −D … (Npad−1)+D.
            Npad = cld(N, _CORR_W) * _CORR_W
            lim = Int16(1) << (adc_bits - 1)                       # 12-bit signed range ±2048
            meas = zeros(Int16, 2 * Npad)                          # [re₀, im₀, re₁, im₁, …]
            meas[1:2N] .= rand(-lim:(lim - Int16(1)), 2N)
            ext = zeros(Int8, Npad + 2D); csin = zeros(Int8, Npad); ccos = zeros(Int8, Npad)

            g = SUITE["correlation"][label]
            g["fused"] = @benchmarkable correlate_epl_fused(
                $meas, $ceng1, $reng, $(Val(_CORR_W)), $(Val(D)),
            ) evals = 1 samples = 1000
            g["unfused"] = @benchmarkable correlate_epl_unfused!(
                $ext, $csin, $ccos, $meas, $code_gen, $tbl, $fs, $freq, $N,
                $(Val(_CORR_W)), $(Val(D)),
            ) evals = 1 samples = 1000
        end
    end
end
