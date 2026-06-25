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
# Correlation signal path — FUSED vs UNFUSED.
#
# Models one tracking correlation with a realistic front-end word layout:
#   measurement — COMPLEX, Int16, re/im interleaved as [reₙ, imₙ, …] (one 12-bit ADC sample
#                 per component, |m| ≲ 2048); a single Vector{Int16} of length 2·N.
#   carrier     — Int8, small configurable amplitude (the complex local carrier cos + j·sin)
#   code        — Int8 ±1 replica (a CBOC/BOC replica can reach ±2)
# The correlator output is the complex sum of measurement × code × conj(cos + j·sin):
#   I (real) = Σ code·(mᵣ·cos + mᵢ·sin),   Q (imag) = Σ code·(mᵢ·cos − mᵣ·sin).
# Per-sample products exceed Int16 (2048·8·2 = 32768), so both arms accumulate into Int32.
# The hot multiply uses `vpmaddwd` (Int16×Int16 → Int32 pairwise multiply-add) where available:
# it forms the product at Int32 width (no overflow, no Int16 staging) and its pairwise reduction
# halves the accumulator lane count, so on AVX-512 the kernel is allocation-free where the plain
# widen-to-Int32 path spills the wide `Vec{64,Int32}` accumulators. The interleaved measurement
# is deinterleaved into mᵣ/mᵢ once per chunk; `mᵣ·code` / `mᵢ·code` (≤ ±4096, fit Int16) feed the
# four pmaddwd terms. Platforms without `vpmaddwd` (NEON/scalar) fall back to a plain Int32
# multiply at the backend's native width; both paths give identical results.
#
#   UNFUSED — materialise the code and carrier into preallocated Int8 buffers (gen_code! +
#             generate_carrier!), then stream the buffers + measurement through the
#             correlation loop. Extra memory passes; buffers preallocated once.
#   FUSED   — never materialise code/carrier: drive the value-based code/carrier engines
#             (`code_engine` + `carrier_engine`) straight into the correlation loop. Both are
#             loop-invariant and built once; the hot loop holds 4 isbits states W samples apart
#             (built with `Val(4)` for code, `carrier_state(eng, kW)` for carrier) so four
#             DDA/NCO carry chains overlap (full ILP, fastest). Each `code_lookup`/`carrier_lookup`
#             returns an Int8 SIMD chunk in registers; the states are renewed by value with
#             `code_advance`/`carrier_advance` — no heap, nothing materialised. The < 4·W tail is
#             finished by single-wide engines (code `Val(1)` positioned at the tail's start
#             sample; the carrier engine is interleave-agnostic and reused with `carrier_advance(…, 1)`),
#             and the final < W samples are absorbed by zero-padding `meas` to a W multiple
#             (padding · anything = 0) — so no slow scalar remainder.
#
# GPS L1 C/A, 5 MHz, 1 ms ⇒ N = 5000 samples. Carrier amplitude and ADC bit depth are
# configurable (the `carrier_amplitude` / `adc_bits` locals below). The plan (CodeReplicaLUT),
# SinCosTable, the value-based engines, the warm unfused generator and the buffers are all built
# once in setup (a receiver reuses them); the timed region only creates the cheap isbits states
# and runs the correlation loop.
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

# `pmaddwd` exists on AVX2 (256-bit) and AVX-512BW (512-bit) — i.e. whenever the Int8
# code/carrier backend is AVX2 (W=32) or AVX-512 (W=64). NEON has no direct equivalent
# (its `smull`/`smlal` widening multiply is what the portable Int32 path lowers to), and
# the scalar fallback (W=1) has none either, so both use `_accum_portable!`.
const _HAS_VPMADDWD = _CORR_W == 32 || _CORR_W == 64

@inline _wide16(v::Vec{W,Int8}) where {W} = convert(Vec{W,Int16}, v)
@inline _wide32(v::Vec{W,Int8}) where {W} = convert(Vec{W,Int32}, v)
@inline _wide32(v::Vec{W,Int16}) where {W} = convert(Vec{W,Int32}, v)

# Lower / upper half of a vector (free: selects one register half of the wider vector).
@inline _lo(v::Vec{W,T}) where {W,T} = shufflevector(v, Val(ntuple(i -> i - 1, Val(W ÷ 2))))
@inline _hi(v::Vec{W,T}) where {W,T} = shufflevector(v, Val(ntuple(i -> i - 1 + W ÷ 2, Val(W ÷ 2))))

# Deinterleave [re, im, …] → real (even lanes) / imag (odd lanes), Vec{2W,Int16} → Vec{W,Int16}.
@inline _re(v::Vec{M,Int16}) where {M} = shufflevector(v, Val(ntuple(i -> 2 * (i - 1), Val(M ÷ 2))))
@inline _im(v::Vec{M,Int16}) where {M} = shufflevector(v, Val(ntuple(i -> 2 * (i - 1) + 1, Val(M ÷ 2))))

# pmaddwd: a,b Int16 → Int32, dst[j] = a[2j]·b[2j] + a[2j+1]·b[2j+1] (each product formed at
# Int32 width, so no Int16 overflow). Reached via llvmcall as SIMD.jl does not expose it; one
# variant per register width. Only *called* on the matching backend (so the avx512 form is
# never lowered on an AVX2 host and vice versa), but defining both unconditionally is harmless.
@inline function _vpmaddwd(a::Vec{32,Int16}, b::Vec{32,Int16})   # AVX-512BW (512-bit)
    Vec(Base.llvmcall(
        (raw"""
         declare <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16>, <32 x i16>)
         define <16 x i32> @entry(<32 x i16> %a, <32 x i16> %b) #0 {
         top:
           %r = call <16 x i32> @llvm.x86.avx512.pmaddw.d.512(<32 x i16> %a, <32 x i16> %b)
           ret <16 x i32> %r
         }
         attributes #0 = { "target-features"="+avx512f,+avx512bw,+avx512vl" }
         """, "entry"),
        NTuple{16,Base.VecElement{Int32}},
        Tuple{NTuple{32,Base.VecElement{Int16}},NTuple{32,Base.VecElement{Int16}}},
        a.data, b.data))
end
@inline function _vpmaddwd(a::Vec{16,Int16}, b::Vec{16,Int16})   # AVX2 (256-bit)
    Vec(Base.llvmcall(
        (raw"""
         declare <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16>, <16 x i16>)
         define <8 x i32> @entry(<16 x i16> %a, <16 x i16> %b) #0 {
         top:
           %r = call <8 x i32> @llvm.x86.avx2.pmadd.wd(<16 x i16> %a, <16 x i16> %b)
           ret <8 x i32> %r
         }
         attributes #0 = { "target-features"="+avx2" }
         """, "entry"),
        NTuple{8,Base.VecElement{Int32}},
        Tuple{NTuple{16,Base.VecElement{Int16}},NTuple{16,Base.VecElement{Int16}}},
        a.data, b.data))
end

# W-wide Int16 multiply-add → W/2 Int32 lanes (two pmaddwd over the register halves, concatenated).
@inline _madd(a::Vec{W,Int16}, b::Vec{W,Int16}) where {W} = shufflevector(
    _vpmaddwd(_lo(a), _lo(b)), _vpmaddwd(_hi(a), _hi(b)), Val(ntuple(i -> i - 1, Val(W ÷ 2))))

# Complex pmaddwd accumulation. `v` loads 2W interleaved Int16 (W complex samples); deinterleave
# to mᵣ/mᵢ (once), form mᵣ·code / mᵢ·code in Int16 (≤ ±4096), then four pmaddwd against the
# widened carrier accumulate the complex output in Int32:
#   I += code·(mᵣ·cos + mᵢ·sin),   Q += code·(mᵢ·cos − mᵣ·sin).
macro _accum_madd!(AI, AQ, v, sinv, cosv, code)
    quote
        vv = $(esc(v)); cd = _wide16($(esc(code)))
        mcr = _re(vv) * cd; mci = _im(vv) * cd
        cs = _wide16($(esc(cosv))); sn = _wide16($(esc(sinv)))
        $(esc(AI)) += _madd(mcr, cs) + _madd(mci, sn)
        $(esc(AQ)) += _madd(mci, cs) - _madd(mcr, sn)
    end
end

# Portable complex accumulation (any width): widen everything to Int32 and multiply there —
# no overflow, no intrinsic. Same I/Q expressions as the pmaddwd path.
macro _accum_portable!(AI, AQ, v, sinv, cosv, code)
    quote
        vv = $(esc(v)); cd = _wide32($(esc(code)))
        mcr = _wide32(_re(vv)) * cd; mci = _wide32(_im(vv)) * cd
        cs = _wide32($(esc(cosv))); sn = _wide32($(esc(sinv)))
        $(esc(AI)) += mcr * cs + mci * sn
        $(esc(AQ)) += mci * cs - mcr * sn
    end
end

# ── FUSED kernels: code × carrier in registers, consumed in place (no array materialised) ──
# The value-based engines (`code_engine`/`carrier_engine`) are loop-invariant and built once in
# setup, then passed in: dispatch specialises each kernel on the concrete engine types (a function
# barrier), so the renew-by-value state stepping stays allocation-free even where the outer
# inference is abstract. `ceng4`/`ceng1` are the K=4 (4-way interleaved) and K=1 (tail) code
# engines; `reng` is the single, interleave-agnostic carrier engine (its `carrier_advance` takes
# an explicit chunk count, so the same engine drives both the 4-way main loop and the 1-wide tail).
#
# pmaddwd path (AVX2 W=32 / AVX-512 W=64): Int32 accumulators with halved lane count (W÷2).
function correlate_fused_madd(meas::Vector{Int16}, ceng4, ceng1, reng, N::Int, ::Val{W}) where {W}
    AI = zero(Vec{W ÷ 2,Int32}); AQ = zero(Vec{W ÷ 2,Int32})
    n4 = (N ÷ (4W)) * (4W)
    # 4-way interleaved main loop: four code/carrier streams W samples apart → carry chains overlap.
    cs1 = code_state(ceng4, 0); cs2 = code_state(ceng4, 1); cs3 = code_state(ceng4, 2); cs4 = code_state(ceng4, 3)
    rs1 = carrier_state(reng, 0); rs2 = carrier_state(reng, W); rs3 = carrier_state(reng, 2W); rs4 = carrier_state(reng, 3W)
    n = 0
    @inbounds for _ in 1:(N ÷ (4W))
        s1, co1 = carrier_lookup(reng, rs1)
        @_accum_madd! AI AQ vload(Vec{2W,Int16}, meas, 2n + 1)        s1 co1 code_lookup(ceng4, cs1)
        s2, co2 = carrier_lookup(reng, rs2)
        @_accum_madd! AI AQ vload(Vec{2W,Int16}, meas, 2(n + W) + 1)  s2 co2 code_lookup(ceng4, cs2)
        s3, co3 = carrier_lookup(reng, rs3)
        @_accum_madd! AI AQ vload(Vec{2W,Int16}, meas, 2(n + 2W) + 1) s3 co3 code_lookup(ceng4, cs3)
        s4, co4 = carrier_lookup(reng, rs4)
        @_accum_madd! AI AQ vload(Vec{2W,Int16}, meas, 2(n + 3W) + 1) s4 co4 code_lookup(ceng4, cs4)
        n += 4W
        cs1 = code_advance(ceng4, cs1); cs2 = code_advance(ceng4, cs2)
        cs3 = code_advance(ceng4, cs3); cs4 = code_advance(ceng4, cs4)
        rs1 = carrier_advance(reng, rs1, 4); rs2 = carrier_advance(reng, rs2, 4)
        rs3 = carrier_advance(reng, rs3, 4); rs4 = carrier_advance(reng, rs4, 4)
    end
    # < 4·W tail, single-wide engines positioned at sample n4; sub-W remainder rides on meas's pad.
    cst = code_state(ceng1, n4 ÷ W); rst = carrier_state(reng, n4)
    @inbounds while n < length(meas) ÷ 2
        sv, cosv = carrier_lookup(reng, rst)
        @_accum_madd! AI AQ vload(Vec{2W,Int16}, meas, 2n + 1) sv cosv code_lookup(ceng1, cst)
        n += W
        cst = code_advance(ceng1, cst); rst = carrier_advance(reng, rst, 1)
    end
    (sum(AI), sum(AQ))
end

# Portable fallback (NEON / scalar): plain Int32 multiply at the backend's width.
function correlate_fused_portable(meas::Vector{Int16}, ceng4, ceng1, reng, N::Int, ::Val{W}) where {W}
    AI = zero(Vec{W,Int32}); AQ = zero(Vec{W,Int32})
    n4 = (N ÷ (4W)) * (4W)
    cs1 = code_state(ceng4, 0); cs2 = code_state(ceng4, 1); cs3 = code_state(ceng4, 2); cs4 = code_state(ceng4, 3)
    rs1 = carrier_state(reng, 0); rs2 = carrier_state(reng, W); rs3 = carrier_state(reng, 2W); rs4 = carrier_state(reng, 3W)
    n = 0
    @inbounds for _ in 1:(N ÷ (4W))
        s1, co1 = carrier_lookup(reng, rs1)
        @_accum_portable! AI AQ vload(Vec{2W,Int16}, meas, 2n + 1)        s1 co1 code_lookup(ceng4, cs1)
        s2, co2 = carrier_lookup(reng, rs2)
        @_accum_portable! AI AQ vload(Vec{2W,Int16}, meas, 2(n + W) + 1)  s2 co2 code_lookup(ceng4, cs2)
        s3, co3 = carrier_lookup(reng, rs3)
        @_accum_portable! AI AQ vload(Vec{2W,Int16}, meas, 2(n + 2W) + 1) s3 co3 code_lookup(ceng4, cs3)
        s4, co4 = carrier_lookup(reng, rs4)
        @_accum_portable! AI AQ vload(Vec{2W,Int16}, meas, 2(n + 3W) + 1) s4 co4 code_lookup(ceng4, cs4)
        n += 4W
        cs1 = code_advance(ceng4, cs1); cs2 = code_advance(ceng4, cs2)
        cs3 = code_advance(ceng4, cs3); cs4 = code_advance(ceng4, cs4)
        rs1 = carrier_advance(reng, rs1, 4); rs2 = carrier_advance(reng, rs2, 4)
        rs3 = carrier_advance(reng, rs3, 4); rs4 = carrier_advance(reng, rs4, 4)
    end
    cst = code_state(ceng1, n4 ÷ W); rst = carrier_state(reng, n4)
    @inbounds while n < length(meas) ÷ 2
        sv, cosv = carrier_lookup(reng, rst)
        @_accum_portable! AI AQ vload(Vec{2W,Int16}, meas, 2n + 1) sv cosv code_lookup(ceng1, cst)
        n += W
        cst = code_advance(ceng1, cst); rst = carrier_advance(reng, rst, 1)
    end
    (sum(AI), sum(AQ))
end

# Build the value-based engines once (a receiver reuses them across integrations): K=4 / K=1 code
# engines over the plan + the carrier NCO engine. `code_engine`'s parameterless `default_backend()`
# const-folds to the host backend, so the engine types are concrete and the kernels stay 0-alloc.
@inline function _fused_engines(plan, tbl, fs, fc, freq)
    (GNSSSignals.code_engine(plan, fs, fc, Val(4)),
     GNSSSignals.code_engine(plan, fs, fc, Val(1)),
     carrier_engine(tbl; frequency = freq, sampling_frequency = fs))
end

# ── UNFUSED kernels: materialise code + carrier into preallocated buffers, then correlate.
# The code generator is prebuilt (warm) and reused so the timed region allocates nothing —
# the comparison is purely register-fusion vs the extra memory passes of materialise-then-read.
function correlate_unfused_madd!(code::Vector{Int8}, csin::Vector{Int8}, ccos::Vector{Int8},
                                 meas::Vector{Int16}, code_gen, tbl, fs, freq, N::Int, ::Val{W}) where {W}
    gen_code!(view(code, 1:N), code_gen)
    generate_carrier!(view(csin, 1:N), view(ccos, 1:N), tbl; frequency = freq, sampling_frequency = fs)
    AI = zero(Vec{W ÷ 2,Int32}); AQ = zero(Vec{W ÷ 2,Int32}); n = 0
    @inbounds while n < length(code)
        @_accum_madd! AI AQ vload(Vec{2W,Int16}, meas, 2n + 1) vload(Vec{W,Int8}, csin, n + 1) vload(Vec{W,Int8}, ccos, n + 1) vload(Vec{W,Int8}, code, n + 1)
        n += W
    end
    (sum(AI), sum(AQ))
end

function correlate_unfused_portable!(code::Vector{Int8}, csin::Vector{Int8}, ccos::Vector{Int8},
                                     meas::Vector{Int16}, code_gen, tbl, fs, freq, N::Int, ::Val{W}) where {W}
    gen_code!(view(code, 1:N), code_gen)
    generate_carrier!(view(csin, 1:N), view(ccos, 1:N), tbl; frequency = freq, sampling_frequency = fs)
    AI = zero(Vec{W,Int32}); AQ = zero(Vec{W,Int32}); n = 0
    @inbounds while n < length(code)
        @_accum_portable! AI AQ vload(Vec{2W,Int16}, meas, 2n + 1) vload(Vec{W,Int8}, csin, n + 1) vload(Vec{W,Int8}, ccos, n + 1) vload(Vec{W,Int8}, code, n + 1)
        n += W
    end
    (sum(AI), sum(AQ))
end

# Register the benchmarks only where the value-based code engine exists (skipped on a
# baseline without the code-LUT API, so AirspeedVelocity can still diff against it).
if isdefined(GNSSSignals, :code_engine)
    let
        fs = 5e6           # 5 MHz sampling
        fc = 1.023e6       # GPS L1 C/A chip rate
        freq = 1234.0      # residual carrier / Doppler to wipe off (Hz)
        N = 5000           # 1 ms at 5 MHz
        carrier_amplitude = 8   # carrier ∈ ±amplitude (configurable, stored Int8)
        adc_bits = 12           # measurement is a 12-bit ADC sample, stored Int16
        plan = GNSSSignals.CodeReplicaLUT(_GPSL1(), 1)
        code_gen = GNSSSignals.CodeGeneratorLUT(plan, fs, fc)   # warm generator for the unfused fill
        tbl = SinCosTable(Int8; steps = 64, amplitude = carrier_amplitude)

        # Value-based fused engines, built once (loop-invariant): K=4 + K=1 code engines + carrier.
        ceng4, ceng1, reng = _fused_engines(plan, tbl, fs, fc, freq)

        # Zero-pad to a whole number of W-wide chunks (Npad complex samples) so the fused tail
        # and the unfused loop never read past the end and the sub-W remainder contributes 0.
        # The measurement is complex with re/im interleaved, so it holds 2·Npad Int16 words.
        Npad = cld(N, _CORR_W) * _CORR_W
        lim = Int16(1) << (adc_bits - 1)                        # 12-bit signed range ±2048
        meas = zeros(Int16, 2 * Npad)                           # [re₀, im₀, re₁, im₁, …]
        meas[1:2N] .= rand(-lim:(lim - Int16(1)), 2N)
        code = zeros(Int8, Npad); csin = zeros(Int8, Npad); ccos = zeros(Int8, Npad)

        g = SUITE["correlation"]["GPSL1CA 5MHz 1ms"]
        if _HAS_VPMADDWD
            g["fused"] = @benchmarkable correlate_fused_madd(
                $meas, $ceng4, $ceng1, $reng, $N, $(Val(_CORR_W)),
            ) evals = 1 samples = 1000
            g["unfused"] = @benchmarkable correlate_unfused_madd!(
                $code, $csin, $ccos, $meas, $code_gen, $tbl, $fs, $freq, $N, $(Val(_CORR_W)),
            ) evals = 1 samples = 1000
        else
            g["fused"] = @benchmarkable correlate_fused_portable(
                $meas, $ceng4, $ceng1, $reng, $N, $(Val(_CORR_W)),
            ) evals = 1 samples = 1000
            g["unfused"] = @benchmarkable correlate_unfused_portable!(
                $code, $csin, $ccos, $meas, $code_gen, $tbl, $fs, $freq, $N, $(Val(_CORR_W)),
            ) evals = 1 samples = 1000
        end
    end
end
