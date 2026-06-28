# ─────────────────────────────────────────────────────────────────────────────
# Fast LUT-based code generation (`CodeReplicaLUT` plan + `gen_code!` method).
#
# This file vendors the GNSSSignalsLUT.jl machinery as an internal submodule
# `CodeLUT` (a verbatim copy of its permute / kernel / iterate / modulation
# sources plus its backend-selection + `CodeTable` definition) and adds a thin
# adapter at GNSSSignals top level: a `CodeReplicaLUT` plan (built once per
# signal/PRN) plus `gen_code!` / `gen_code` methods dispatching on it, matching
# the `gen_code!` signature.
#
# The LUT resampler bakes the BOC/TMBOC subcarrier (and short secondary codes)
# into an expanded ±1 Int8 table and resamples it with a single AVX-512 `vpermb`
# / AVX2 `vpshufb` sliding-window permute over a drift-free integer DDA — or, once the
# baked table is heavily oversampled (so consecutive samples repeat a chip), a broadcast
# run-fill that matches the original `gen_code!`'s store-bound speed instead of paying a
# permute per window. The baked table is Int8: ±1 for BPSK/BOC/TMBOC, or a multi-level
# integer approximation of the sqrt-power amplitudes for CBOC (Galileo E1B); cosine-BOC is
# unsupported. Requires sub-chip oversampling (`sampling_frequency ≥ code_frequency · subchip_factor`).
# ─────────────────────────────────────────────────────────────────────────────

# `GNSSSignals.CodeLUT` — internal submodule vendoring the GNSSSignalsLUT.jl SIMD code
# resampler. Not exported and not part of the public GNSSSignals API; use `CodeReplicaLUT`
# with `gen_code!` / `gen_code` instead. All names (`CodeTable`, `LOC`, `BOC`, `TMBOC`,
# `CBOC`, `ModulatedCode`, `code_replica`, `generate_code!`, `generate_code`, `AVX512`, `AVX2`,
# `Portable`, `default_backend`, …) live inside this submodule so they do not clash with
# GNSSSignals' own `LOC` / `BOC` / `TMBOC` / `CBOC` / `Modulation`. (Plain comment, not a docstring,
# so Documenter's checkdocs doesn't require this internal module in the manual.)
module CodeLUT

using SIMD

# ---- backends ----
abstract type Backend end
struct AVX512   <: Backend end   # vpermb over a 64-chip sliding window (W = 64)
struct AVX2     <: Backend end   # vpshufb over two independent 16-chip windows (W = 32)
struct Neon     <: Backend end   # tbl1 over a single 16-chip window (W = 16, AArch64)
struct Portable <: Backend end   # scalar fallback (any CPU)

backend_name(::AVX512)   = "AVX-512"
backend_name(::AVX2)     = "AVX2"
backend_name(::Neon)     = "NEON"
backend_name(::Portable) = "portable"

# Window width / SIMD lane count per backend.
_vwidth(::AVX512)   = Val(64)
_vwidth(::AVX2)     = Val(32)
_vwidth(::Neon)     = Val(16)
_vwidth(::Portable) = Val(1)

"""
    CodeTable(chips::AbstractVector{<:Integer})

Holds a GNSS spreading code of length `L` as `Int8` chips (values are stored verbatim;
pass ±1 for a standard correlation replica). Internally keeps a copy padded with its own
first 63 chips so any 64-chip window `chips[base : base+63]` is a single contiguous load.
"""
struct CodeTable
    chips::Vector{Int8}    # length L
    padded::Vector{Int8}   # length L + WINDOW_PAD
    length::Int
end

# vpermb reads a 64-chip window; the last valid base is L-1, reading up to index L+62
# (0-based). Pad by 63 so that load is always in-bounds.
const WINDOW_PAD = 63

function CodeTable(chips::AbstractVector{<:Integer})
    L = length(chips)
    L > 0 || throw(ArgumentError("code length must be positive"))
    c = Int8.(chips)
    padded = vcat(c, c[1:min(WINDOW_PAD, L)])
    # if L < WINDOW_PAD, repeat until we have L + WINDOW_PAD entries
    while length(padded) < L + WINDOW_PAD
        padded = vcat(padded, c[1:min(L + WINDOW_PAD - length(padded), L)])
    end
    CodeTable(c, padded, L)
end

Base.length(t::CodeTable) = t.length

include("code_lut/permute.jl")
include("code_lut/kernel.jl")
include("code_lut/iterate.jl")
include("code_lut/modulation.jl")
include("code_lut/generator.jl")

# ---- backend selection ----
@static if Sys.ARCH in (:x86_64, :i686)
    function default_backend()
        HOST_FEATURES.avx512vbmi ? AVX512() :
        HOST_FEATURES.avx2       ? AVX2()   : Portable()
    end
elseif Sys.ARCH === :aarch64
    default_backend() = Neon()
else
    default_backend() = Portable()
end
# Length-aware: AVX2/NEON widen the phase vector to Int32 for tables > typemax(Int16), so
# they address anything up to typemax(Int32) (slower than the Int16 path, but ~20× over
# scalar). Only fall back to Portable for the (unreachable for GNSS) even-longer tables.
function default_backend(table::CodeTable)
    be = default_backend()
    (be isa Union{AVX2,Neon} && table.length > typemax(Int32)) ? Portable() : be
end

end # module CodeLUT

# ─────────────────────────────────────────────────────────────────────────────
# Adapter: `CodeReplicaLUT` plan + `gen_code!` / `gen_code` methods at top level.
# ─────────────────────────────────────────────────────────────────────────────

# Per-length Int8 scratch buffers for non-Int8 output (e.g. Int16), keyed by
# length, to avoid a per-call allocation.
const _GEN_CODE_LUT_SCRATCH = Dict{Int,Vector{Int8}}()

# Strip Unitful to a plain Float64 Hz value; pass plain numbers through.
@inline _to_hz(x::Frequency) = Float64(ustrip(u"Hz", x))
@inline _to_hz(x) = Float64(x)

# Map a GNSSSignals modulation to a CodeLUT modulation. ±1 for LOC/BOC/TMBOC; CBOC bakes a
# multi-level Int8 integer approximation of its sqrt-power amplitudes (`cboc_amplitudes`).
# Errors on cosine BOC and code-factor n ≠ 1. The `cboc_amplitudes` argument is ignored by
# every modulation except CBOC (the generic two-arg method below forwards to the one-arg form).
_codelut_modulation(m, ::Tuple{Integer,Integer}) = _codelut_modulation(m)
function _codelut_modulation(m::LOC)
    CodeLUT.LOC()
end
function _codelut_modulation(m::BOCsin)
    m.n == 1 || error("CodeReplicaLUT supports only code-factor n==1 BOC; use gen_code!")
    CodeLUT.BOC(Int(m.m))
end
function _codelut_modulation(m::TMBOC)
    (m.boc1.n == 1 && m.boc2.n == 1 && m.boc1.m == 1) ||
        error("CodeReplicaLUT supports only TMBOC with BOC(1,1) base and code-factor n==1; use gen_code!")
    CodeLUT.TMBOC(Int(m.boc2.m), collect(Bool, m.pattern))
end
function _codelut_modulation(m::CBOC, cboc_amplitudes::Tuple{Integer,Integer})
    (m.boc1.n == 1 && m.boc2.n == 1) ||
        error("CodeReplicaLUT supports only code-factor n==1 CBOC; use gen_code!")
    a1, a2 = cboc_amplitudes
    (a1 > 0 && a2 > 0) ||
        error("CodeReplicaLUT cboc_amplitudes must be positive integers; got $cboc_amplitudes")
    Int(a1) + Int(a2) <= typemax(Int8) ||
        error("CodeReplicaLUT cboc_amplitudes must satisfy a1 + a2 ≤ $(typemax(Int8)) (Int8 table); got $cboc_amplitudes")
    CodeLUT.CBOC(Int(m.boc1.m), Int(m.boc2.m), Int8(a1), Int8(a2))
end
function _codelut_modulation(m::BOCcos)
    error("CodeReplicaLUT does not support cosine-phased BOC (Int8/±1 only); use gen_code!")
end

"""
    CodeReplicaLUT(signal::AbstractGNSSSignal, prn::Integer; cboc_amplitudes = (19, 6))

A reusable, rate-independent plan for fast LUT-based code generation of
`signal`/`prn`. Building it bakes the fully-modulated replica (primary ×
subcarrier × short secondary) into an expanded `Int8` `CodeTable` once — this is
the expensive step. The resulting plan is **immutable and read-only**, so it can
be shared across threads (build once per `(signal, prn)`, reuse it for every
integration; a receiver holds one plan per tracked channel).

Pass the plan to [`gen_code!`](@ref) or [`CodeGeneratorLUT`](@ref) with a sampling
frequency to resample it. The baked table is independent of the sampling and code
frequency, so a single plan serves any rate.

# CBOC (Galileo E1B)
CBOC signals are supported via an **Int8 integer approximation** of the two
sqrt-power subcarrier amplitudes. `cboc_amplitudes = (a1, a2)` sets the integer
amplitudes of the `BOC(1,1)` and `BOC(6,1)` components respectively; the baked
sub-chip values are the composite `±(a1 ± a2)`. For a correlation replica only the
*ratio* `a1/a2` matters (overall scale is irrelevant): the subcarrier correlation
is exactly the dot product of the amplitude vectors, so the match to the float
spec is the cosine between `(a1, a2)` and the true `(sqrt(10/11), sqrt(1/11)) ∝
(sqrt(10), 1) = (3.162, 1)`. The default `(19, 6)` (ratio `3.167`, baked values
`±25, ±13`) reproduces the float subcarrier to ~0 dB correlation loss; coarser
choices trade accuracy for smaller magnitudes (`(3, 1)` → −0.001 dB; `(2, 1)` →
−0.108 dB). Any `a1, a2 > 0` with `a1 + a2 ≤ 127` is accepted. Because the
amplitudes only scale the replica, the **sign pattern is identical to the float
`gen_code!`** for any `a1 > a2 > 0`. The argument is ignored for non-CBOC signals.

# Limitations (vs the plain-`signal` `gen_code!`)
- **Int8 output.** Non-CBOC signals are ±1; CBOC is the multi-level integer
  approximation above. Cosine-phased BOC raises an error at construction (use
  `gen_code!` with the signal directly).
- **Sub-chip oversampling required:** `sampling_frequency ≥
  code_frequency · subchip_factor` (else `gen_code!` raises an error).
- **Integer-chip phase only:** `start_phase` + the `start_index_shift`
  contribution are rounded to the nearest *primary* chip; the fractional
  sub-chip residual is dropped (up to ~1 sub-chip vs the plain `gen_code!`).
  `start_phase = 0.0, start_index_shift = 0` gives `phase = 0` exactly.
- **High-oversampling run-fill is approximate to ≤ a couple of samples.** Above
  ~8× sub-chip oversampling the resampler broadcast-fills runs of identical chips
  (matching the original `gen_code!`'s speed there) using a fixed-point samples-
  per-chip DDA. Its chip boundaries are exact for any single fill and drift at
  most a couple of samples only over a *very* long continued stream (~10⁶ chips) —
  the same order as the permute path's own rate-quantisation drift, and well
  inside the integer-chip-phase rounding above.
"""
struct CodeReplicaLUT{S<:AbstractGNSSSignal}
    signal::S
    prn::Int
    mc::CodeLUT.ModulatedCode   # rate-INDEPENDENT expanded ±1 table, built once
end

function CodeReplicaLUT(signal::AbstractGNSSSignal, prn::Integer;
                       cboc_amplitudes::Tuple{Integer,Integer} = (19, 6))
    modulation = _codelut_modulation(get_modulation(signal), cboc_amplitudes)
    Lp = get_code_length(signal)
    primary = Int8[get_code_at_index(signal, i, prn) for i in 0:Lp-1]
    sec = get_secondary_code(signal)
    Ls = secondary_code_length(sec)
    # secondary_value returns `true` for NoSecondaryCode (→ +1) and the stored
    # ±1 chip otherwise; Int8(true) == 1, so this maps both conventions to ±1.
    secondary = Ls == 1 ? Int8[1] : Int8[Int8(secondary_value(sec, prn, s)) for s in 0:Ls-1]
    # max_bake = typemax(Int16) keeps baked tables AVX2-addressable; longer
    # secondaries (e.g. L1C-P's 1800-chip overlay) stay unbaked.
    mc = CodeLUT.code_replica(primary, modulation; secondary = secondary, max_bake = typemax(Int16))
    CodeReplicaLUT{typeof(signal)}(signal, Int(prn), mc)
end

# ─────────────────────────────────────────────────────────────────────────────
# `CodeGeneratorLUT`: a mutable, *continuing* generator built once per (plan, rate).
#
# The DDA setup (fixed-point step + four stream inits, ~40 ns) runs ONCE in the
# `CodeGeneratorLUT` constructor. Each `gen_code!(out, gen)` then fills the NEXT `length(out)`
# samples from the carried state and saves the state advanced by exactly `length(out)`, so
# concatenating consecutive fills equals one big generation. This amortises the init across
# every 1 ms integration — the per-sample loop already beats the original `gen_code!`.
# ─────────────────────────────────────────────────────────────────────────────

"""
    CodeGeneratorLUT

Mutable, stateful, *continuing* LUT code generator built from a [`CodeReplicaLUT`](@ref)
plan and a sampling/code rate via its constructor (`CodeGeneratorLUT(plan, fs, fc)`). It
holds the resampler's DDA state (the `CodeLUT` rel/scalar-base streams, or the phase-based
equivalent for AVX2/Portable), the precomputed step ratio, the chosen backend, and a
secondary-period counter for signals with a non-baked secondary (e.g. GPS L5I's NH10).

Build it once (the DDA init runs here), then call [`gen_code!`](@ref)`(out, gen)` per
integration — that hot path does no rate setup and no DDA re-init, just the windowed
permute loop with a single-stream + scalar tail. It can also be iterated to yield `Vec{W,Int8}`
chunks (advancing the state) for fused, allocation-free correlation against a carrier.

Not thread-safe (it mutates its DDA state); use one generator per stream/channel.
"""
mutable struct CodeGeneratorLUT{S<:AbstractGNSSSignal,G<:CodeLUT.CodeGeneratorAny}
    const plan::CodeReplicaLUT{S}
    const engine::G                 # backend-specific continuing DDA generator
    const secondary::Vector{Int8}   # residual (non-baked) secondary; [1] if none
    const period_subchips::Int      # sub-chips per primary period
    const subchip_factor::Int
    const step_num::Int             # sub-chip step over the *sub-chip* table
    const step_den::Int
    const phase_sub::Int            # initial sub-chip phase offset
    n_abs::Int                      # absolute sample offset of the next sample to emit
end

"""
    CodeGeneratorLUT(plan::CodeReplicaLUT, sampling_frequency,
                     code_frequency = get_code_frequency(plan.signal);
                     start_phase = 0.0, start_index_shift = 0) -> CodeGeneratorLUT

Construct a continuing code generator over `plan` at the given rate. The one-time DDA setup
(fixed-point step + stream init) runs here; afterwards [`gen_code!`](@ref)`(out, gen)` fills
successive chunks with no per-call setup, and `for v in gen` yields `Vec{W,Int8}` chunks
(advancing the state). Use it for seamless block-to-block continuation (tracking) or
register-fused iteration; for a one-shot fill, the [`gen_code!`](@ref)`(out, plan, …)`
method is simpler.

Same oversampling requirement and integer-chip phase limitation as the plan
[`gen_code!`](@ref) method.
"""
function CodeGeneratorLUT(
    plan::CodeReplicaLUT,
    sampling_frequency,
    code_frequency = get_code_frequency(plan.signal);
    start_phase = 0.0,
    start_index_shift::Integer = 0,
)
    mc = plan.mc
    fc = _to_hz(code_frequency)
    fs = _to_hz(sampling_frequency)
    P = mc.subchip_factor
    fs < fc * P && error(
        "CodeGeneratorLUT needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz); use gen_code! with the signal directly.",
    )
    # Integer primary-chip phase (matches gen_code!'s start_phase_including_shift; the
    # fractional sub-chip residual is dropped — see the CodeReplicaLUT docstring).
    eff_chips = start_phase + start_index_shift * fc / fs
    phase = round(Int, eff_chips)
    # Parameterless default_backend() const-folds to the host's concrete backend; the
    # table-aware overload's length guard (fall back to Portable for tables larger than
    # typemax(Int32)) is dead code for GNSS codes and would erase that inference, boxing
    # the runtime-typed engine on every construction.
    backend = CodeLUT.default_backend()
    # Resample the baked sub-chip table at fc·P; phase is scaled to sub-chips.
    cps = (fc * P) / fs
    sn, sd = CodeLUT._fixed_point_step(cps)   # ← fixed-point step (one multiply + round)
    phase_sub = phase * P
    engine = CodeLUT.make_generator(mc.table, sn, sd; phase = phase_sub, backend = backend)
    return _wrap_code_generator_lut(plan, engine, mc.secondary, mc.period_subchips, P, sn, sd, phase_sub)
end

# Function barrier. `make_generator` returns a small Union over the phase type (Int16 vs
# Int32, chosen by _phase_type from the runtime table length), so the engine is Union-typed
# at the call site. Splitting on the concrete engine type `G` here keeps the
# CodeGeneratorLUT construction type-stable — without it the runtime-typed engine is boxed
# into the struct on every construction (the one-shot/threaded allocation hot spot).
function _wrap_code_generator_lut(
    plan::CodeReplicaLUT{S}, engine::G, secondary, period_subchips,
    subchip_factor, step_num, step_den, phase_sub,
) where {S,G<:CodeLUT.CodeGeneratorAny}
    CodeGeneratorLUT{S,G}(
        plan, engine, secondary, period_subchips, subchip_factor, step_num, step_den, phase_sub, 0,
    )
end

"""
    gen_code!(sampled_code, gen::CodeGeneratorLUT) -> sampled_code

Fill `sampled_code` with the **next** `length(sampled_code)` samples of the resampled
fully-modulated ±1 replica, continuing seamlessly from `gen`'s current state (no
`rationalize`, no DDA re-init — this is the hot path). The state is advanced by exactly
`length(sampled_code)`, so concatenating the outputs of consecutive calls equals one big
generation. Any non-baked secondary (e.g. GPS L5I's NH10) is applied per primary period
across call boundaries via the carried period counter.

`sampled_code` is most naturally `Int8` (written directly); other integer element types go
through a cached `Int8` scratch buffer + broadcast convert.
"""
function gen_code!(sampled_code::AbstractVector, gen::CodeGeneratorLUT)
    N = length(sampled_code)
    if eltype(sampled_code) == Int8
        CodeLUT.fill_continue!(sampled_code, gen.engine)
        _apply_secondary_continue!(sampled_code, gen)
    else
        scratch = get!(() -> Vector{Int8}(undef, N), _GEN_CODE_LUT_SCRATCH, N)
        CodeLUT.fill_continue!(scratch, gen.engine)
        _apply_secondary_continue!(scratch, gen)
        sampled_code .= scratch
    end
    gen.n_abs += N
    return sampled_code
end

# Apply the non-baked secondary as a per-primary-period sign flip over `out` (whose first
# sample is absolute sample `gen.n_abs`). Constant over a period → contiguous range negate.
# No-op when the secondary is baked / absent.
@inline function _apply_secondary_continue!(out::AbstractVector{<:Integer}, gen::CodeGeneratorLUT)
    sec = gen.secondary
    length(sec) <= 1 && return out
    Ls = length(sec); per = gen.period_subchips
    sn = gen.step_num; sd = gen.step_den
    n0 = gen.n_abs; N = length(out); ps = gen.phase_sub
    # Absolute sample n (0-based) maps to sub-chip floor(n·sn/sd) + ps; period p spans
    # sub-chips [p·per, (p+1)·per). First sample of period p: smallest n with sub-chip ≥ p·per.
    # The first emitted sample is absolute n0, so window index = n - n0.
    sub0 = (n0 * sn) ÷ sd + ps
    p = (sub0) ÷ per                # period index of the first emitted sample
    @inbounds while true
        T0 = p * per - ps
        n_start = T0 <= 0 ? 0 : cld(T0 * sd, sn)      # absolute first sample of period p
        n_start - n0 >= N && break
        T1 = (p + 1) * per - ps
        n_end = min(T1 <= 0 ? 0 : cld(T1 * sd, sn), n0 + N)  # absolute end (exclusive)
        if sec[mod(p, Ls) + 1] == -1
            lo = max(n_start, n0) - n0 + 1            # 1-based into out
            hi = n_end - n0                            # 1-based inclusive
            @simd for n in lo:hi
                out[n] = -out[n]
            end
        end
        p += 1
    end
    out
end

"""
    gen_code!(sampled_code, plan::CodeReplicaLUT, sampling_frequency,
              code_frequency = get_code_frequency(plan.signal),
              start_phase = 0.0, start_index_shift = 0, PHASET = Int32)

Convenience one-shot: fills `sampled_code` once from the requested integer-chip phase of the
given rate, paying a fresh DDA setup each call. Fine for one-off use; for repeated
integrations build `gen = CodeGeneratorLUT(plan, fs, fc)` once and call `gen_code!(out, gen)`
per integration to amortise the init.

Unlike the continuing [`CodeGeneratorLUT`](@ref), this drives the *immutable* one-shot
windowed fill (`CodeLUT.generate_code!`) directly — no mutable generator is built,
so the fill is allocation-free (the continuing generator's mutable DDA state is unnecessary
when filling exactly once).

`sampled_code` is most naturally `Int8` (written directly); other integer element types go
through a cached `Int8` scratch buffer + broadcast convert. Returns `sampled_code`.
"""
function gen_code!(
    sampled_code::AbstractVector,
    plan::CodeReplicaLUT,
    sampling_frequency,
    code_frequency = get_code_frequency(plan.signal),
    start_phase = 0.0,
    start_index_shift::Integer = 0,
    PHASET = Int32,
)
    mc = plan.mc
    fc = _to_hz(code_frequency)
    fs = _to_hz(sampling_frequency)
    P = mc.subchip_factor
    fs < fc * P && error(
        "CodeReplicaLUT gen_code! needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz); use gen_code! with the signal directly.",
    )
    # Integer primary-chip phase (matches CodeGeneratorLUT; the fractional sub-chip residual
    # is dropped — see the CodeReplicaLUT docstring).
    phase = round(Int, start_phase + start_index_shift * fc / fs)
    # Parameterless default_backend() const-folds to a concrete backend, keeping
    # generate_code! type-stable (the table-aware overload would erase that inference).
    backend = CodeLUT.default_backend()
    if eltype(sampled_code) == Int8
        CodeLUT.generate_code!(sampled_code, mc;
            code_frequency = fc, sampling_frequency = fs, phase = phase, backend = backend)
    else
        N = length(sampled_code)
        scratch = get!(() -> Vector{Int8}(undef, N), _GEN_CODE_LUT_SCRATCH, N)
        CodeLUT.generate_code!(scratch, mc;
            code_frequency = fc, sampling_frequency = fs, phase = phase, backend = backend)
        sampled_code .= scratch
    end
    return sampled_code
end

# ─────────────────────────────────────────────────────────────────────────────
# Value-based code engine: top-level adapter over `CodeLUT.code_engine` for a
# `CodeReplicaLUT` plan. Build one engine per interleave factor `K` (`Val(K)`), hold K
# isbits `code_state`s (W apart), and drive them with `code_lookup` / `code_advance` — the
# allocation-free, register-resident counterpart to filling an array with `gen_code!`, and
# the code-side partner for SinCosLUT's `carrier_engine`. (See CodeLUT's iterate.jl for the
# per-backend note: AVX-512 ≈ 8× the AVX2/NEON rate, an ISA limit, not a tuning miss.)
# ─────────────────────────────────────────────────────────────────────────────

using .CodeLUT: code_state, code_lookup, code_advance, code_width

"""
    code_engine(plan::CodeReplicaLUT, sampling_frequency,
                code_frequency = get_code_frequency(plan.signal), Val(K);
                start_phase = 0.0, start_index_shift = 0) -> CodeLUT.CodeEngine

Build a loop-invariant, value-based code engine over `plan` for a `K`-way interleaved
fused loop. Pair with `K` states `code_state(eng, k)` (`k = 0..K-1`, `W` samples apart) and
drive each with [`code_lookup`](@ref) / [`code_advance`](@ref); nothing is materialised or
heap-allocated. The code-side counterpart to SinCosLUT's `carrier_engine`. Same oversampling
requirement and integer-chip phase limitation as the plan [`gen_code!`](@ref) method; a
non-baked secondary (e.g. GPS L5I's NH10) or `sampling_frequency < code_frequency·subchip_factor`
raises an error (use [`gen_code!`](@ref) / [`CodeGeneratorLUT`](@ref) instead).
"""
function code_engine(plan::CodeReplicaLUT, sampling_frequency, code_frequency, ::Val{K};
                     start_phase = 0.0, start_index_shift::Integer = 0) where {K}
    mc = plan.mc
    fc = _to_hz(code_frequency)
    fs = _to_hz(sampling_frequency)
    P = mc.subchip_factor
    fs < fc * P && error(
        "code_engine needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz); use gen_code! / CodeGeneratorLUT.",
    )
    length(mc.secondary) > 1 && error(
        "code_engine does not support a non-baked secondary (e.g. GPS L5I NH10); use gen_code! / CodeGeneratorLUT.",
    )
    # Integer primary-chip phase (matches gen_code!'s start_phase_including_shift; the
    # fractional sub-chip residual is dropped — see the CodeReplicaLUT docstring).
    eff_chips = start_phase + start_index_shift * fc / fs
    phase = round(Int, eff_chips)
    # Resample the baked sub-chip table at fc·P, phase scaled to sub-chips. Parameterless
    # default_backend() const-folds to the host's concrete backend (keeps the engine type
    # inferable; the table-aware overload would box it).
    CodeLUT.code_engine(mc.table, (fc * P) / fs, Val(K);
                        phase = phase * P, backend = CodeLUT.default_backend())
end

code_engine(plan::CodeReplicaLUT, sampling_frequency, vk::Val; kwargs...) =
    code_engine(plan, sampling_frequency, get_code_frequency(plan.signal), vk; kwargs...)
