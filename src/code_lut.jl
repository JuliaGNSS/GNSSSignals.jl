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
# / AVX2 `vpshufb` sliding-window permute over a drift-free integer DDA. It is
# Int8/±1 only (no CBOC / cosine-BOC) and requires sub-chip oversampling
# (`sampling_frequency ≥ code_frequency · subchip_factor`).
# ─────────────────────────────────────────────────────────────────────────────

"""
    GNSSSignals.CodeLUT

Internal submodule vendoring the GNSSSignalsLUT.jl SIMD code resampler. Not
exported and not part of the public GNSSSignals API; use [`CodeReplicaLUT`](@ref)
with [`gen_code!`](@ref) / [`gen_code`](@ref) instead. All names (`CodeTable`, `LOC`, `BOC`, `TMBOC`,
`ModulatedCode`, `code_replica`, `generate_code!`, `generate_code`, `AVX512`,
`AVX2`, `Portable`, `default_backend`, …) live inside this submodule so they do
not clash with GNSSSignals' own `LOC` / `BOC` / `TMBOC` / `Modulation`.
"""
module CodeLUT

using SIMD

# ---- backends ----
abstract type Backend end
struct AVX512   <: Backend end   # vpermb over a 64-chip sliding window (W = 64)
struct AVX2     <: Backend end   # vpshufb over two independent 16-chip windows (W = 32)
struct Portable <: Backend end   # scalar fallback (any CPU)

backend_name(::AVX512)   = "AVX-512"
backend_name(::AVX2)     = "AVX2"
backend_name(::Portable) = "portable"

# Window width / SIMD lane count per backend.
_vwidth(::AVX512)   = Val(64)
_vwidth(::AVX2)     = Val(32)
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
else
    default_backend() = Portable()
end
# Length-aware: AVX2 widens its phase vector to Int32 for tables > typemax(Int16), so it
# addresses anything up to typemax(Int32) (slower than the Int16 path, but ~20× over
# scalar). Only fall back to Portable for the (unreachable for GNSS) even-longer tables.
function default_backend(table::CodeTable)
    be = default_backend()
    (be isa AVX2 && table.length > typemax(Int32)) ? Portable() : be
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

# Map a GNSSSignals modulation to a CodeLUT modulation (Int8/±1 only). Errors on
# unsupported modulations (CBOC, cosine BOC, code-factor n ≠ 1).
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
function _codelut_modulation(m::CBOC)
    error("CodeReplicaLUT does not support CBOC (Int8/±1 only); use gen_code!")
end
function _codelut_modulation(m::BOCcos)
    error("CodeReplicaLUT does not support cosine-phased BOC (Int8/±1 only); use gen_code!")
end

"""
    CodeReplicaLUT(signal::AbstractGNSSSignal, prn::Integer)

A reusable, rate-independent plan for fast LUT-based code generation of
`signal`/`prn`. Building it bakes the fully-modulated ±1 replica (primary ×
subcarrier × short secondary) into an expanded `Int8` `CodeTable` once — this is
the expensive step. The resulting plan is **immutable and read-only**, so it can
be shared across threads (build once per `(signal, prn)`, reuse it for every
integration; a receiver holds one plan per tracked channel).

Pass the plan to [`gen_code!`](@ref) or [`gen_code`](@ref) with a sampling
frequency to resample it. The baked table is independent of the sampling and code
frequency, so a single plan serves any rate.

# Limitations (vs the plain-`signal` `gen_code!`)
- **Int8 / ±1 only.** CBOC and cosine-phased BOC raise an error at construction
  (use `gen_code!` with the signal directly).
- **Sub-chip oversampling required:** `sampling_frequency ≥
  code_frequency · subchip_factor` (else `gen_code!` raises an error).
- **Integer-chip phase only:** `start_phase` + the `start_index_shift`
  contribution are rounded to the nearest *primary* chip; the fractional
  sub-chip residual is dropped (up to ~1 sub-chip vs the plain `gen_code!`).
  `start_phase = 0.0, start_index_shift = 0` gives `phase = 0` exactly.
"""
struct CodeReplicaLUT{S<:AbstractGNSSSignal}
    signal::S
    prn::Int
    mc::CodeLUT.ModulatedCode   # rate-INDEPENDENT expanded ±1 table, built once
end

function CodeReplicaLUT(signal::AbstractGNSSSignal, prn::Integer)
    modulation = _codelut_modulation(get_modulation(signal))
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
# The expensive DDA setup (rationalize + four `_init_rel` streams, ~550 ns) runs ONCE in
# the `gen_code` constructor. Each `gen_code!(out, gen)` then fills the NEXT `length(out)`
# samples from the carried state and saves the state advanced by exactly `length(out)`, so
# concatenating consecutive fills equals one big generation. This amortises the init across
# every 1 ms integration — the per-sample loop already beats the original `gen_code!`.
# ─────────────────────────────────────────────────────────────────────────────

"""
    CodeGeneratorLUT

Mutable, stateful, *continuing* LUT code generator built from a [`CodeReplicaLUT`](@ref)
plan and a sampling/code rate via [`gen_code`](@ref). It holds the resampler's DDA state
(the `CodeLUT` rel/scalar-base streams, or the phase-based equivalent for AVX2/Portable),
the precomputed step ratio, the chosen backend, and a secondary-period counter for signals
with a non-baked secondary (e.g. GPS L5I's NH10).

Build it once (the DDA init runs here), then call [`gen_code!`](@ref)`(out, gen)` per
integration — that hot path does no `rationalize` and no DDA re-init, just the windowed
permute loop with a fully vectorised tail. It can also be iterated to yield `Vec{W,Int8}`
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
    gen_code(plan::CodeReplicaLUT, sampling_frequency,
             code_frequency = get_code_frequency(plan.signal);
             start_phase = 0.0, start_index_shift = 0) -> CodeGeneratorLUT

Construct a [`CodeGeneratorLUT`](@ref): a continuing code generator over `plan` at the
given rate. **This is where the one-time DDA setup runs** (`rationalize` + stream init);
afterwards [`gen_code!`](@ref)`(out, gen)` fills successive chunks with no per-call setup.

Same oversampling requirement and integer-chip phase limitation as the plan
[`gen_code!`](@ref) method. The fast repeated path is: build the generator once outside
your integration loop, then `gen_code!(out, gen)` per integration.
"""
function gen_code(
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
        "gen_code (CodeReplicaLUT) needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz); use gen_code! with the signal directly.",
    )
    # Integer primary-chip phase (matches gen_code!'s start_phase_including_shift; the
    # fractional sub-chip residual is dropped — see the CodeReplicaLUT docstring).
    eff_chips = start_phase + start_index_shift * fc / fs
    phase = round(Int, eff_chips)
    backend = CodeLUT.default_backend(mc.table)
    # Resample the baked sub-chip table at fc·P; phase is scaled to sub-chips.
    cps = (fc * P) / fs
    sn, sd = CodeLUT._fixed_point_step(cps)   # ← fixed-point step (one multiply + round)
    phase_sub = phase * P
    engine = CodeLUT.make_generator(mc.table, sn, sd; phase = phase_sub, backend = backend)
    CodeGeneratorLUT{typeof(plan.signal),typeof(engine)}(
        plan, engine, mc.secondary, mc.period_subchips, P, sn, sd, phase_sub, 0,
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

Convenience one-shot: builds a [`CodeGeneratorLUT`](@ref) (paying the DDA init) and fills
`sampled_code` once from phase zero of the requested rate. Fine for one-off use; for
repeated integrations build `gen = gen_code(plan, fs, fc)` once and call
`gen_code!(out, gen)` per integration to amortise the init.

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
    gen = gen_code(plan, sampling_frequency, code_frequency;
                   start_phase = start_phase, start_index_shift = start_index_shift)
    gen_code!(sampled_code, gen)
end

# ---- iteration: yield Vec{W,Int8} chunks, advancing the generator state ----
CodeLUT.gen_width(gen::CodeGeneratorLUT) = CodeLUT.gen_width(gen.engine)

Base.IteratorSize(::Type{<:CodeGeneratorLUT}) = Base.SizeUnknown()
Base.eltype(::Type{CodeGeneratorLUT{S,G}}) where {S,G} = eltype(CodeLUT.GeneratorChunks{G})

@inline function Base.iterate(gen::CodeGeneratorLUT, ::Nothing = nothing)
    length(gen.secondary) > 1 &&
        error("iterating a CodeGeneratorLUT with a non-baked secondary is not supported; use gen_code!")
    W = CodeLUT.gen_width(gen.engine)
    chunks = CodeLUT.GeneratorChunks(gen.engine, typemax(Int))
    r = iterate(chunks)
    r === nothing && return nothing
    vec, _ = r
    gen.n_abs += W
    (vec, nothing)
end
