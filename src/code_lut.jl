# ─────────────────────────────────────────────────────────────────────────────
# Fast LUT-based code generation (THE `gen_code!` generator + continuing fill + value engine).
#
# This file vendors the GNSSSignalsLUT.jl machinery as an internal submodule
# `CodeLUT` (a verbatim copy of its permute / kernel / iterate / modulation
# sources plus its backend-selection + `CodeTable` definition) and adds a thin
# adapter at GNSSSignals top level: the `gen_code!(out::Vector{Int8}, signal, prn, …)`
# generator (resampling PRN `prn`'s baked column of the signal's embedded `SignalLUT`),
# plus the continuing immutable [`code_engine`](@ref)`(signal, prn, fs, fc)` →
# [`CodeFillEngine`](@ref) / [`code_state`](@ref) / value-threaded `gen_code!(out, eng, st)`,
# and the fused value-based [`code_engine`](@ref)`(signal, prn, fs, fc, Val(K))`. All read
# `signal.lut`; there is no rate-independent plan object — the baked table lives in the signal.
#
# The LUT resampler bakes the BOC/TMBOC subcarrier (and short secondary codes)
# into an expanded ±1 Int8 table and resamples it with a single AVX-512 `vpermb`
# / AVX2 `vpshufb` sliding-window permute over a drift-free integer DDA — or, once the
# baked table is heavily oversampled (so consecutive samples repeat a chip), an exact
# boundary fill that splat-stores one chip run per store at the original `gen_code!`'s
# store-bound speed instead of paying a permute per window. The baked table is Int8: ±1 for BPSK/BOC/TMBOC, or a multi-level
# integer approximation of the sqrt-power amplitudes for CBOC (Galileo E1B); cosine-BOC is
# unsupported. Requires sub-chip oversampling (`sampling_frequency ≥ code_frequency · subchip_factor`).
# ─────────────────────────────────────────────────────────────────────────────

# `GNSSSignals.CodeLUT` — internal submodule vendoring the GNSSSignalsLUT.jl SIMD code
# resampler. Not exported and not part of the public GNSSSignals API; use the public
# `gen_code!` / `gen_code` / `code_engine` instead. All names (`CodeTable`, `LOC`, `BOC`, `TMBOC`,
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
struct CodeTable{V<:AbstractVector{Int8}}
    chips::V    # length L
    padded::V   # length L + WINDOW_PAD
    length::Int
end

# vpermb reads a 64-chip window; the last valid base is L-1, reading up to index L+62
# (0-based). Pad by 63 so that load is always in-bounds.
const WINDOW_PAD = 63

# NOTE on zero-copy view tables: the parametric struct's auto-generated outer constructor
# `CodeTable(chips, padded, length)` infers `V` from the arguments, so passing unit-stride
# `SubArray{Int8}` column views (e.g. a column of `SignalLUT.padded`) yields a view-backed
# `CodeTable{<:SubArray}` with NO data copied — the SIMD `VecRange{W}` loads in the permute
# kernels read straight from the column (verified to lower to the same contiguous vload as a
# `Vector{Int8}`). Passing `Vector{Int8}` (the bake path below) yields the original owning table.

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
    # Select from the RUNTIME-detected feature set (populated in `__init__`), NOT the
    # precompile-time `HOST_FEATURES` const: a pkgimage baked on a VBMI host but loaded on a
    # non-VBMI CPU must demote AVX512 -> AVX2/Portable instead of emitting `vpermb` and SIGILL
    # (#104). This trades the old const-fold for a single Ref load per call — off the hot path
    # (called once per gen_code!/engine build, not per sample), and the returned value is a
    # small `Union{AVX512,AVX2,Portable}` the callers branch on statically.
    default_backend() = _select_backend(RUNTIME_FEATURES[])

    # Re-validate CPU features on the machine actually running, once per session. `HOST_FEATURES`
    # was baked at precompile and may over-claim on a relocated/shared pkgimage; refreshing the
    # runtime Ref here is what makes `default_backend()` demote away from AVX-512 on a non-VBMI
    # host (see permute.jl / #104).
    __init__() = _refresh_host_features!()
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
# Adapter: the public `gen_code!` / `gen_code` generators, the continuing immutable
# `CodeFillEngine` (value-threaded `CodeFillState`), and the value-based `code_engine` —
# all reading the signal's embedded `SignalLUT`. Int8-only output, no scratch Dict, no plan.
# ─────────────────────────────────────────────────────────────────────────────

# Strip a Unitful `Frequency` to a plain Float64 Hz value. The public entry points constrain
# `sampling_frequency`/`code_frequency` to `Frequency`, so a bare number raises a MethodError
# there rather than being silently misread as Hz (issue #105); no plain-number method here.
@inline _to_hz(x::Frequency) = Float64(ustrip(u"Hz", x))

# Split a real primary-chip start phase into the sub-chip integer offset `θ_int` and the
# fixed-point fractional sub-chip residual `rem0 ∈ [0, step_den)` the kernels consume.
#
#   eff_chips = start_phase + start_index_shift·fc/fs   (real primary chips)
#   θ         = eff_chips · P                            (real sub-chips)
#   θ_int     = floor(θ)                                 → kernel phase offset (does mod L)
#   rem0      = round((θ − θ_int)·step_den)              ∈ [0, step_den)
#
# `floor` (not round) is correct for negative phases too (e.g. start_index_shift = -1); the
# fractional part is always in [0,1). The rounding edge (rem0 == step_den) carries into θ_int.
# The sub-sample phase is tracked at 2^-30-sub-chip precision. rem0 = 0 ⇒ integer-sub-chip phase.
@inline function _subchip_phase_split(start_phase, start_index_shift, fc, fs, P, step_den::Int)
    θ = (start_phase + start_index_shift * fc / fs) * P
    θ_int = floor(Int, θ)
    rem0 = round(Int, (θ - θ_int) * step_den)
    if rem0 >= step_den
        rem0 -= step_den
        θ_int += 1
    end
    (θ_int, rem0)
end

# Map a GNSSSignals modulation to a CodeLUT modulation. ±1 for LOC/BOC/TMBOC; CBOC bakes a
# multi-level Int8 integer approximation of its sqrt-power amplitudes (`cboc_amplitudes`).
# Errors on cosine BOC and code-factor n ≠ 1. The `cboc_amplitudes` argument is ignored by
# every modulation except CBOC (the generic two-arg method below forwards to the one-arg form).
# Used by `build_signal_lut` (eager bake at signal construction).
_codelut_modulation(m, ::Tuple{Integer,Integer}) = _codelut_modulation(m)
function _codelut_modulation(m::LOC)
    CodeLUT.LOC()
end
function _codelut_modulation(m::BOCsin)
    m.n == 1 || error("the embedded LUT supports only code-factor n==1 BOC")
    CodeLUT.BOC(Int(m.m))
end
function _codelut_modulation(m::TMBOC)
    (m.boc1.n == 1 && m.boc2.n == 1 && m.boc1.m == 1) ||
        error("the embedded LUT supports only TMBOC with BOC(1,1) base and code-factor n==1")
    CodeLUT.TMBOC(Int(m.boc2.m), collect(Bool, m.pattern))
end
function _codelut_modulation(m::CBOC, cboc_amplitudes::Tuple{Integer,Integer})
    (m.boc1.n == 1 && m.boc2.n == 1) ||
        error("the embedded LUT supports only code-factor n==1 CBOC")
    a1, a2 = cboc_amplitudes
    (a1 > 0 && a2 > 0) ||
        error("cboc_amplitudes must be positive integers; got $cboc_amplitudes")
    Int(a1) + Int(a2) <= typemax(Int8) ||
        error("cboc_amplitudes must satisfy a1 + a2 ≤ $(typemax(Int8)) (Int8 table); got $cboc_amplitudes")
    # `boc2_sign` selects the relative phase of the BOC(m2,1) component: CBOC(+) for Galileo
    # E1B (`+1`), CBOC(−) for Galileo E1C (`-1`). CodeLUT.CBOC bakes `a1·BOC(m1,1) + a2·BOC(m2,1)`,
    # so a negative second amplitude reproduces the anti-phase composite.
    a2_signed = m.boc2_sign < 0 ? Int8(-Int(a2)) : Int8(a2)
    CodeLUT.CBOC(Int(m.boc1.m), Int(m.boc2.m), Int8(a1), a2_signed)
end
function _codelut_modulation(m::BOCcos)
    error("the embedded LUT does not support cosine-phased BOC (Int8/±1 only)")
end

"""
    build_signal_lut(modulation, codes, sec::SecondaryCode; cboc_amplitudes = (19, 6))
        -> SignalLUT

Bake every PRN's fully-modulated replica into one `SignalLUT` matrix, eagerly at signal
construction. For each PRN it bakes the padded table via
`CodeLUT.code_replica(primary, mod; secondary, max_bake = typemax(Int16))` and stacks the
per-PRN `mc.table.padded` as column `prn`. The shared metadata (P, table_length,
period_subchips) is read from the first PRN's `mc`.

`codes` is the signal's primary code matrix (`Int16`, ±1 per chip, unpadded); the primary
±1 for PRN `prn` is `sign(codes[:, prn])`. `sec` is the signal's [`SecondaryCode`](@ref):
short secondaries are BAKED into each column, a long overlay (GPS L1C-P, 1800 chips) is kept
RESIDUAL and stored in `SignalLUT.secondary` (per-PRN for a `PerPRNSecondaryCode`, one shared
column otherwise) for runtime application.

**Errors** for modulations the LUT can't bake (cosine BOC, code-factor n ≠ 1, …) — these are
unimplemented for the embedded LUT, so the signal fails to construct rather than silently
omitting its LUT. `_codelut_modulation` raises the specific error.
"""
function build_signal_lut(modulation, codes::AbstractMatrix, sec::SecondaryCode;
                          cboc_amplitudes::Tuple{Integer,Integer} = (19, 6))
    clmod = _codelut_modulation(modulation, cboc_amplitudes)   # throws on unimplemented modulation
    num_prns = size(codes, 2)
    num_prns >= 1 ||
        throw(ArgumentError("build_signal_lut: code matrix has no PRNs (got $(size(codes)))"))
    Lp = size(codes, 1)
    # First PRN: build the reference `mc`, size the padded matrix, read the shared metadata.
    sec1 = _residual_secondary_for_prn(sec, 1)
    mc1 = CodeLUT.code_replica(_primary_col(codes, 1), clmod; secondary = sec1, max_bake = typemax(Int16))
    coltot = length(mc1.table.padded)        # table_length + WINDOW_PAD
    padded = Matrix{Int8}(undef, coltot, num_prns)
    @inbounds padded[:, 1] .= mc1.table.padded
    # Residual secondary store: `mc1.secondary` is `Int8[1]` when baked/none. When residual and
    # PER-PRN (L1C-P), keep one column per PRN; otherwise one shared column reused for all PRNs.
    Lsr = length(mc1.secondary)
    residual_perprn = Lsr > 1 && sec isa PerPRNSecondaryCode
    secmat = Matrix{Int8}(undef, Lsr, residual_perprn ? num_prns : 1)
    @inbounds secmat[:, 1] .= mc1.secondary
    for prn in 2:num_prns
        secp = _residual_secondary_for_prn(sec, prn)
        mc = CodeLUT.code_replica(_primary_col(codes, prn), clmod; secondary = secp, max_bake = typemax(Int16))
        # All PRNs of one signal share P / table length, so the columns match.
        length(mc.table.padded) == coltot ||
            error("build_signal_lut: PRN $prn baked table length $(length(mc.table.padded)) ≠ $coltot")
        @inbounds padded[:, prn] .= mc.table.padded
        residual_perprn && (@inbounds secmat[:, prn] .= mc.secondary)
    end
    SignalLUT(padded, mc1.subchip_factor, secmat, mc1.table.length, mc1.period_subchips)
end

# Primary ±1 column for PRN `prn` from the (Int16) code matrix — sign of each stored chip.
@inline _primary_col(codes::AbstractMatrix, prn::Integer) =
    Int8[Int8(sign(codes[i, prn])) for i in 1:size(codes, 1)]

# Residual ±1 secondary vector for PRN `prn`: `Int8[1]` when there is no secondary, else the
# per-PRN/-shared ±1 chips.
function _residual_secondary_for_prn(sec::SecondaryCode, prn::Integer)
    Ls = secondary_code_length(sec)
    Ls == 1 ? Int8[1] : Int8[Int8(secondary_value(sec, prn, s)) for s in 0:Ls-1]
end

# Build a transient, view-backed `CodeLUT.ModulatedCode` over PRN `prn`'s column of `lut.padded`
# (ZERO-COPY: `chips`/`padded` are unit-stride `SubArray{Int8}` views into the matrix — the SIMD
# `VecRange{W}` loads read them directly). Reuses the shared metadata + the per-PRN residual
# secondary view. Fed to the shared resample core / continuing fill engine / value engine.
@inline function _modulated_code_view(lut::SignalLUT, prn::Int)
    coltot = size(lut.padded, 1)
    # Use matching `UnitRange` row slices so both views share the same `SubArray` type `V`
    # (the struct's single backing-vector type parameter); the auto outer constructor requires
    # identical `V` for `chips` and `padded`.
    padded_col = @view lut.padded[1:coltot, prn]
    chips_col = @view lut.padded[1:lut.table_length, prn]
    table = CodeLUT.CodeTable(chips_col, padded_col, lut.table_length)
    CodeLUT.ModulatedCode(table, lut.subchip_factor, _signal_lut_secondary(lut, prn), lut.period_subchips)
end

"""
    gen_code!(sampled_code::AbstractVector{Int8}, signal, prn, sampling_frequency,
              code_frequency = get_code_frequency(signal),
              start_phase = 0.0, start_index_shift = 0) -> sampled_code

Generate the sampled spreading code for PRN `prn` of `signal` in-place, by resampling PRN
`prn`'s fully-modulated baked replica from the signal's embedded `SignalLUT`. This is THE code
generator: it wraps PRN `prn`'s baked column in a zero-copy view-backed table and drives the
SIMD windowed-permute / boundary-fill resampler.

# Output
`sampled_code` must be `Vector{Int8}` (or any `AbstractVector{Int8}`). Non-CBOC signals are
±1; CBOC (Galileo E1B) is the multi-level Int8 *integer approximation* of its sqrt-power
subcarrier amplitudes (sign pattern matches the float spec). There is no scratch buffer and
the read-only `signal.lut` is shared, so this is thread-safe (write your own Int8 buffer).

# Phase
`start_phase` (primary chips) and the `start_index_shift` contribution are honored at full
sub-sample resolution: the phase is split into an integer sub-chip offset plus a
`2^-30`-sub-chip fixed-point fractional residual seeded into the DDA. `start_phase = 0.0,
start_index_shift = 0` gives phase `0` exactly; negative `start_index_shift` is handled
correctly. Any non-baked secondary (e.g. GPS L5I's NH10) is applied per primary period.

# Rate quantization
The chips-per-sample rate itself is quantized to the DDA's fixed `2^-30` sub-chip grid, so
within a fill the phase drifts from the exact requested rate by at most `N · 2^-31`
sub-chips over `N` samples. This is negligible for closed-loop tracking (≈ 1e-4 chips over
a 200k-sample epoch, and every call re-anchors to the exact float `start_phase`), but a
second-scale constant-rate fill (simulation / open-loop snapshot) accrues ~0.01–0.02 chips
per second at 50 MHz. Split constant-rate fills longer than ~1e7 samples into segments and
re-anchor each segment via `start_phase`; see the usage docs for an example.

# Requirements
`sampling_frequency ≥ code_frequency · subchip_factor` (else an error is raised). The embedded
`SignalLUT` is always present (signals whose modulation the LUT can't bake fail to construct),
so there is no missing-LUT case.

# Examples
```julia-repl
julia> using Unitful: MHz
julia> buffer = zeros(Int8, 4000)
julia> gen_code!(buffer, GPSL1CA(), 1, 4MHz)
```
"""
function gen_code!(
    sampled_code::AbstractVector{Int8},
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(signal),
    start_phase = 0.0,
    start_index_shift::Integer = 0,
)
    mc = _modulated_code_view(signal.lut, Int(prn))
    _gen_code_lut_core!(
        sampled_code, mc, _to_hz(sampling_frequency), _to_hz(code_frequency),
        start_phase, start_index_shift,
    )
end

# Shared one-shot resample core, driven from a rate-independent `CodeLUT.ModulatedCode`
# (baked table + metadata: subchip_factor P, residual secondary, period_subchips), whose
# table VIEWS into PRN `prn`'s column of the signal's embedded `SignalLUT.padded` matrix
# (zero-copy). Honors the `fs ≥ fc·P` oversampling check, fractional sub-chip `start_phase` /
# `start_index_shift` via `_subchip_phase_split`, host backend selection, and the non-baked
# secondary application (carried inside `CodeLUT.generate_code!(out, mc; …)`). Int8 only.
@inline function _gen_code_lut_core!(
    sampled_code::AbstractVector{Int8},
    mc::CodeLUT.ModulatedCode,
    fs::Float64,
    fc::Float64,
    start_phase,
    start_index_shift::Integer,
)
    P = mc.subchip_factor
    fs < fc * P && error(
        "gen_code! needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz).",
    )
    # Split the real primary-chip start phase into an integer sub-chip offset (θ_int) and a
    # fractional residual `rem0` at 2^-30-sub-chip precision (see `_subchip_phase_split`).
    sd = CodeLUT._STEP_DEN
    phase_sub, rem0 = _subchip_phase_split(start_phase, start_index_shift, fc, fs, P, sd)
    # Parameterless default_backend() returns the host's concrete backend as a small
    # Union{AVX512,AVX2,Portable} (a runtime Ref read since #104, so no longer const-folded, but
    # union-split at the barrier below); the table-aware overload would erase even that.
    backend = CodeLUT.default_backend()
    CodeLUT.generate_code!(sampled_code, mc;
        code_frequency = fc, sampling_frequency = fs, phase_sub = phase_sub, rem0 = rem0, backend = backend)
    return sampled_code
end

# ─────────────────────────────────────────────────────────────────────────────
# Value-based continuing fill: an immutable `CodeFillEngine` built once per (signal, prn, rate)
# plus an isbits `CodeFillState` threaded by the caller. The DDA setup (fixed-point step +
# stream init, ~40 ns) runs ONCE in `code_engine`; each `gen_code!(out, eng, st)` then fills
# the NEXT `length(out)` samples from `st` and *returns* the state advanced by exactly
# `length(out)`, so concatenating consecutive fills equals one big generation. No hidden /
# mutable state — the caller threads `st` — and the isbits state keeps the steady-state fill
# allocation-free. This amortises the init across every 1 ms integration; for fused,
# register-resident correlation use the value-based [`code_engine`](@ref)`(signal, prn, fs, fc,
# Val(K))` / [`code_lookup`](@ref) API instead.
# ─────────────────────────────────────────────────────────────────────────────

"""
    CodeFillEngine

Immutable, loop-invariant LUT code-fill engine built from a `(signal, prn)` (reading
`signal.lut`) and a sampling/code rate via [`code_engine`](@ref)`(signal, prn, fs, fc)`. It
holds the chosen backend's resampler config (precomputed step ratio + DDA deltas), the
residual (non-baked) secondary (e.g. GPS L5I's NH10), and the initial phase used to seed
[`code_state`](@ref).

Build it once (the DDA setup runs here), pair it with a [`code_state`](@ref) seed, then call
[`gen_code!`](@ref)`(out, eng, st)` per integration — that hot path does no rate setup and no
DDA re-init, just the windowed permute loop with a single-stream + scalar tail, and returns
the advanced state. The engine is read-only and shareable; per-stream/channel position lives
entirely in the threaded state, so one engine can drive many independent streams.
"""
struct CodeFillEngine{E<:CodeLUT.FillEngineAny}
    engine::E                 # backend-specific value-based fill engine
    secondary::Vector{Int8}   # residual (non-baked) secondary; [1] if none
    period_subchips::Int      # sub-chips per primary period
    subchip_factor::Int
    step_num::Int             # sub-chip step over the *sub-chip* table
    step_den::Int
    phase_sub::Int            # initial sub-chip phase offset
    rem0::Int                 # initial fractional sub-chip residual (for the secondary boundaries)
end

"""
    CodeFillState

Immutable, isbits position state for a [`CodeFillEngine`](@ref): the backend DDA state plus
the absolute sample offset of the next sample to emit (used to place a non-baked secondary's
sign flips across call boundaries). Seed it with [`code_state`](@ref)`(eng)` and thread the
value returned by [`gen_code!`](@ref)`(out, eng, st)` into the next call.
"""
struct CodeFillState{St}
    dda::St          # backend DDA state (CodeState512 / CodeStatePhase / BoundaryState)
    n_abs::Int       # absolute sample offset of the next sample to emit
end

"""
    code_engine(signal::AbstractGNSSSignal, prn::Integer, sampling_frequency,
                code_frequency = get_code_frequency(signal);
                start_phase = 0.0, start_index_shift = 0,
                backend = default_backend()) -> CodeFillEngine

Build an immutable, continuing code-fill engine over PRN `prn` of `signal` at the given rate
(reading `signal.lut`; no `Val(K)` — that overload builds the fused register-resident engine
instead). The one-time DDA setup (fixed-point step + stream init) runs here; afterwards seed a
[`code_state`](@ref)`(eng)` and call [`gen_code!`](@ref)`(out, eng, st) -> st` to fill
successive chunks with no per-call setup, threading the returned state for seamless
block-to-block continuation (tracking). Non-baked secondaries (e.g. GPS L5I's NH10) are
supported. For a one-shot fill the [`gen_code!`](@ref)`(out, signal, prn, …)` method is simpler.

Same oversampling requirement and fractional sub-chip phase support as the one-shot
[`gen_code!`](@ref) — including its rate quantization: the engine locks the rate to the
DDA's `2^-30` sub-chip grid at construction and never re-anchors to the requested float
rate, so the phase drifts by at most `N · 2^-31` sub-chips over the `N` samples generated
since the engine's `start_phase` anchor. Closed-loop tracking is unaffected (the engine is
rebuilt — re-anchoring the phase — on every Doppler update); for open-loop constant-rate
streaming, rebuild the engine with a freshly computed `start_phase` about every 1e7
samples. `backend` defaults to the host's best SIMD backend; pass a weaker one
(e.g. `CodeLUT.AVX2()`/`CodeLUT.Portable()`) to force it — handy for testing the non-default
paths on a given CPU. Forcing a backend the CPU does not support is invalid.
"""
function code_engine(
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_code_frequency(signal);
    start_phase = 0.0,
    start_index_shift::Integer = 0,
    backend::CodeLUT.Backend = CodeLUT.default_backend(),
)
    mc = _modulated_code_view(signal.lut, Int(prn))
    fc = _to_hz(code_frequency)
    fs = _to_hz(sampling_frequency)
    P = mc.subchip_factor
    fs < fc * P && error(
        "code_engine needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz); use gen_code! with the signal directly.",
    )
    # `backend` defaults to the parameterless `default_backend()`, a small
    # `Union{AVX512,AVX2,Portable}` (a runtime Ref read since #104) union-split at the
    # `_wrap_code_fill_engine` barrier — the table-aware overload's length guard is dead code
    # for GNSS codes and would additionally box the engine. It is exposed so tests can force a
    # *weaker* backend than the host (e.g. AVX2/Portable on an AVX-512 runner) — forcing a
    # backend the CPU lacks would SIGILL.
    # Resample the baked sub-chip table at fc·P. Split the real primary-chip start phase into
    # an integer sub-chip offset `phase_sub` (θ_int) and a fractional residual `rem0`.
    cps = (fc * P) / fs
    sn, sd = CodeLUT._fixed_point_step(cps)   # ← fixed-point step (one multiply + round)
    phase_sub, rem0 = _subchip_phase_split(start_phase, start_index_shift, fc, fs, P, sd)
    engine = CodeLUT.make_fill_engine(mc.table, sn, sd; phase = phase_sub, rem0 = rem0, backend = backend)
    # `mc.secondary` is a view into `signal.lut.secondary`; materialise it as an owning Vector
    # so the engine does not alias the (otherwise read-only) embedded matrix.
    return _wrap_code_fill_engine(engine, collect(mc.secondary), mc.period_subchips, P, sn, sd, phase_sub, Int(rem0))
end

# Function barrier. `make_fill_engine` returns a small Union over the phase type (Int16 vs
# Int32, chosen by _phase_type from the runtime table length), so the engine is Union-typed
# at the call site. Splitting on the concrete engine type `E` here keeps the CodeFillEngine
# construction type-stable — without it the runtime-typed engine is boxed into the struct on
# every construction (the one-shot/threaded allocation hot spot).
function _wrap_code_fill_engine(
    engine::E, secondary, period_subchips,
    subchip_factor, step_num, step_den, phase_sub, rem0,
) where {E<:CodeLUT.FillEngineAny}
    CodeFillEngine{E}(
        engine, secondary, period_subchips, subchip_factor, step_num, step_den, phase_sub, rem0,
    )
end

"""
    gen_code!(sampled_code::AbstractVector{Int8}, eng::CodeFillEngine, st::CodeFillState)
        -> CodeFillState

Fill `sampled_code` with the **next** `length(sampled_code)` samples of the resampled
fully-modulated replica, continuing seamlessly from `st` (no `rationalize`, no DDA re-init —
this is the hot path), and **return** the state advanced by exactly `length(sampled_code)`.
Concatenating the outputs of consecutive calls — threading the returned state — equals one big
generation. Any non-baked secondary (e.g. GPS L5I's NH10) is applied per primary period across
call boundaries via the state's absolute sample offset. Int8 output only.

Because the engine's rate is fixed on the DDA's `2^-30` sub-chip grid, the phase drift bound
of `N · 2^-31` sub-chips (see [`code_engine`](@ref)) applies to the *total* `N` samples
generated since the engine's phase anchor, not per call.
"""
function gen_code!(sampled_code::AbstractVector{Int8}, eng::CodeFillEngine, st::CodeFillState)
    N = length(sampled_code)
    dda = CodeLUT.fill_continue!(sampled_code, eng.engine, st.dda)
    _apply_secondary_continue!(sampled_code, eng, st.n_abs)
    return CodeFillState(dda, st.n_abs + N)
end

# Apply the non-baked secondary as a per-primary-period sign flip over `out` (whose first
# sample is absolute sample `n0`). Constant over a period → contiguous range negate.
# No-op when the secondary is baked / absent.
@inline function _apply_secondary_continue!(out::AbstractVector{<:Integer}, eng::CodeFillEngine, n0::Int)
    sec = eng.secondary
    length(sec) <= 1 && return out
    # Same period-walk as the one-shot fill, generalised by the absolute sample offset `n0`
    # (see `CodeLUT._apply_secondary!`); n0 = 0 reproduces the one-shot version exactly.
    CodeLUT._apply_secondary!(out, sec, eng.period_subchips, eng.step_num, eng.step_den,
                              eng.phase_sub, eng.rem0, n0)
end

# ─────────────────────────────────────────────────────────────────────────────
# Value-based code engine: top-level adapter over `CodeLUT.code_engine` for PRN `prn` of a
# signal (reading `signal.lut`). Build one engine per interleave factor `K` (`Val(K)`), hold K
# isbits `code_state`s (W apart), and drive them with `code_lookup` / `code_advance` — the
# allocation-free, register-resident counterpart to filling an array with `gen_code!`, and the
# code-side partner for SinCosLUT's `carrier_engine`. (See CodeLUT's iterate.jl for the
# per-backend note: AVX-512 ≈ 8× the AVX2/NEON rate, an ISA limit, not a tuning miss.)
# ─────────────────────────────────────────────────────────────────────────────

using .CodeLUT: code_state, code_lookup, code_advance, code_width

"""
    code_state(eng::CodeFillEngine) -> CodeFillState

Initial isbits fill state for `eng` (DDA seeded at the engine's start phase, absolute offset
0). Thread the value returned by [`gen_code!`](@ref)`(out, eng, st)` into the next call.
"""
@inline CodeLUT.code_state(eng::CodeFillEngine) = CodeFillState(CodeLUT.fill_state(eng.engine), 0)

"""
    code_engine(signal::AbstractGNSSSignal, prn::Integer, sampling_frequency,
                code_frequency = get_code_frequency(signal), Val(K);
                start_phase = 0.0, start_index_shift = 0) -> CodeLUT.CodeEngine

Build a loop-invariant, value-based code engine over PRN `prn` of `signal` (reading
`signal.lut`) for a `K`-way interleaved fused loop. Pair with `K` states `code_state(eng, k)`
(`k = 0..K-1`, `W` samples apart) and drive each with [`code_lookup`](@ref) /
[`code_advance`](@ref); nothing is materialised or heap-allocated. The code-side counterpart
to SinCosLUT's `carrier_engine`. Same oversampling requirement, fractional sub-chip phase
support, and `2^-30` rate quantization (drift ≤ `N · 2^-31` sub-chips over `N` samples) as
[`gen_code!`](@ref); a non-baked secondary (e.g. GPS L5I's NH10) or
`sampling_frequency < code_frequency·subchip_factor` raises an error (use the array-filling
[`gen_code!`](@ref) / continuing [`code_engine`](@ref) without `Val(K)` instead).
"""
function code_engine(signal::AbstractGNSSSignal, prn::Integer, sampling_frequency::Frequency, code_frequency::Frequency, ::Val{K};
                     start_phase = 0.0, start_index_shift::Integer = 0) where {K}
    mc = _modulated_code_view(signal.lut, Int(prn))
    fc = _to_hz(code_frequency)
    fs = _to_hz(sampling_frequency)
    P = mc.subchip_factor
    fs < fc * P && error(
        "code_engine(…, Val(K)) needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz); use gen_code! / code_engine without Val(K).",
    )
    any(!=(Int8(1)), mc.secondary) && error(
        "code_engine(…, Val(K)) does not support a non-baked secondary (e.g. GPS L5I NH10); use gen_code! / code_engine without Val(K).",
    )
    # Split the real primary-chip start phase into an integer sub-chip offset (θ_int) and a
    # fractional residual `rem0` (see `_subchip_phase_split`).
    sd = CodeLUT._STEP_DEN
    phase_sub, rem0 = _subchip_phase_split(start_phase, start_index_shift, fc, fs, P, sd)
    # Resample the baked sub-chip table at fc·P. Parameterless default_backend() returns the
    # host's concrete backend as a small Union (a runtime Ref read since #104), union-split into
    # the inferable per-backend engine; the table-aware overload would box it instead.
    CodeLUT.code_engine(mc.table, (fc * P) / fs, Val(K);
                        phase = phase_sub, rem0 = rem0, backend = CodeLUT.default_backend())
end

code_engine(signal::AbstractGNSSSignal, prn::Integer, sampling_frequency::Frequency, vk::Val; kwargs...) =
    code_engine(signal, prn, sampling_frequency, get_code_frequency(signal), vk; kwargs...)
