# ─────────────────────────────────────────────────────────────────────────────
# Fast LUT-based code generation: `gen_code!`, the continuing fill engine, and the value engine.
#
# Vendors the GNSSSignalsLUT.jl machinery as the internal submodule `CodeLUT` and adds a thin
# adapter at GNSSSignals top level: `gen_code!(out::Vector{Int8}, signal, prn, …)` (resamples PRN
# `prn`'s baked column of the signal's embedded `SignalLUT`), the continuing immutable
# [`code_engine`](@ref)`(signal, prn, fs, fc)` → [`CodeFillEngine`](@ref) / [`code_state`](@ref) /
# value-threaded `gen_code!(out, eng, st)`, and the fused value-based
# [`code_engine`](@ref)`(signal, prn, fs, fc, Val(K))`. All read `signal.lut`.
#
# The resampler bakes the BOC/TMBOC subcarrier (and short secondary codes) into an expanded Int8
# table and resamples it with a single AVX-512 `vpermb` / AVX2 `vpshufb` sliding-window permute over
# a drift-free integer DDA — or, once the baked table is heavily oversampled (consecutive samples
# repeat a chip), an exact boundary fill that splat-stores one chip run per store instead of
# paying a permute per window. The table is ±1 for BPSK/BOC/TMBOC, or a multi-level integer
# approximation of the sqrt-power amplitudes for CBOC (Galileo E1B). Sine BOC(m,n) is baked at
# 2m/n sub-chips per chip, cosine BOC(m,n) at 4m/n quarter-sub-chips (the quarter-period shift
# halves the grid); a BOC whose count is not a whole number cannot be baked. Requires sub-chip
# oversampling (`sampling_frequency ≥ code_frequency · subchip_factor`).
# ─────────────────────────────────────────────────────────────────────────────

include("code_lut/CodeLUT.jl")

# ─────────────────────────────────────────────────────────────────────────────
# Adapter: the public `gen_code!` / `gen_code` generators, the continuing immutable
# `CodeFillEngine` (value-threaded `CodeFillState`), and the value-based `code_engine`.
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# `SignalLUT`: per-signal embedded LUT holding all PRNs' baked code tables as one matrix.
#
# Defined before the signal structs that embed it (`lut::SignalLUT`); built and consumed by
# this file (`build_signal_lut`, `gen_code!`). The metadata is identical across a
# signal's PRNs, so one `padded` matrix (column `prn` = that PRN's baked table) plus the
# shared metadata is stored. `gen_code!` resamples it zero-copy.
# ─────────────────────────────────────────────────────────────────────────────
struct SignalLUT
    padded::Matrix{Int8}      # (table_length + WINDOW_PAD) × num_prns; column prn = that PRN's padded baked table
    subchip_factor::Int       # P
    # Residual non-baked secondary applied per primary period, as a (Ls × n) matrix; a single
    # row of 1 when none/baked. A long per-PRN overlay (GPS L1C-P's 1800 chips) is stored per
    # column and applied at runtime, since baking it would grow the matrix ~1800×. A shared
    # secondary (GPS L5I's NH10) stores one column reused for all PRNs.
    secondary::Matrix{Int8}
    table_length::Int         # L·P (length of one PRN's baked table, i.e. column length minus WINDOW_PAD)
    period_subchips::Int      # sub-chips per primary period (Lp·P), for secondary application
    # RMS amplitude of the baked code, `sqrt(mean(code^2))` over one primary table. Exactly
    # `1` for ±1 codes (BPSK/BOC/TMBOC); for the multi-level CBOC Int8 approximation it is
    # `sqrt(a1^2 + a2^2)` (≈ 19.92 for E1B's baked (19, 6)) and tracks any `cboc_amplitudes`
    # override. PRN-independent. Read via `get_code_amplitude`; divide correlation results
    # by it to get a unit-power scale.
    code_amplitude::Float64
end

# Residual secondary column for PRN `prn` (a shared secondary has one column for all PRNs).
@inline _signal_lut_secondary(lut::SignalLUT, prn::Int) =
    @view lut.secondary[:, size(lut.secondary, 2) == 1 ? 1 : prn]

# Strip a Unitful `Frequency` to Float64 Hz. The entry points constrain frequencies to `Frequency`,
# so a bare number raises a MethodError instead of being silently misread as Hz (#105).
@inline _to_hz(x::Frequency) = Float64(ustrip(u"Hz", x))

# Split a real primary-chip start phase into the sub-chip integer offset `θ_int` and the
# fixed-point fractional sub-chip residual `rem0 ∈ [0, step_den)` the kernels consume.
#
#   eff_chips = start_phase + start_index_shift·fc/fs   (real primary chips)
#   θ         = eff_chips · P                            (real sub-chips)
#   θ_int     = floor(θ)                                 → kernel phase offset (does mod L)
#   rem0      = round((θ − θ_int)·step_den)              ∈ [0, step_den)
#
# `floor` (not round) keeps the fractional part in [0,1) for negative phases too; the rounding
# edge (rem0 == step_den) carries into θ_int. Sub-sample phase has 2^-30-sub-chip precision.
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

# Shared preamble of the three LUT resample paths (`_gen_code_lut_core!`, `code_engine`,
# `code_engine(…, Val(K))`): applies the `fs ≥ fc·P` guard and splits the start phase via
# `_subchip_phase_split`. Sharing it keeps the guard and phase convention identical, so gen_code!
# and code_engine stay byte-for-byte equivalent (#131). `what` names the caller in the error and
# `hint` is an optional trailing redirect. (`_STEP_DEN` is the DDA denominator every path's
# `_fixed_point_step` also returns.)
@inline function _lut_resample_setup(
    mc::CodeLUT.ModulatedCode, fs::Float64, fc::Float64, start_phase, start_index_shift::Integer;
    what::AbstractString, hint::AbstractString = "",
)
    P = mc.subchip_factor
    fs < fc * P && error(
        "$what needs sampling_frequency ≥ code_frequency·subchip_factor (=$(fc * P) Hz)$hint.",
    )
    _subchip_phase_split(start_phase, start_index_shift, fc, fs, P, CodeLUT._STEP_DEN)
end

# Map a GNSSSignals modulation to a CodeLUT modulation. ±1 for BPSK/BOCsin/BOCcos/TMBOC; CBOC
# bakes a multi-level Int8 approximation of its sqrt-power amplitudes, derived from the CBOC's
# `boc1_power` (see `_cboc_int_amplitudes`) unless `cboc_amplitudes` is an explicit override
# (ignored by every other modulation). ArgumentError where the sub-chip count is not a whole
# number (see the file header). Used by `build_signal_lut`.
_codelut_modulation(m, ::Union{Nothing,Tuple{Integer,Integer}}) = _codelut_modulation(m)
function _codelut_modulation(m::BPSK)
    CodeLUT.BPSK()
end
_codelut_modulation(m::BOCsin) = CodeLUT.BOC(m.m, m.n)
# The composites bake their components' sine subcarriers at a common resolution; a cosine
# component would need the quarter grid and is not supported.
function _sine_components(m, what)
    m.boc1 isa BOCsin && m.boc2 isa BOCsin || throw(ArgumentError(
        "the code lookup table bakes a $what of sine-phased BOC components only; got " *
        "$(nameof(typeof(m.boc1))) and $(nameof(typeof(m.boc2)))"))
end
function _codelut_modulation(m::TMBOC)
    _sine_components(m, "TMBOC")
    CodeLUT.TMBOC(m.boc1.m, m.boc2.m, m.boc1.n, collect(Bool, m.pattern))
end
# Best integer amplitude pair (a1, a2) for baking a CBOC subcarrier: the continued-fraction
# convergent of the target amplitude ratio √(boc1_power) : √(1-boc1_power) whose sum stays within
# the Int8 table budget a1 + a2 ≤ 127. The E1 default (boc1_power = 10/11, ratio √10) yields (19, 6).
function _cboc_int_amplitudes(boc1_power::Real)
    ratio = sqrt(boc1_power / (1 - boc1_power))          # target a1 : a2
    h2, h1 = 0, 1                                         # convergent numerators   (h_{-2}, h_{-1})
    k2, k1 = 1, 0                                         # convergent denominators (k_{-2}, k_{-1})
    x = float(ratio)
    best = (1, 1)
    for _ in 1:64                                         # depth cap; the budget breaks far sooner
        a = floor(Int, x)
        h0 = a * h1 + h2
        k0 = a * k1 + k2
        h0 + k0 <= typemax(Int8) || break
        (h0 >= 1 && k0 >= 1) && (best = (h0, k0))
        h2, h1 = h1, h0
        k2, k1 = k1, k0
        frac = x - a
        frac > 1e-9 || break                             # exhausted the expansion
        x = 1 / frac
    end
    best
end

function _codelut_modulation(m::CBOC,
                             cboc_amplitudes::Union{Nothing,Tuple{Integer,Integer}} = nothing)
    _sine_components(m, "CBOC")
    a1, a2 = cboc_amplitudes === nothing ? _cboc_int_amplitudes(m.boc1_power) : cboc_amplitudes
    (a1 > 0 && a2 > 0) ||
        error("cboc_amplitudes must be positive integers; got $cboc_amplitudes")
    Int(a1) + Int(a2) <= typemax(Int8) ||
        error("cboc_amplitudes must satisfy a1 + a2 ≤ $(typemax(Int8)) (Int8 table); got $cboc_amplitudes")
    # `boc2_sign`: CBOC(+) for Galileo E1B, CBOC(−) for E1C. CodeLUT.CBOC bakes
    # `a1·BOC(m1,n) + a2·BOC(m2,n)`, so a negative a2 gives the anti-phase composite.
    a2_signed = m.boc2_sign < 0 ? Int8(-Int(a2)) : Int8(a2)
    CodeLUT.CBOC(m.boc1.m, m.boc2.m, m.boc1.n, Int8(a1), a2_signed)
end
_codelut_modulation(m::BOCcos) = CodeLUT.BOCcos(m.m, m.n)

"""
    build_signal_lut(modulation, codes, sec::SecondaryCode; cboc_amplitudes = nothing)
        -> SignalLUT

Bake every PRN's fully-modulated replica into one `SignalLUT` matrix, eagerly at signal
construction: each PRN's padded table, from `CodeLUT.code_replica(primary, mod; secondary,
max_bake = typemax(Int16))`, becomes column `prn`. The shared metadata (P, table_length,
period_subchips) is read from the first PRN.

`codes` is the signal's primary code matrix (`Int8`, unpadded); the ±1 primary for PRN `prn` is
`sign(codes[:, prn])`. `sec` is the signal's [`SecondaryCode`](@ref): short secondaries are
BAKED into each column, a long overlay (GPS L1C-P, 1800 chips) is kept RESIDUAL in
`SignalLUT.secondary` (per-PRN for a `PerPRNSecondaryCode`, else one shared column) for runtime
application.

**Throws an `ArgumentError`** for modulations the LUT can't bake — a BOC(m,n) whose chip spans
no whole number of sub-chips (`2m/n` for sine, `4m/n` for cosine phase), or a TMBOC/CBOC with a
cosine-phased component — so the signal fails to construct.

For `CBOC`, `cboc_amplitudes = nothing` derives the Int8 amplitude pair from the modulation's
`boc1_power` (see `_cboc_int_amplitudes`; the E1 10/11 split yields `(19, 6)`). Pass an explicit
`(a1, a2)` tuple to override.
"""
function build_signal_lut(modulation, codes::AbstractMatrix, sec::SecondaryCode;
                          cboc_amplitudes::Union{Nothing,Tuple{Integer,Integer}} = nothing)
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
    # `mc1.secondary` is `Int8[1]` when baked/none. A residual PER-PRN secondary (L1C-P) keeps one
    # column per PRN; otherwise one shared column.
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
    # Per-sample RMS amplitude (see `SignalLUT.code_amplitude`): PRN-independent (the primary only
    # flips signs), so computed from one PRN's table. Not necessarily PRN 1: undefined PRN slots
    # are all-zero columns (BeiDou B2b_I defines only PRN 6..58) and would give amplitude 0, so
    # take the first column that carries a code. Widen to `Int` before squaring: `abs2` on `Int8`
    # overflows for CBOC's ±25.
    ref_prn = something(
        findfirst(prn -> any(!iszero, @view padded[1:mc1.table.length, prn]), 1:num_prns),
        1,
    )
    tbl = @view padded[1:mc1.table.length, ref_prn]
    code_amplitude = sqrt(sum(x -> abs2(Int(x)), tbl) / length(tbl))
    SignalLUT(
        padded,
        mc1.subchip_factor,
        secmat,
        mc1.table.length,
        mc1.period_subchips,
        code_amplitude,
    )
end

"""
$(SIGNATURES)

Per-sample RMS amplitude of `signal`'s sampled code replica — `sqrt(mean(code^2))`
over one primary code period.

Returns exactly `1` for every ±1 code (BPSK/BOC/TMBOC). For the multi-level CBOC
Int8 approximation (Galileo E1B/E1C) it is `sqrt(a1^2 + a2^2)` for the baked integer
sub-carrier amplitudes `(a1, a2)` — ≈ `19.92` for E1B's default `(19, 6)` — and
reflects any explicit `cboc_amplitudes` override passed to `build_signal_lut`.

Divide a correlation result by this to recover a unit-power amplitude scale; without it a
CBOC prompt is scaled up by the code's integer amplitude relative to a BPSK one.

# Examples
```julia-repl
julia> get_code_amplitude(GPSL1CA())
1.0

julia> get_code_amplitude(GalileoE1B(CBOC))   # sqrt(19^2 + 6^2)
19.924858845171276
```
"""
@inline get_code_amplitude(signal::AbstractGNSSSignal) = signal.lut.code_amplitude

# Primary ±1 column for PRN `prn` from the code matrix — sign of each stored chip.
@inline _primary_col(codes::AbstractMatrix, prn::Integer) =
    Int8[Int8(sign(codes[i, prn])) for i in 1:size(codes, 1)]

# Residual ±1 secondary vector for PRN `prn`: `Int8[1]` when there is no secondary, else the
# per-PRN/-shared ±1 chips.
function _residual_secondary_for_prn(sec::SecondaryCode, prn::Integer)
    Ls = secondary_code_length(sec)
    Ls == 1 ? Int8[1] : Int8[Int8(secondary_value(sec, prn, s)) for s in 0:Ls-1]
end

# Transient, zero-copy `CodeLUT.ModulatedCode` over PRN `prn`'s column of `lut.padded`
# (`chips`/`padded` are unit-stride `SubArray{Int8}` views), with the shared metadata and the
# per-PRN residual secondary view.
@inline function _modulated_code_view(lut::SignalLUT, prn::Int)
    coltot = size(lut.padded, 1)
    # Matching `UnitRange` slices so both views share one `SubArray` type `V`, as the
    # constructor requires.
    padded_col = @view lut.padded[1:coltot, prn]
    chips_col = @view lut.padded[1:lut.table_length, prn]
    table = CodeLUT.CodeTable(chips_col, padded_col, lut.table_length)
    CodeLUT.ModulatedCode(table, lut.subchip_factor, _signal_lut_secondary(lut, prn), lut.period_subchips)
end

"""
$(SIGNATURES)

Generate a sampled code signal for a given PRN.

Allocates a new `Vector{Int8}` of `num_samples` samples. Thin wrapper over the in-place
[`gen_code!`](@ref); see it for the Int8 output semantics (±1, or the CBOC integer
approximation), fractional-phase support, and the
`sampling_frequency ≥ code_frequency · subchip_factor` requirement.

# Arguments
- `sampling_frequency`: Sampling frequency (must be ≥ code frequency · subchip_factor)
- `code_frequency`: Code frequency the replica runs at (default: [`get_chip_rate`](@ref))
- `start_phase`: Initial code phase in chips (default: 0.0)
- `start_index`: Index offset (default: 0)

# Examples
```julia-repl
julia> using Unitful: MHz
julia> sampled_code = gen_code(4000, GPSL1CA(), 1, 4MHz)
julia> length(sampled_code)
4000
```
"""
function gen_code(
    num_samples::Integer,
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_chip_rate(signal),
    start_phase = 0.0,
    start_index::Integer = 0,
)
    code = zeros(Int8, num_samples)
    gen_code!(code, signal, prn, sampling_frequency, code_frequency, start_phase, start_index)
end

"""
    gen_code!(sampled_code::AbstractVector{Int8}, signal, prn, sampling_frequency,
              code_frequency = get_chip_rate(signal),
              start_phase = 0.0, start_index_shift = 0) -> sampled_code

Generate the sampled spreading code for PRN `prn` of `signal` in-place, by resampling PRN
`prn`'s fully-modulated baked replica from the signal's embedded `SignalLUT` with the SIMD
windowed-permute / boundary-fill resampler.

# Output
`sampled_code` must be `Vector{Int8}` (or any `AbstractVector{Int8}`). Non-CBOC signals are
±1; CBOC (Galileo E1B) is the multi-level Int8 *integer approximation* of its sqrt-power
subcarrier amplitudes (sign pattern matches the float spec). Thread-safe: no scratch buffer,
and the shared `signal.lut` is read-only.

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
`sampling_frequency ≥ code_frequency · subchip_factor` (else an error is raised).

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
    code_frequency::Frequency = get_chip_rate(signal),
    start_phase = 0.0,
    start_index_shift::Integer = 0,
)
    mc = _modulated_code_view(signal.lut, Int(prn))
    _gen_code_lut_core!(
        sampled_code, mc, _to_hz(sampling_frequency), _to_hz(code_frequency),
        start_phase, start_index_shift,
    )
end

# One-shot resample core over a `CodeLUT.ModulatedCode` (see `_modulated_code_view`). The
# non-baked secondary is applied inside `CodeLUT.generate_code!(out, mc; …)`.
@inline function _gen_code_lut_core!(
    sampled_code::AbstractVector{Int8},
    mc::CodeLUT.ModulatedCode,
    fs::Float64,
    fc::Float64,
    start_phase,
    start_index_shift::Integer,
)
    phase_sub, rem0 = _lut_resample_setup(mc, fs, fc, start_phase, start_index_shift; what = "gen_code!")
    # Parameterless default_backend() const-folds on non-AVX-512 builds; on AVX-512 it is a small
    # Union{AVX512,AVX2,Portable} (runtime Ref read for the #104 demotion), union-split at the
    # barrier below. The table-aware overload would erase even that inference.
    backend = CodeLUT.default_backend()
    CodeLUT.generate_code!(sampled_code, mc;
        code_frequency = fc, sampling_frequency = fs, phase_sub = phase_sub, rem0 = rem0, backend = backend)
    return sampled_code
end

# ─────────────────────────────────────────────────────────────────────────────
# Value-based continuing fill: an immutable `CodeFillEngine` built once per (signal, prn, rate)
# plus an isbits `CodeFillState` threaded by the caller. The DDA setup (~40 ns) runs once in
# `code_engine`; each `gen_code!(out, eng, st)` fills the next `length(out)` samples and returns
# the advanced state, so consecutive fills concatenate to one big generation. The isbits state
# keeps the steady-state fill allocation-free. For fused, register-resident correlation use
# [`code_engine`](@ref)`(signal, prn, fs, fc, Val(K))` / [`code_lookup`](@ref) instead.
# ─────────────────────────────────────────────────────────────────────────────

"""
    CodeFillEngine

Immutable, loop-invariant LUT code-fill engine built by
[`code_engine`](@ref)`(signal, prn, fs, fc)`. It holds the backend's resampler config
(precomputed step ratio + DDA deltas), the residual (non-baked) secondary (e.g. GPS L5I's
NH10), and the initial phase used to seed [`code_state`](@ref).

Build it once, pair it with a [`code_state`](@ref) seed, then call
[`gen_code!`](@ref)`(out, eng, st)` per integration: no rate setup or DDA re-init, just the
windowed permute loop with a single-stream + scalar tail. The engine is read-only and shareable;
stream position lives in the threaded state, so one engine can drive many independent streams.
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
                code_frequency = get_chip_rate(signal);
                start_phase = 0.0, start_index_shift = 0,
                backend = default_backend()) -> CodeFillEngine

Build an immutable, continuing code-fill engine over PRN `prn` of `signal` at the given rate
(the `Val(K)` overload builds the fused register-resident engine instead). The one-time DDA
setup runs here; afterwards seed a [`code_state`](@ref)`(eng)` and call
[`gen_code!`](@ref)`(out, eng, st) -> st` to fill successive chunks with no per-call setup,
threading the returned state for seamless continuation. Non-baked secondaries (e.g. GPS L5I's
NH10) are supported. For a one-shot fill the [`gen_code!`](@ref)`(out, signal, prn, …)` method is
simpler.

Same oversampling requirement and fractional sub-chip phase support as the one-shot
[`gen_code!`](@ref), including its rate quantization: the engine locks the rate to the DDA's
`2^-30` sub-chip grid and never re-anchors to the requested float rate, so the phase drifts by
at most `N · 2^-31` sub-chips over the `N` samples since the `start_phase` anchor. Closed-loop
tracking is unaffected (the engine is rebuilt on every Doppler update); for open-loop
constant-rate streaming, rebuild it with a fresh `start_phase` about every 1e7 samples.
`backend` defaults to the host's best SIMD backend; pass a weaker one (e.g.
`CodeLUT.AVX2()`/`CodeLUT.Portable()`) to test the non-default paths. Forcing a backend the CPU
does not support is invalid.
"""
function code_engine(
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_chip_rate(signal);
    start_phase = 0.0,
    start_index_shift::Integer = 0,
    backend::CodeLUT.Backend = CodeLUT.default_backend(),
)
    mc = _modulated_code_view(signal.lut, Int(prn))
    fc = _to_hz(code_frequency)
    fs = _to_hz(sampling_frequency)
    phase_sub, rem0 = _lut_resample_setup(mc, fs, fc, start_phase, start_index_shift;
        what = "code_engine", hint = "; use gen_code! with the signal directly")
    # The parameterless `default_backend()` const-folds on non-AVX-512 builds, else is a small
    # `Union{AVX512,AVX2,Portable}` (runtime Ref read for #104) union-split at the
    # `_wrap_code_fill_engine` barrier; the table-aware overload's length guard is dead code for
    # GNSS codes and would box the engine. `backend` is exposed so tests can force a weaker one.
    # Resample the baked sub-chip table at fc·P (fixed-point step: one multiply + round).
    sn, sd = CodeLUT._fixed_point_step((fc * mc.subchip_factor) / fs)
    engine = CodeLUT.make_fill_engine(mc.table, sn, sd; phase = phase_sub, rem0 = rem0, backend = backend)
    # `mc.secondary` views `signal.lut.secondary`; copy so the engine does not alias it.
    return _wrap_code_fill_engine(engine, collect(mc.secondary), mc.period_subchips, mc.subchip_factor, sn, sd, phase_sub, Int(rem0))
end

# Function barrier. `make_fill_engine` returns a small Union over the phase type (Int16 vs
# Int32, from the runtime table length); splitting on the concrete `E` here keeps construction
# type-stable instead of boxing the engine on every construction.
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
fully-modulated replica, continuing seamlessly from `st` (the hot path: no `rationalize`, no
DDA re-init), and **return** the state advanced by exactly `length(sampled_code)`. Consecutive
calls threading the returned state concatenate to one big generation. Any non-baked secondary
(e.g. GPS L5I's NH10) is applied per primary period across call boundaries via the state's
absolute sample offset.

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
    # Same period-walk as the one-shot fill, offset by `n0` (see `CodeLUT._apply_secondary!`).
    CodeLUT._apply_secondary!(out, sec, eng.period_subchips, eng.step_num, eng.step_den,
                              eng.phase_sub, eng.rem0, n0)
end

# ─────────────────────────────────────────────────────────────────────────────
# Value-based code engine: adapter over `CodeLUT.code_engine`. Build one engine per interleave
# factor `K` (`Val(K)`), hold K isbits `code_state`s (W apart), and drive them with
# `code_lookup` / `code_advance` — the allocation-free, register-resident counterpart to
# `gen_code!` and partner of SinCosLUT's `carrier_engine`. (Per-backend rates: see iterate.jl.)
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
                code_frequency = get_chip_rate(signal), Val(K);
                start_phase = 0.0, start_index_shift = 0) -> CodeLUT.CodeEngine

Build a loop-invariant, value-based code engine over PRN `prn` of `signal` for a `K`-way
interleaved fused loop. Pair with `K` states `code_state(eng, k)` (`k = 0..K-1`, `W` samples
apart) and drive each with [`code_lookup`](@ref) / [`code_advance`](@ref); nothing is
materialised or heap-allocated. The code-side counterpart to SinCosLUT's `carrier_engine`. Same
oversampling requirement, fractional sub-chip phase support, and `2^-30` rate quantization
(drift ≤ `N · 2^-31` sub-chips over `N` samples) as [`gen_code!`](@ref). A non-baked secondary
(e.g. GPS L5I's NH10) or `sampling_frequency < code_frequency·subchip_factor` raises an error
(use [`gen_code!`](@ref) / [`code_engine`](@ref) without `Val(K)` instead).
"""
function code_engine(signal::AbstractGNSSSignal, prn::Integer, sampling_frequency::Frequency, code_frequency::Frequency, ::Val{K};
                     start_phase = 0.0, start_index_shift::Integer = 0) where {K}
    mc = _modulated_code_view(signal.lut, Int(prn))
    fc = _to_hz(code_frequency)
    fs = _to_hz(sampling_frequency)
    phase_sub, rem0 = _lut_resample_setup(mc, fs, fc, start_phase, start_index_shift;
        what = "code_engine(…, Val(K))", hint = "; use gen_code! / code_engine without Val(K)")
    any(!=(Int8(1)), mc.secondary) && error(
        "code_engine(…, Val(K)) does not support a non-baked secondary (e.g. GPS L5I NH10); use gen_code! / code_engine without Val(K).",
    )
    # Resample the baked sub-chip table at fc·P. Parameterless default_backend(): see
    # `_gen_code_lut_core!`; the table-aware overload would box the engine.
    CodeLUT.code_engine(mc.table, (fc * mc.subchip_factor) / fs, Val(K);
                        phase = phase_sub, rem0 = rem0, backend = CodeLUT.default_backend())
end

code_engine(signal::AbstractGNSSSignal, prn::Integer, sampling_frequency::Frequency, vk::Val; kwargs...) =
    code_engine(signal, prn, sampling_frequency, get_chip_rate(signal), vk; kwargs...)
