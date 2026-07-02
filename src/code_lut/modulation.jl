# Subcarrier (BOC/TMBOC) + secondary-code support, by *baking the modulation into the
# resampled table*.
#
# A BOC subcarrier and a secondary (tiered) code are themselves periodic ±1 patterns, so
# the fully-modulated replica is still a ±1 sequence — just at a finer rate. We expand
# the primary code into `P` sub-chips per chip (P = subcarrier factor), with each
# sub-chip carrying the subcarrier sign, and resample that expanded ±1 table in a single
# windowed-permute pass — no separate `multiply_with_subcarrier!`/secondary loop.
#
#   LOC (BPSK):   P = 1,  sub-chip sign = +1
#   BOC(m,1) sin: P = 2m, sub-chip k sign = iseven(k) ? +1 : -1
#   BOC(m,1) cos: P = 4m, sub-chip k sign = iseven((k+m)÷2) ? +1 : -1 (sine shifted a ¼ cycle;
#                 the shift straddles the sine grid, so quarter-sub-chips align every boundary)
#   TMBOC(m2):    P = 2·m2, BOC(1,1) at most chip positions, BOC(m2,1) at `pattern`
#                 positions (the high-rate component sets the resolution)
#   CBOC(m1,m2):  P = 2·lcm(m1,m2), sub-chip value = a1·BOC(m1,1) ± a2·BOC(m2,1) — a multi-
#                 level *integer approximation* (e.g. ±(a1±a2)) of the irrational
#                 sqrt-power CBOC amplitudes. The expanded table holds those Int8 values
#                 verbatim; the permute/run-fill backends are value-agnostic (they gather
#                 chips by position, not by value), so the ±1 machinery resamples them
#                 unchanged. Galileo E1B is CBOC(1,6); the default (a1,a2)=(19,6) ≈
#                 (sqrt(10/11), sqrt(1/11)) (ratio 3.167 ≈ √10) and the caller may pick another.
#
# Resample at `chip_frequency · P` (the sub-chip rate); this needs `fs ≥ chip_frequency·P`
# (sub-chip oversampling ≥ 1), the same window-span condition as the plain code.
#
# Secondary codes multiply whole primary periods. A short secondary is *baked* into the
# table (tile ×Ls); a long one (e.g. L1C-P's 1800-chip overlay) is applied as a per-
# primary-period sign flip — constant over a whole period, so just a range negate.

# ---- modulation descriptors ----
abstract type Modulation end
struct LOC <: Modulation end                    # BPSK, no subcarrier
struct BOC <: Modulation                        # sine-phased BOC(m, 1)
    m::Int
end
BOC() = BOC(1)
struct BOCcos <: Modulation                      # cosine-phased BOC(m, 1)
    m::Int
end
BOCcos() = BOCcos(1)
struct TMBOC <: Modulation                      # BOC(1,1) + BOC(m2,1) at `pattern` positions
    m2::Int
    pattern::Vector{Bool}                       # pattern[(pos mod end)+1] == true → use BOC(m2,1)
end
# Composite BOC: a1·BOC(m1,1) + a2·BOC(m2,1) with *integer* amplitudes (an Int8 approximation
# of the irrational sqrt-power CBOC amplitudes). Galileo E1B is CBOC(m1=1, m2=6); the default
# (a1,a2)=(19,6) approximates (sqrt(10/11), sqrt(1/11)) (ratio 3.167 ≈ √10).
struct CBOC <: Modulation
    m1::Int
    m2::Int
    a1::Int8                                    # amplitude of the BOC(m1,1) component
    a2::Int8                                    # amplitude of the BOC(m2,1) component
end

subchip_factor(::LOC)   = 1
subchip_factor(b::BOC)  = 2 * b.m
subchip_factor(b::BOCcos) = 4 * b.m            # quarter-sub-chips: the ¼-cycle cosine shift halves the grid
subchip_factor(t::TMBOC) = 2 * t.m2
subchip_factor(c::CBOC) = 2 * lcm(c.m1, c.m2)

# Sub-carrier value for sub-chip `k` (0-based) of a chip at primary position `pos`, given P.
# ±1 for BPSK/BOC/TMBOC; a multi-level integer for CBOC (see below).
@inline _sc_sign(::LOC, k, pos, P)  = Int8(1)
@inline _sc_sign(b::BOC, k, pos, P) = iseven(k) ? Int8(1) : Int8(-1)
# Cosine BOC(m,1) at P = 4m: the sub-carrier is the sine BOC shifted by a quarter cycle
# (`get_subcarrier_code(BOCcos, φ) == get_subcarrier_code(BOCsin, φ+¼)`). At the sub-chip
# midpoint φ = pos + (k+½)/(4m), `floor((φ+¼)·2m)` reduces to `pos·2m + (k+m)÷2`; `pos·2m` is
# even, so the sign is position-independent: `iseven((k+m)÷2)`. Matches the float spec exactly.
@inline _sc_sign(b::BOCcos, k, pos, P) = iseven((k + b.m) ÷ 2) ? Int8(1) : Int8(-1)
@inline function _sc_sign(t::TMBOC, k, pos, P)
    if @inbounds t.pattern[mod(pos, length(t.pattern)) + 1]   # BOC(m2,1): flips every sub-chip
        iseven(k) ? Int8(1) : Int8(-1)
    else                                                       # BOC(1,1): +1 first half, -1 second
        k < P ÷ 2 ? Int8(1) : Int8(-1)
    end
end
# CBOC returns the composite *value* a1·b1 + a2·b2 (e.g. ±(a1±a2)), not a bare sign — the Int8
# table carries it verbatim. `div(k·2m, P)` is the BOC(m,1) sub-carrier half-period index at the
# composite resolution P, mirroring the float `get_subcarrier_code`'s `floor(phase·2m)` so the
# sign at every sub-chip matches the spec (`a1 > a2 > 0` ⇒ same sign as the sqrt-power version).
@inline function _sc_sign(c::CBOC, k, pos, P)
    b1 = iseven(div(k * 2 * c.m1, P)) ? Int8(1) : Int8(-1)
    b2 = iseven(div(k * 2 * c.m2, P)) ? Int8(1) : Int8(-1)
    c.a1 * b1 + c.a2 * b2
end

"""
    ModulatedCode

A primary code with its subcarrier (and, if short enough, its secondary code) baked into
a register-resident ±1 `CodeTable`. Build with [`code_replica`](@ref); resample with
`generate_code!(out, mc; code_frequency, sampling_frequency, …)` (a single windowed-
permute pass). `code_frequency` is the **primary chip rate**; the table is resampled at
`code_frequency · subchip_factor`.
"""
struct ModulatedCode{T<:CodeTable,S<:AbstractVector{Int8}}
    table::T
    subchip_factor::Int        # P: output is resampled at chip_frequency · P
    secondary::S               # residual secondary to apply per primary period ([1] if baked/none)
    period_subchips::Int       # sub-chips in one primary period (Lp · P)
end
ModulatedCode(table::T, P, secondary::S, period_subchips) where {T<:CodeTable,S<:AbstractVector{Int8}} =
    ModulatedCode{T,S}(table, P, secondary, period_subchips)

Base.length(mc::ModulatedCode) = length(mc.table)

"""
    code_replica(primary_chips, modulation; secondary=Int8[1], max_bake=typemax(Int16))

Build a [`ModulatedCode`](@ref) from a primary ±1 `primary_chips` vector and a
`modulation` (`LOC()`, `BOC(m)`, `BOCcos(m)`, `TMBOC(m2, pattern)`, or `CBOC(m1, m2, a1, a2)`).
The table is ±1 for all but `CBOC`; for `CBOC` it holds the multi-level integer composite
`a1·BOC(m1,1) ± a2·BOC(m2,1)` (still Int8, so all backends resample it unchanged).
`secondary` is a ±1 overlay
that multiplies whole primary periods; it is baked into the table when
`length(secondary)·length(primary)·subchip_factor ≤ max_bake`, otherwise applied per
period (a cheap range-negate) at generation time.

`max_bake` defaults to `typemax(Int16)` so a baked table never exceeds the AVX2 backend's
addressable length — e.g. GPS L5I's NH10 secondary is *not* baked (10230·10 > 32767), so
the table stays 10230 and AVX2 runs at full speed instead of falling back to scalar. Pass
`max_bake = typemax(Int)` to force baking (single-pass, AVX-512/Portable only) when you
don't need AVX2.
"""
function code_replica(primary::AbstractVector{<:Integer}, modulation::Modulation;
                      secondary::AbstractVector{<:Integer} = Int8[1], max_bake::Integer = typemax(Int16))
    Lp = length(primary); P = subchip_factor(modulation); Ls = length(secondary)
    bake_sec = Ls * Lp * P <= max_bake
    nperiods = bake_sec ? Ls : 1
    expanded = Vector{Int8}(undef, nperiods * Lp * P)
    i = 1
    @inbounds for s in 0:nperiods-1
        sec = bake_sec ? Int8(secondary[s + 1]) : Int8(1)
        for c in 0:Lp-1
            pc = Int8(primary[c + 1]) * sec
            for k in 0:P-1
                expanded[i] = pc * _sc_sign(modulation, k, c, P); i += 1
            end
        end
    end
    residual = bake_sec ? Int8[1] : Int8.(secondary)
    ModulatedCode(CodeTable(expanded), P, residual, Lp * P)
end

"""
    generate_code!(out, mc::ModulatedCode; code_frequency, sampling_frequency, phase=0, backend=…)

Fill `out` with the fully-modulated replica (primary × subcarrier × secondary). `out` is
resampled from the baked table at `code_frequency · subchip_factor`; requires
`sampling_frequency ≥ code_frequency · subchip_factor`. `phase` is an integer primary-chip
offset.
"""
function generate_code!(out::AbstractVector{<:Integer}, mc::ModulatedCode;
                        code_frequency::Real, sampling_frequency::Real,
                        phase::Integer = 0, phase_sub::Integer = Int(phase) * mc.subchip_factor,
                        rem0::Integer = 0, backend::Backend = default_backend(mc.table))
    # `phase_sub` is the integer sub-chip start offset (θ_int); `rem0` the fixed-point
    # fractional sub-chip residual. Default `phase_sub = phase·P, rem0 = 0` reproduces the
    # original integer primary-chip phase. The secondary only needs the integer sub-chip
    # offset (the fractional part never changes which primary period a sample belongs to).
    generate_code!(out, mc.table;
        code_frequency = code_frequency * mc.subchip_factor, sampling_frequency = sampling_frequency,
        phase = phase_sub, rem0 = rem0, backend = backend)
    any(!=(Int8(1)), mc.secondary) &&
        _apply_secondary!(out, mc, code_frequency, sampling_frequency, Int(phase_sub), _RemT(rem0))
    out
end

# Multiply each whole primary period by its secondary chip. The secondary is constant over
# a period (period_subchips sub-chips), so we negate contiguous sample ranges — a handful
# of vectorisable negates, not a per-sample multiply. `phase_sub` is the integer sub-chip
# start offset (θ_int); `rem0` the fractional sub-chip residual seeded into the DDA (both must
# match the baked-table lookup, or the sign flip lands on the wrong sample at a period edge).
function _apply_secondary!(out, mc::ModulatedCode, code_frequency, sampling_frequency, phase_sub,
                           rem0::_RemT = _RemT(0))
    sn, sd = _fixed_point_step(code_frequency * mc.subchip_factor / sampling_frequency)
    _apply_secondary!(out, mc.secondary, mc.period_subchips, sn, sd, Int(phase_sub), Int(rem0), 0)
end

# Period-walk core: multiply each whole primary period covered by `out` by its secondary chip,
# where `out`'s first sample is absolute sample `n0`. `per` = sub-chips per primary period,
# `sn/sd` the fixed-point sub-chip step, `phase_sub` the integer sub-chip start offset (θ_int),
# `r0` the fractional sub-chip residual. `n0 = 0` is the one-shot fill; a positive `n0`
# continues an ongoing fill (see `_apply_secondary_continue!`).
function _apply_secondary!(out, secondary::AbstractVector, per::Int, sn::Int, sd::Int,
                           phase_sub::Int, r0::Int, n0::Int)
    Ls = length(secondary); N = length(out)
    # absolute sample n (0-based) maps to sub-chip floor((n·sn + rem0)/sd) + phase_sub; period p
    # spans sub-chips [p·per, (p+1)·per). First sample of period p: smallest n with that sub-chip
    # ≥ p·per, i.e. n ≥ cld(T·sd − rem0, sn) (T·sd > rem0 whenever T ≥ 1, since rem0 < sd).
    # Start at the period the first emitted sample (absolute n0) falls in — NOT period 0: a
    # negative start_index_shift puts n0 in a negative period, and start_phase ≥ one code period
    # in a period ≥ 1; starting at 0 would apply the wrong secondary chip to the leading samples.
    #
    # OVERFLOW SAFETY (issues #102 / #108). Two products used to overflow here:
    #   • the seed `n0·sn` overflows Int64 once the absolute sample count n0 exceeds
    #     typemax(Int64) ÷ sn (≈ 2^34 for GPS L5I at 20.46 MHz — ~14 min of streaming); the
    #     seed p then wraps hugely negative and the period-walk spins ~10^6+ iterations, so
    #     gen_code! effectively never returns (#102);
    #   • the per-period products `T·sd` / `Tn·sd` (with sd = _STEP_DEN = 2^30) overflow native
    #     `Int` on 32-bit Julia for any real T ≥ 2, so cld(…) turns negative and the negate loop
    #     writes out-of-bounds (#108).
    # Fix: (a) form the seed in Int128; (b) re-anchor the effective sub-chip phase modulo the
    # secondary-code period Q = Ls·per sub-chips, so the running period index p — and hence
    # every per-period product — stays bounded no matter how long streaming continues; (c) do
    # all per-period arithmetic in explicit Int64 (never native Int), like the core DDA kernels.
    # Re-anchoring shifts the effective phase by whole multiples of Q = Ls·per sub-chips, i.e.
    # p by whole multiples of Ls, so `mod(p, Ls)` (the selected secondary chip) is unchanged and
    # the output is byte-identical to the exact walk for every fill length and continuation n0.
    snw = Int64(sn); sdw = Int64(sd); perw = Int64(per)
    # Fold n0 (and its fractional residual r0) into an effective integer sub-chip phase at local
    # sample 0: sub(n0+j) = floor((j·sn + R0)/sd) + phase_sub with R0 = n0·sn + r0 (Int128 seed,
    # the only widening needed — R0 ≥ 0, so its ÷/% are non-negative).
    R0 = Int128(n0) * snw + Int64(r0)
    rd = Int64(R0 % sdw)                       # fractional residual at local 0, 0 ≤ rd < sd
    qd = Int64(R0 ÷ sdw)                       # whole sub-chips absorbed (≤ n0, fits Int64)
    Q  = Int64(Ls) * perw                      # secondary-code period, in sub-chips
    phase0 = Int64(phase_sub) + qd             # effective integer sub-chip phase at local 0
    phase_eff = phase0 - fld(phase0, Q) * Q    # ∈ [0, Q): strips whole multiples of Ls from p
    # Walk periods in the local frame (first emitted sample = local sample 0, residual rd, phase
    # phase_eff). p starts in [0, Ls) and only advances over the periods this fill spans, so
    # p·per and T·sd stay ≤ ~Q·sd ≪ typemax(Int64).
    sub0 = rd ÷ sdw + phase_eff                # = phase_eff (rd < sd)
    p = fld(sub0, perw)
    @inbounds while true
        T = p * perw - phase_eff
        n_start = T <= 0 ? Int64(0) : cld(T * sdw - rd, snw)     # local first sample of period p
        n_start >= N && break
        Tn = (p + 1) * perw - phase_eff
        n_end = min(Tn <= 0 ? Int64(0) : cld(Tn * sdw - rd, snw), Int64(N))  # local end (exclusive)
        if secondary[mod(p, Ls) + 1] == -1
            lo = n_start + Int64(1)                    # 1-based into out (n_start ≥ 0 here)
            hi = n_end                                  # 1-based inclusive
            @simd for n in lo:hi
                out[n] = -out[n]
            end
        end
        p += 1
    end
    out
end
