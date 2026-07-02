# Tests for the vendored `GNSSSignals.CodeLUT` SIMD code resampler (src/code_lut/*).
#
# All names live inside the internal `CodeLUT` submodule; reach them as `CL.X`. The kernel
# uses a *fixed-point* DDA: the chips-per-sample rate is `step_num / 2^_B` (the public
# `cps::Real` API rounds `cps·2^_B` to `step_num`), so the reference uses `div(step_num·m, 2^_B)`
# — NOT the standalone package's `rationalize`-based num/den reference.

using Random

const CL = GNSSSignals.CodeLUT
const _SD = CL._STEP_DEN          # = 2^_B, the fixed-point DDA denominator

# Fixed-point scalar reference: sample n (1-based) → chip index, then chip.
#   index(n) = mod(div(step_num·(n-1), 2^_B) + phase, L);  out[n] = chips[index+1]
function _ref(chips, step_num, phase, n)
    L = length(chips)
    [chips[mod(div(step_num * (i - 1), _SD) + phase, L) + 1] for i in 1:n]
end
# Fractional-phase scalar reference: the DDA computes the SHIFTED stream
#   index(n) = mod(((step_num·(n-1) + rem0) >> _B) + phase, L)   (n = 1-based; rem0 ∈ [0, 2^_B))
# i.e. _ref with a fixed-point fractional sub-chip residual `rem0` seeded into the remainder.
function _ref_rem0(chips, step_num, rem0, phase, n)
    L = length(chips)
    [chips[mod(((step_num * Int64(i - 1) + Int64(rem0)) >> CL._B) + phase, L) + 1] for i in 1:n]
end
# step_num that the public `cps::Real` API derives for a given rate.
_step_num(cps) = round(Int, cps * _SD)

_width(be) = be isa CL.AVX512 ? 64 : be isa CL.AVX2 ? 32 : be isa CL.Neon ? 16 : 1

# Backends to exercise: Portable always; AVX-512 / AVX2 only when the host supports them;
# NEON only on aarch64 (Apple-Silicon CI). HOST_FEATURES is x86-only (`@static`-guarded),
# so guard the access for ARM/other archs, where only Portable (+ NEON on aarch64) runs.
_hasfeat(f) = isdefined(CL, :HOST_FEATURES) && getfield(CL.HOST_FEATURES, f)
const _BACKENDS = (CL.Portable(),
                   (_hasfeat(:avx512vbmi) ? (CL.AVX512(),) : ())...,
                   (_hasfeat(:avx2)       ? (CL.AVX2(),)   : ())...,
                   (Sys.ARCH === :aarch64 ? (CL.Neon(),)   : ())...)

@testset "CodeLUT" begin
    @testset "CodeTable padding" begin
        ct = CL.CodeTable(Int8[1, -1, 1, 1, -1])
        @test length(ct) == 5
        @test ct.chips == Int8[1, -1, 1, 1, -1]
        @test length(ct.padded) >= 5 + CL.WINDOW_PAD
        # window wraps around the code boundary (padding repeats the head)
        @test ct.padded[6:10] == ct.chips[1:5]
        @test_throws ArgumentError CL.CodeTable(Int8[])
    end

    @testset "resampling correctness (all backends == fixed-point reference)" begin
        rng = MersenneTwister(1)
        for L in (1023, 10230, 257, 64)                  # power-of-two and not, short and long
            chips = rand(rng, Int8[-1, 1], L)
            ct = CL.CodeTable(chips)
            for cps in (0.2046, 0.5, 0.999, 1 / 3, 0.07)
                sn = _step_num(cps)
                for phase in (0, 3, L - 1)
                    n = 4096
                    ref = _ref(chips, sn, phase, n)
                    for be in _BACKENDS
                        out = Vector{Int8}(undef, n)
                        CL.generate_code!(out, ct, cps; phase = phase, backend = be)
                        @test out == ref
                    end
                end
            end
        end
    end

    @testset "non-stride-aligned sizes (tail handling)" begin
        chips = rand(MersenneTwister(3), Int8[-1, 1], 1023)
        ct = CL.CodeTable(chips); cps = 0.2046; sn = _step_num(cps)
        for n in (1, 31, 32, 33, 63, 64, 65, 127, 128, 129, 1000, 4097)
            ref = _ref(chips, sn, 0, n)
            for be in _BACKENDS
                out = Vector{Int8}(undef, n)
                CL.generate_code!(out, ct, cps; backend = be)
                @test out == ref
            end
        end
    end

    @testset "explicit step_numerator/step_denominator (2^_B denominator)" begin
        chips = rand(MersenneTwister(5), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        sn = _step_num(1 / 3); n = 8192
        ref = _ref(chips, sn, 0, n)
        for be in _BACKENDS
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, sn, _SD; backend = be)
            @test out == ref
        end
        # denominator must satisfy 0 < den ≤ 2^_B and num ≤ den (oversampling)
        out = Vector{Int8}(undef, 64)
        @test_throws ArgumentError CL.generate_code!(out, ct, _SD + 1, _SD)
        @test_throws ArgumentError CL.generate_code!(out, ct, _SD + 1, _SD - 1)
    end

    @testset "step_denominator != 2^_B rejected (issue #109)" begin
        # Every kernel hardcodes the 2^_B denominator (`>> _B`, `& (2^_B-1)`), so any other
        # denominator must be rejected loudly rather than silently produce garbage. Here
        # 1/2 chips-per-sample expressed as 1/2 (unsupported) must throw, while the exact
        # equivalent 2^(_B-1)/2^_B (supported) must produce the correct resampled sequence.
        chips = rand(MersenneTwister(17), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        n = 4096; sn_half = _SD >> 1   # numerator for exactly 0.5 chips/sample at den = 2^_B
        for be in _BACKENDS
            out = Vector{Int8}(undef, n)
            # Unsupported denominator (2 ≠ 2^_B): all three entry points must throw.
            @test_throws ArgumentError CL.generate_code!(out, ct, 1, 2; backend = be)
            @test_throws ArgumentError CL.code_engine(ct, 1, 2, Val(1); backend = be)
            @test_throws ArgumentError CL.make_fill_engine(ct, 1, 2; backend = be)
            # Supported denominator (== 2^_B): correct resampled sequence.
            CL.generate_code!(out, ct, sn_half, _SD; backend = be)
            @test out == _ref(chips, sn_half, 0, n)
        end
    end

    @testset "real-frequency interface" begin
        chips = rand(MersenneTwister(42), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        n = 20000
        sn = _step_num(1.023e6 / 5e6)
        ref = _ref(chips, sn, 0, n)
        for be in _BACKENDS
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct; code_frequency = 1.023e6, sampling_frequency = 5e6, backend = be)
            @test out == ref
        end
    end

    @testset "value-based engine (K=1) matches array fill" begin
        chips = rand(MersenneTwister(7), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for be in _BACKENDS, cps in (0.2046, 1 / 3, 0.999)
            W = _width(be)
            n = 4W * 80
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, cps; backend = be)
            eng = CL.code_engine(ct, cps, Val(1); backend = be)
            @test CL.code_width(eng) == W
            collected = Int8[]
            st = CL.code_state(eng)
            for _ in 1:(n ÷ W)
                v = CL.code_lookup(eng, st)
                append!(collected, [v[j] for j in 1:W])
                st = CL.code_advance(eng, st)
            end
            @test collected == out[1:length(collected)]
        end
    end

    @testset "value-based engine (K=4 interleaved) matches array fill" begin
        chips = rand(MersenneTwister(11), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for be in _BACKENDS, (cps, ph) in ((0.2046, 0), (0.999, 0), (1 / 3, 9))
            W = _width(be)
            n = 4W * 80
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, cps; phase = ph, backend = be)
            # 4 streams W samples apart, each advancing 4·W samples per step.
            eng = CL.code_engine(ct, cps, Val(4); phase = ph, backend = be)
            collected = Vector{Int8}(undef, n)
            sts = ntuple(k -> CL.code_state(eng, k - 1), 4)
            for step in 0:(n ÷ (4W) - 1)
                base = step * 4W
                for k in 1:4
                    v = CL.code_lookup(eng, sts[k])
                    for j in 1:W
                        collected[base + (k - 1) * W + j] = v[j]
                    end
                end
                sts = ntuple(k -> CL.code_advance(eng, sts[k]), 4)
            end
            @test collected == out
        end
    end

    @testset "continuing fill engine (make_fill_engine + fill_continue!)" begin
        # Concatenating consecutive fills equals one big generation.
        chips = rand(MersenneTwister(13), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for be in _BACKENDS, cps in (0.2046, 1 / 3)
            sn = _step_num(cps)
            n = 5000
            ref = _ref(chips, sn, 0, n)
            eng = CL.make_fill_engine(ct, cps; backend = be); st = CL.fill_state(eng)
            got = Int8[]
            for m in (777, 1, 1024, 64, n - 777 - 1 - 1024 - 64)   # uneven chunk sizes
                buf = Vector{Int8}(undef, m)
                st = CL.fill_continue!(buf, eng, st)
                append!(got, buf)
            end
            @test got == ref
        end
    end

    @testset "high-oversampling boundary-fill path" begin
        # At high oversampling `make_fill_engine` / `generate_code!` switch from the windowed
        # permute to the boundary fill (`_boundary_*`), whose exact magic-reciprocal chip
        # boundaries are byte-identical to the `floor(step_num·n/2^_B)` fixed-point reference
        # for ANY fill length or continuation chunking. Sweep `m = samples/chip` across the
        # store-width ladder (SW = 16/32/64) AND the `m > 62` long-run EXTRAS variant, with
        # non-power-of-two rates (`m`/`m+1` mixed runs), several phases/lengths, and chunked
        # continuation.
        rng = MersenneTwister(99)
        for L in (1023, 257, 64)
            chips = rand(rng, Int8[-1, 1], L); ct = CL.CodeTable(chips)
            # cps chosen to hit a spread of samples-per-chip incl. exact (power-of-two) and
            # fractional, the ladder rungs (≈16,20,24,32,48,64), and the generic m>64 path.
            for cps in (1 / 8, 1 / 16, 0.0613, 1 / 32, 0.019, 1 / 64, 0.013, 0.004)
                sn = _step_num(cps)
                for phase in (0, 7, L - 1)
                    n = 4 * 1024 + 37                     # not stride-aligned
                    ref = _ref(chips, sn, phase, n)
                    for be in _BACKENDS
                        # one-shot
                        out = Vector{Int8}(undef, n)
                        CL.generate_code!(out, ct, cps; phase = phase, backend = be)
                        @test out == ref
                        # continuing, uneven chunks across run boundaries
                        eng = CL.make_fill_engine(ct, cps; phase = phase, backend = be); st = CL.fill_state(eng)
                        got = Int8[]
                        for chunk in (1, 333, 64, 2048, n - 1 - 333 - 64 - 2048)
                            buf = Vector{Int8}(undef, chunk)
                            st = CL.fill_continue!(buf, eng, st)
                            append!(got, buf)
                        end
                        @test got == ref
                    end
                end
            end
        end
    end

    @testset "boundary kernel forced: exact vs oracle across every regime" begin
        # Drive `_generate_boundary!` DIRECTLY (bypassing the `_use_boundary` gate) so every
        # store-width branch (SW = 16/32/64, EXTRAS), the power-of-two shift path and the
        # magic-reciprocal path, code wraps, wind-down clamping and the < SW scalar tail are
        # all exercised — even at rates the dispatcher would send to the permute kernels.
        # The kernel is exact for ANY rate, so the gate can never change the output.
        rng = MersenneTwister(4242)
        for L in (64, 127, 1023, 10230, 767250)
            chips = rand(rng, Int8[-1, 1], L); ct = CL.CodeTable(chips)
            for cps in (0.95, 0.51, 0.5, 1 / 3, 0.25, 0.126, 0.125, 1 / 12.3, 0.0625,
                        1 / 17, 0.02, 1 / 64, 1 / 64.7, 1 / 65, 1 / 128, 1 / 1023, 1 / 5000),
                n in (1, 15, 16, 17, 63, 64, 65, 200, 4096, 12345),
                phase in (0, L - 1), rem0 in (0, 1, _SD ÷ 3, _SD - 1)

                sn = _step_num(cps)
                out = Vector{Int8}(undef, n)
                CL._generate_boundary!(out, ct, sn, _SD, phase, CL._RemT(rem0))
                @test out == _ref_rem0(chips, sn, rem0, phase, n)
            end
        end
    end

    @testset "fractional sub-chip start phase (rem0) — all backends == shifted reference" begin
        # A fixed-point fractional sub-chip offset `rem0 ∈ [0, 2^_B)` seeds the DDA's running
        # remainder so it tracks the SHIFTED stream `((step_num·n + rem0) >> _B + phase)`. Every
        # backend must match the shifted reference byte-exactly — including the
        # high-oversampling boundary fill, whose exact reciprocal has no tolerance carve-out
        # (the old run-fill needed ≤2 samples of slack here). rem0 = 0 ⇒ the original
        # integer-phase output.
        rng = MersenneTwister(123)
        for L in (1023, 257, 64)
            chips = rand(rng, Int8[-1, 1], L); ct = CL.CodeTable(chips)
            for cps in (0.2046, 0.5, 0.999, 1 / 3, 0.07, 1 / 8, 1 / 16, 1 / 32, 0.004)
                sn = _step_num(cps)
                for phase in (0, 3, L - 1)
                    for rem0 in (0, 1, sn ÷ 2, sn - 1, _SD ÷ 4, _SD ÷ 2, _SD - 1)
                        rem0 >= _SD && continue
                        n = 4096
                        ref = _ref_rem0(chips, sn, rem0, phase, n)
                        for be in _BACKENDS
                            out = Vector{Int8}(undef, n)
                            CL.generate_code!(out, ct, sn, _SD; phase = phase, rem0 = rem0, backend = be)
                            @test out == ref                          # byte-exact in BOTH regimes
                        end
                    end
                end
            end
        end
        # rem0 must satisfy 0 ≤ rem0 < step_denominator.
        cterr = CL.CodeTable(rand(rng, Int8[-1, 1], 1023)); out = Vector{Int8}(undef, 64)
        @test_throws ArgumentError CL.generate_code!(out, cterr, _step_num(0.5), _SD; rem0 = _SD)
        @test_throws ArgumentError CL.generate_code!(out, cterr, _step_num(0.5), _SD; rem0 = -1)
    end

    @testset "fractional phase — continuing generator + value engine" begin
        rng = MersenneTwister(321)
        for L in (1023, 257)
            chips = rand(rng, Int8[-1, 1], L); ct = CL.CodeTable(chips)
            for cps in (0.2046, 1 / 3, 0.999, 1 / 16)
                sn = _step_num(cps)
                for phase in (0, 5), rem0 in (0, sn ÷ 2, _SD ÷ 4, _SD - 1)
                    n = 5000
                    ref = _ref_rem0(chips, sn, rem0, phase, n)
                    for be in _BACKENDS
                        # continuing generator across uneven chunks
                        eng = CL.make_fill_engine(ct, sn, _SD; phase = phase, rem0 = rem0, backend = be); st = CL.fill_state(eng)
                        got = Int8[]
                        for m in (777, 1, 1024, 64, n - 777 - 1 - 1024 - 64)
                            buf = Vector{Int8}(undef, m); st = CL.fill_continue!(buf, eng, st); append!(got, buf)
                        end
                        @test got == ref                      # exact in BOTH regimes
                        # value-based engine (always permute-based, any rate)
                        eng = CL.code_engine(ct, sn, _SD, Val(1); phase = phase, rem0 = rem0, backend = be)
                        W = CL.code_width(eng); col = Int8[]; st = CL.code_state(eng)
                        for _ in 1:(n ÷ W)
                            v = CL.code_lookup(eng, st)
                            for j in 1:W; push!(col, v[j]); end
                            st = CL.code_advance(eng, st)
                        end
                        @test col == ref[1:length(col)]
                    end
                end
            end
        end
    end

    @testset "boundary-fill long fills are exact (no drift, no segmentation cap)" begin
        # The old run-fill's rounded reciprocal drifted ~1–2 samples per ~10⁵ chips and needed
        # a per-call accumulator cap (`_RUNFILL_MAXFILL`). The boundary fill's magic-reciprocal
        # boundaries are exact for any dividend a code period can produce, so a multi-million-
        # sample fill is byte-identical to the floor reference with no segmentation at all.
        chips = rand(MersenneTwister(7), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for cps in (1 / 8, 0.0613, 1 / 16)    # pow2 (shift path) and fractional (magic path)
            sn = _step_num(cps); n = 5_000_000
            ref = _ref(chips, sn, 0, n)
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, sn, _SD; backend = CL.Portable())
            @test out == ref
        end
    end

    @testset "modulation (BOC / TMBOC / secondary)" begin
        # Independent reference: build the FULL primary×secondary×subcarrier ±1 table and
        # index it at the sub-chip fixed-point rate `a/b` (sub-chips/sample, ≤ 1). The
        # resampler is driven with code_frequency = a, sampling_frequency = b·P, so the baked
        # sub-chip table (rate code_frequency·P) is resampled at (a·P)/(b·P) = a/b.
        function modref(primary, modulation, secondary, a, b, phase, n)
            Lp = length(primary); Ls = length(secondary); P = CL.subchip_factor(modulation)
            full = Int8[Int8(primary[c+1]) * Int8(secondary[s+1]) * CL._sc_sign(modulation, k, c, P)
                        for s in 0:Ls-1 for c in 0:Lp-1 for k in 0:P-1]
            Lf = length(full); psub = phase * P; sn = _step_num(a / b)
            [full[mod(div(sn * (i - 1), _SD) + psub, Lf) + 1] for i in 1:n]
        end
        rng = MersenneTwister(2)
        Lp = 31
        mods = [CL.LOC(), CL.BOC(1), CL.BOC(2), CL.TMBOC(3, [k in (0, 4, 6) for k in 0:10]),
                CL.CBOC(1, 6, Int8(19), Int8(6))]   # multi-level table (values ±25, ±13)
        for modulation in mods, Ls in (1, 4), (a, b) in ((1, 1), (3, 5), (7, 9))
            P = CL.subchip_factor(modulation)
            primary = rand(rng, Int8[-1, 1], Lp)
            secondary = rand(rng, Int8[-1, 1], Ls)
            n = 4 * 64 * 6
            for phase in (0, 5), bake in (true, false)
                mc = CL.code_replica(primary, modulation; secondary = secondary,
                                     max_bake = bake ? typemax(Int) : 0)
                ref = modref(primary, modulation, secondary, a, b, phase, n)
                for be in _BACKENDS
                    be isa Union{CL.AVX2,CL.Neon} && length(mc) > typemax(Int16) && continue
                    out = Vector{Int8}(undef, n)
                    CL.generate_code!(out, mc; code_frequency = a, sampling_frequency = b * P,
                                      phase = phase, backend = be)
                    @test out == ref
                end
            end
        end
    end

    @testset "AVX2 long tables (Int32 phase)" begin
        # Tables > typemax(Int16) must work on AVX2 via the widened Int32 phase, byte-exact.
        for L in (40000, 122760)
            chips = rand(MersenneTwister(L), Int8[-1, 1], L)
            ct = CL.CodeTable(chips); n = 4 * 32 * 50
            cps = 7 / 9; sn = _step_num(cps)
            ref = _ref(chips, sn, 0, n)
            port = Vector{Int8}(undef, n); CL.generate_code!(port, ct, cps; backend = CL.Portable())
            @test port == ref
            if _hasfeat(:avx2)
                out = Vector{Int8}(undef, n); CL.generate_code!(out, ct, cps; backend = CL.AVX2())
                @test out == ref
                @test CL.default_backend(ct) isa Union{CL.AVX512, CL.AVX2}   # not Portable
                # value-based engine path too (long table → widened Int32 phase)
                eng = CL.code_engine(ct, cps, Val(1); backend = CL.AVX2())
                W = CL.code_width(eng); col = Int8[]; st = CL.code_state(eng)
                for _ in 1:(n ÷ W)
                    v = CL.code_lookup(eng, st)
                    for j in 1:W; push!(col, v[j]); end
                    st = CL.code_advance(eng, st)
                end
                @test col == ref[1:length(col)]
            end
            if Sys.ARCH === :aarch64   # NEON: same Int32-phase long-table path, W = 16
                out = Vector{Int8}(undef, n); CL.generate_code!(out, ct, cps; backend = CL.Neon())
                @test out == ref
                @test CL.default_backend(ct) isa CL.Neon   # not Portable
                eng = CL.code_engine(ct, cps, Val(1); backend = CL.Neon())
                W = CL.code_width(eng); col = Int8[]; st = CL.code_state(eng)
                for _ in 1:(n ÷ W)
                    v = CL.code_lookup(eng, st)
                    for j in 1:W; push!(col, v[j]); end
                    st = CL.code_advance(eng, st)
                end
                @test col == ref[1:length(col)]
            end
        end
        # A normal-length table passes the AVX2 / NEON length check (lengths > typemax(Int32)
        # throw, but that is not allocatable here).
        @test CL._check_windowed_length(CL.CodeTable(ones(Int8, 8)), CL.AVX2())
        @test CL._check_windowed_length(CL.CodeTable(ones(Int8, 8)), CL.Neon())
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# Value-based code engine: `code_engine` / `code_state` / `code_lookup` /
# `code_advance` over a `(signal, prn)` (reading `signal.lut`) — the allocation-free,
# register-resident counterpart to filling an array with `gen_code!`. A K-way interleaved loop
# (built with `Val(K)`) holds K states `W` samples apart and advances each by K chunks per
# step; the concatenated `Vec{W,Int8}` lookups must reproduce the `gen_code!` array byte-exactly.
# ─────────────────────────────────────────────────────────────────────────────
@testset "code_engine (value-based)" begin
    # SIMD width of the host's default backend (what the engine uses internally).
    Wdef = _width(CL.default_backend())

    # Baked-secondary signals only. Excluded: GPSL5I (non-baked NH10) and GPSL1C_P
    # (its 1800-chip overlay is too long to bake) — both yield a residual secondary
    # and are covered by the error testset below.
    sigs = [GPSL1CA(), GPSL1C_D(), GalileoE1B_BOC11(), GalileoE1B()]   # E1B: CBOC, multi-level long table
    for signal in sigs
        prn = 1
        P = signal.lut.subchip_factor
        fc = get_code_frequency(signal)
        # A couple of rates, all satisfying fs ≥ fc·P.
        for fs in (fc * P, fc * P * 2.5)
            eng = code_engine(signal, prn, fs, fc, Val(1))
            W = code_width(eng)
            @test W == Wdef
            nsteps = 6
            num_samples = nsteps * 4 * W

            # Byte-exact oracle: the embedded-LUT one-shot gen_code!.
            oracle = Vector{Int8}(undef, num_samples)
            gen_code!(oracle, signal, prn, fs, fc)

            # K=1: one stream reproduces the oracle chunk by chunk.
            out = Vector{Int8}(undef, num_samples)
            st = code_state(eng); i = 1
            for _ in 1:(num_samples ÷ W)
                v = code_lookup(eng, st)
                for j in 1:W
                    out[i + j - 1] = v[j]
                end
                i += W
                st = code_advance(eng, st)
            end
            @test out == oracle

            # K=4 interleaved: 4 states W apart, each advancing 4·W samples per step.
            eng4 = code_engine(signal, prn, fs, fc, Val(4))
            @test code_width(eng4) == W
            out4 = Vector{Int8}(undef, num_samples)
            sts = ntuple(k -> code_state(eng4, k - 1), 4)
            for step in 0:(nsteps - 1)
                base = step * 4W
                for k in 1:4
                    v = code_lookup(eng4, sts[k])
                    for j in 1:W
                        out4[base + (k - 1) * W + j] = v[j]
                    end
                end
                sts = ntuple(k -> code_advance(eng4, sts[k]), 4)
            end
            @test out4 == oracle
        end
    end

    @testset "allocation-free steady iteration" begin
        signal = GPSL1CA(); prn = 1
        P = signal.lut.subchip_factor
        fc = get_code_frequency(signal)
        fs = fc * P
        eng = code_engine(signal, prn, fs, fc, Val(1))
        # Drive the value-based loop inside a barrier so the engine type is concrete at
        # the @allocated site — both lookup and the renew-by-value advance must be 0-alloc.
        function drive(eng, n)
            W = code_width(eng); st = code_state(eng); a = Int32(0)
            for _ in 1:(n ÷ W)
                a += Int32(sum(code_lookup(eng, st)))
                st = code_advance(eng, st)
            end
            a
        end
        drive(eng, 6 * 4 * code_width(eng))    # compile
        # 0-alloc value-based stepping is verified on Julia ≥ 1.11; 1.10's weaker inference
        # leaks a small box around the isbits state, so only assert where it's reliable.
        if VERSION >= v"1.11"
            @test (@allocated drive(eng, 6 * 4 * code_width(eng))) == 0
        end
    end

    @testset "errors" begin
        # Non-baked secondary (GPS L5I NH10) → error at engine construction.
        l5i = GPSL5I()
        fc5 = get_code_frequency(l5i)
        @test_throws ErrorException code_engine(l5i, 1, fc5 * l5i.lut.subchip_factor, fc5, Val(1))
        # fs < fc·subchip_factor → error.
        l1cd = GPSL1C_D()   # BOC(1,1) → subchip_factor 2
        P = l1cd.lut.subchip_factor
        @test P > 1
        fc1c = get_code_frequency(l1cd)
        @test_throws ErrorException code_engine(l1cd, 1, fc1c * P / 2, fc1c, Val(1))
    end

    # End-to-end public adapter at HIGH oversampling, where the engine uses the broadcast
    # run-fill rather than the permute. The one-shot `gen_code!(out, signal, prn, …)` and a *warm*
    # continuing `gen_code!(out, eng, st)` must produce the identical fill, and chunked
    # continuation (threading the returned state) must concatenate to one big generation —
    # exercising the run-fill across call boundaries together with the per-primary-period
    # secondary negate (GPS L5I NH10).
    @testset "fill-engine run-fill at high oversampling (L1CA + L5I secondary)" begin
        for (signal, osr) in ((GPSL1CA(), 40), (GPSL5I(), 32))   # L5I carries the NH10 secondary
            prn = 1
            fc = get_code_frequency(signal)
            fs = fc * osr * signal.lut.subchip_factor
            N = 7000
            # one-shot embedded gen_code! (builds + run-fills from phase 0)
            oneshot = Vector{Int8}(undef, N)
            gen_code!(oneshot, signal, prn, fs, fc)
            # warm continuing fill engine, single fill of N
            eng = code_engine(signal, prn, fs, fc)
            warm = Vector{Int8}(undef, N)
            gen_code!(warm, eng, code_state(eng))
            @test warm == oneshot
            # chunked continuation == one big generation (crosses run + secondary-period boundaries)
            st = code_state(eng)
            got = Int8[]
            for chunk in (1, 999, 64, 4096, N - 1 - 999 - 64 - 4096)
                buf = Vector{Int8}(undef, chunk)
                st = gen_code!(buf, eng, st)
                append!(got, buf)
            end
            @test got == oneshot
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Issue #107: engine construction must be ZERO-COPY. `_modulated_code_view` passes the padded
# LUT column as a unit-stride `SubArray{Int8}` view; a concrete `padded::Vector{Int8}` field
# would implicitly `convert`/copy the whole column on every build (~750 KB for GPSL2CL). The
# `padded` field of FillEngine512 / FillEngineBoundary / CodeEngine512 must be parametric
# (`V<:AbstractVector{Int8}`, like PreparedCode) so the view is stored without conversion.
# ─────────────────────────────────────────────────────────────────────────────
@testset "engine construction is zero-copy (issue #107)" begin
    sig = GPSL2CL(); prn = 1                    # 767250-chip code → ~750 KB padded column
    fc = get_code_frequency(sig)
    coltot = size(sig.lut.padded, 1)            # bytes copied if the field is concrete
    # Barriers so the built engine type is concrete at the @allocated site.
    build_fill(s, p, fsx, fcx) = code_engine(s, p, fsx, fcx)
    build_val(s, p, fsx, fcx)  = code_engine(s, p, fsx, fcx, Val(1))
    fs_mod  = fc * 4.9                           # ~4.9 samples/chip → FillEngine512 / CodeEngine512
    fs_high = fc * 12                            # ≥10 samples/chip → FillEngineBoundary
    # A zero-copy build allocates only the small transient view/engine (a few KB); the concrete
    # field copies the whole column. Assert well under a full-column copy.
    bound = 64 * 1024
    @test bound < coltot                         # the copy really is much larger than the bound
    build_fill(sig, prn, fs_mod, fc); build_fill(sig, prn, fs_high, fc); build_val(sig, prn, fs_mod, fc)
    @test (@allocated build_fill(sig, prn, fs_mod,  fc)) < bound   # FillEngine512
    @test (@allocated build_fill(sig, prn, fs_high, fc)) < bound   # FillEngineBoundary
    @test (@allocated build_val(sig, prn, fs_mod,   fc)) < bound   # CodeEngine512
end

# ─────────────────────────────────────────────────────────────────────────────
# Galileo E1B CBOC support in the LUT path. The subcarrier's two sqrt-power amplitudes are
# baked as an Int8 *integer approximation* (default (19,6) ≈ (√(10/11), √(1/11)), ratio √10);
# the permute/run-fill backends carry the resulting multi-level Int8 values verbatim, so the
# resampled replica reproduces the float CBOC spec (identical signs; correlation set by the
# amplitude ratio). The embedded `SignalLUT` bakes the default (19,6); custom amplitudes are
# exercised one level down via `CL.code_replica` / `_codelut_modulation`.
# ─────────────────────────────────────────────────────────────────────────────
@testset "Galileo E1B CBOC (LUT integer approximation)" begin
    signal = GalileoE1B(); prn = 1
    fc = get_code_frequency(signal)
    P = 12                                            # 2·lcm(1,6)
    Lp = get_code_length(signal)
    # Primary ±1 chips for this PRN (sign of the stored Int16 chips), as the embedded bake uses.
    primary = Int8[Int8(sign(get_codes(signal)[i, prn])) for i in 1:Lp]
    # Build a baked CBOC ModulatedCode with given integer amplitudes (the embedded path bakes
    # the default (19,6); here we reach one level down for the amplitude sweep).
    cboc_mc(a1, a2) = CL.code_replica(primary, CL.CBOC(1, 6, Int8(a1), Int8(a2));
                                      max_bake = typemax(Int16))

    @testset "embedded bake + baked composite values" begin
        @test signal.lut.subchip_factor == P
        @test signal.lut.table_length == Lp * P
        embedded = signal.lut.padded[1:signal.lut.table_length, prn]   # the baked sub-chip table
        @test sort(unique(embedded)) == Int8[-25, -13, 13, 25]         # ±(19±6)
        # the embedded default really is (19,6)
        @test embedded == cboc_mc(19, 6).table.chips
    end

    @testset "custom amplitudes scale the table; sign pattern is amplitude-independent" begin
        base = cboc_mc(19, 6).table.chips
        for (amps, vals) in (((2, 1), Int8[-3, -1, 1, 3]), ((3, 1), Int8[-4, -2, 2, 4]))
            chips = cboc_mc(amps...).table.chips
            @test sort(unique(chips)) == vals
            @test sign.(chips) == sign.(base)
        end
    end

    @testset "amplitude validation" begin
        cbocmod = get_modulation(signal)
        for bad in ((0, 1), (1, 0), (-1, 1), (1, -2), (100, 100))   # non-positive or a1+a2 > 127
            @test_throws ErrorException GNSSSignals._codelut_modulation(cbocmod, bad)
        end
    end

    @testset "matches float CBOC (signs exact; correlation tracks the amplitude ratio)" begin
        # At fs = fc·P (integer ratio P) every sample is exactly one sub-chip, so the LUT and
        # the float CBOC align sub-chip-for-sub-chip: the sign pattern is identical for ALL
        # amplitudes, and the normalized correlation equals the cosine between (a1,a2) and the
        # float (√10, 1) — ~1 for (19,6), degrading gracefully for coarser ratios. The float
        # reference is the package's own `get_code` CBOC value at each sub-chip phase.
        fs = fc * P
        N = Lp * P                                   # one full code period
        fcv = GNSSSignals._to_hz(fc); fsv = GNSSSignals._to_hz(fs)
        phase = (0:N-1) .* (fcv / fsv)
        flt = Float32.(get_code.(signal, phase, prn))
        for (amps, mincorr) in (((19, 6), 0.9999), ((3, 1), 0.9998), ((2, 1), 0.987))
            mc = cboc_mc(amps...)
            lut = Vector{Int8}(undef, N)
            CL.generate_code!(lut, mc; code_frequency = fcv, sampling_frequency = fsv,
                              backend = CL.default_backend())
            @test sign.(lut) == Int8.(sign.(flt))    # signs identical to the spec
            rho = sum(Float64.(lut) .* flt) /
                  (sqrt(sum(abs2, Float64.(lut))) * sqrt(sum(abs2, flt)))
            @test rho >= mincorr
        end
    end

    @testset "resampling all backends == fixed-point reference (multi-level table)" begin
        full = cboc_mc(19, 6).table.chips; Lf = length(full)   # default (19,6) → ±25, ±13
        mc = cboc_mc(19, 6)
        for (a, b) in ((1, 1), (3, 5), (7, 9)), phase in (0, 7)
            sn = _step_num(a / b); psub = phase * P; n = 4 * 64 * 8
            ref = [full[mod(div(sn * (i - 1), _SD) + psub, Lf) + 1] for i in 1:n]
            for be in _BACKENDS                       # E1B table > typemax(Int16): AVX2/NEON use Int32 phase
                out = Vector{Int8}(undef, n)
                CL.generate_code!(out, mc; code_frequency = a, sampling_frequency = b * P,
                                  phase = phase, backend = be)
                @test out == ref
            end
        end
    end

    @testset "high-oversampling run-fill: one-shot == warm == chunked (multi-level)" begin
        fs = fc * P * 10                              # broadcast run-fill regime
        N = 9000
        oneshot = Vector{Int8}(undef, N); gen_code!(oneshot, signal, prn, fs, fc)
        eng = code_engine(signal, prn, fs, fc)
        warm = Vector{Int8}(undef, N); gen_code!(warm, eng, code_state(eng))
        @test warm == oneshot
        st = code_state(eng); got = Int8[]
        for chunk in (1, 999, 64, 4096, N - 1 - 999 - 64 - 4096)
            buf = Vector{Int8}(undef, chunk); st = gen_code!(buf, eng, st); append!(got, buf)
        end
        @test got == oneshot
        @test all(in((-25, -13, 13, 25)), oneshot)   # run-fill preserves the multi-level values
    end

    @testset "fs < fc·subchip_factor errors" begin
        @test_throws ErrorException gen_code!(Vector{Int8}(undef, 64), signal, prn, fc * (P - 1), fc)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Characterization test for the unified fixed-point arithmetic (issue #110): the windowed DDA
# carry-advance (`_advance_phase` / `_advance_window_512` and the scalar-tail loops that now
# call them) and the secondary-code period-walk (`_apply_secondary!` one-shot vs
# `_apply_secondary_continue!`). Real secondary-code signals are filled at MODERATE
# oversampling (fs = fc·P·osr with 1 < osr < 2) so every backend takes the windowed permute
# regime — NOT the high-oversampling boundary fast path, whose DDA is separate exact
# arithmetic. Reference = the one-shot full-length `gen_code!`; a chunked continuation whose
# chunk sizes straddle both SIMD-window and secondary-period boundaries must reproduce it
# byte-for-byte, cross-validating the one-shot and continuing copies of each unified routine.
# ─────────────────────────────────────────────────────────────────────────────
@testset "unified DDA + secondary period-walk (issue #110 characterization)" begin
    sigs = [GPSL5I(), GPSL5Q(), GalileoE5aI(), GalileoE1C()]   # all carry a non-baked secondary
    for signal in sigs
        prn = 1
        P = signal.lut.subchip_factor
        fc = get_code_frequency(signal)
        @test length(signal.lut.secondary) > 1                # residual (non-baked) secondary
        for osr in (1.5, 1.75)                                 # windowed permute regime on every backend
            fs = Float64(fc) * P * osr
            N = 40000                                          # spans several secondary periods
            # one-shot: windowed DDA + one-shot `_apply_secondary!` (the reference).
            oneshot = Vector{Int8}(undef, N)
            gen_code!(oneshot, signal, prn, fs, fc)
            # a single warm fill (absolute offset 0) must match the one-shot exactly.
            eng = code_engine(signal, prn, fs, fc)
            warm = Vector{Int8}(undef, N)
            gen_code!(warm, eng, code_state(eng))
            @test warm == oneshot
            # chunked continuation (offset > 0): a scalar-tail DDA advance runs at every chunk
            # edge and `_apply_secondary_continue!` places the sign flips across call
            # boundaries — the concatenation must equal the one big generation.
            st = code_state(eng); got = Int8[]
            for chunk in Iterators.cycle((1, 64, 999, 7, 16384, 63, 4096))
                remaining = N - length(got)
                remaining == 0 && break
                buf = Vector{Int8}(undef, min(chunk, remaining))
                st = gen_code!(buf, eng, st)
                append!(got, buf)
            end
            @test got == oneshot
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# CBOC Int8 amplitudes are DERIVED from the modulation's own `boc1_power` (issue #112),
# not hardcoded to E1's (19, 6). The default 10/11 split must still bake to (19, 6); a
# non-10/11 split must bake amplitudes whose ratio tracks √(boc1_power) : √(1-boc1_power),
# NOT the fixed (19, 6). `cboc_amplitudes` stays an explicit override.
# ─────────────────────────────────────────────────────────────────────────────
@testset "CBOC amplitudes derived from boc1_power (#112)" begin
    # `_codelut_modulation(m::CBOC)` (no override) derives the Int8 pair from `m.boc1_power`.
    amps(m) = let c = GNSSSignals._codelut_modulation(m)
        (Int(c.a1), abs(Int(c.a2)))
    end
    # (a) the E1 default (10/11) still reproduces the historical (19, 6)
    @test amps(GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(6, 1), 10 / 11)) == (19, 6)
    # (b) a 50/50 split must NOT be (19, 6): its ratio → √0.5 / √0.5 = 1
    a1, a2 = amps(GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(6, 1), 0.5))
    @test (a1, a2) != (19, 6)
    @test a1 / a2 ≈ sqrt(0.5) / sqrt(1 - 0.5) rtol = 0.05
    @test a1 + a2 <= typemax(Int8)                      # within the Int8 table budget
    # explicit `cboc_amplitudes` still overrides the derived pair
    ovr = GNSSSignals._codelut_modulation(
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(6, 1), 0.5), (19, 6))
    @test (Int(ovr.a1), Int(ovr.a2)) == (19, 6)
end

# ─────────────────────────────────────────────────────────────────────────────
# End-to-end fractional-phase parity against the fixed-point scalar oracle (`_ref_rem0`),
# which IS the definition of correct embedded-LUT output (the LUT was validated byte-for-sign
# against the original generator in prior PRs). At fractional `start_phase` /
# `start_index_shift` (incl. a negative shift) the resampled output must equal the oracle's
# resampled chip from the baked column, byte-exactly in the permute regime.
# ─────────────────────────────────────────────────────────────────────────────
@testset "fractional start-phase parity vs scalar oracle" begin
    sigs = [GPSL1CA(), GPSL1C_D(), GalileoE1B_BOC11()]   # BPSK, BOC(1,1), BOC(1,1)
    for signal in sigs
        prn = 1
        fc = get_code_frequency(signal)
        P = signal.lut.subchip_factor
        fs = fc * P * 2.5                         # permute regime (sub-chip oversampled), Unitful Hz
        N = 20000
        # The baked sub-chip table for this PRN (the column the embedded gen_code! resamples).
        Lf = signal.lut.table_length
        full = signal.lut.padded[1:Lf, prn]
        fcv = Float64(fc); fsv = Float64(fs)
        sn, sd = GNSSSignals.CodeLUT._fixed_point_step((fcv * P) / fsv)
        for (sp, sis) in ((0.0, 0), (0.25, 0), (0.5, 0), (3.5, 0), (0.0, -1), (1.7, 5))
            lut = Vector{Int8}(undef, N)
            gen_code!(lut, signal, prn, fs, fc, sp, sis)
            # Fixed-point scalar oracle at the same fractional phase split the kernel uses.
            psub, rem0 = GNSSSignals._subchip_phase_split(sp, sis, fcv, fsv, P, sd)
            ref = _ref_rem0(full, sn, rem0, psub, N)
            @test lut == ref
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Deterministic, host-independent coverage of the public adapter and the remaining
# adapter/kernel edge paths. Per-backend KERNEL line coverage is already provided by the many
# `for be in _BACKENDS` loops above (CL.generate_code! / code_engine on every host-supported
# backend) plus the Intel-SDE AVX-512 + macOS-NEON CI jobs — coverage tracks executed lines,
# so no public `backend` kwarg on the one-shot signal API is needed to reach them; the
# continuing `code_engine` does expose `backend`, so we sweep it across `_BACKENDS` here.
# ─────────────────────────────────────────────────────────────────────────────
@testset "adapter all-backend sweep + edge coverage" begin
    @testset "code_engine == one-shot, every supported backend (view-backed table)" begin
        # Forcing each supported backend on the continuing engine must reproduce the default
        # one-shot byte-exactly (the permute output is backend-independent), so the AVX2 /
        # Portable adapter paths over the zero-copy column view are covered regardless of the
        # CI host CPU. Forcing a backend the CPU lacks would SIGILL, so we only sweep _BACKENDS.
        signal = GPSL1CA(); prn = 1
        fc = get_code_frequency(signal)
        fs = fc * 2.5                                     # permute regime (sub-chip oversampled)
        oneshot = Vector{Int8}(undef, 20000); gen_code!(oneshot, signal, prn, fs, fc)
        # Function barrier for the @allocated site so the engine/state types are concrete and
        # the isbits state isn't boxed into the testset's soft scope (mirrors `drive` above).
        fill_once(eng, out, st) = gen_code!(out, eng, st)
        for be in _BACKENDS
            eng = code_engine(signal, prn, fs, fc; backend = be)
            warm = Vector{Int8}(undef, 20000); gen_code!(warm, eng, code_state(eng))
            @test warm == oneshot
            # The threaded fill state is isbits, so the steady-state fill is allocation-free.
            # 0-alloc is reliable on Julia ≥ 1.11 (1.10's weaker inference can leak a small box).
            stw = code_state(eng); fill_once(eng, warm, stw)   # compile
            VERSION >= v"1.11" && @test (@allocated fill_once(eng, warm, stw)) == 0
            # A short fill (no full 4W stride) drives the single-window leftover tail of the
            # AVX2/NEON windowed kernel (up to 3 W-blocks before the scalar remainder).
            short = Vector{Int8}(undef, 100)
            sts = code_state(eng); gen_code!(short, eng, sts)
            @test short == oneshot[1:100]
        end
    end

    @testset "non-baked secondary negate across NH10 periods (GPS L5I)" begin
        # The continuing fill engine applies GPS L5I's non-baked NH10 secondary as a per-period
        # sign flip across call boundaries; cross several periods (incl. the -1 chips at NH10
        # indices 4,5,7,9) and require byte-exact agreement with the one-shot gen_code!.
        sig = GPSL5I(); prn = 1
        fc = get_code_frequency(sig); fs = fc * 2
        N = 230_000                       # > 10 NH10 periods (period = 10230·osr samples)
        oneshot = Vector{Int8}(undef, N); gen_code!(oneshot, sig, prn, fs, fc)
        @test any(==(Int8(-1)), GNSSSignals._signal_lut_secondary(sig.lut, prn))   # negate branch can fire
        eng = code_engine(sig, prn, fs, fc); st = code_state(eng)
        got = Int8[]
        for ch in (1, 99_999, 64, 100_000, N - 1 - 99_999 - 64 - 100_000)
            b = Vector{Int8}(undef, ch); st = gen_code!(b, eng, st); append!(got, b)
        end
        @test got == oneshot
    end

    @testset "_subchip_phase_split rounding edge (rem0 carries into θ_int)" begin
        sd = CL._STEP_DEN
        # A fractional start phase whose sub-chip residual rounds up to a whole sub-chip must
        # carry into the integer offset, leaving rem0 in [0, step_den) — the `rem0 >= step_den`
        # branch.
        θ_int, rem0 = GNSSSignals._subchip_phase_split(5 - 0.5 / sd, 0, 1.0, 1.0, 1, sd)
        @test 0 <= rem0 < sd
        @test θ_int == 5                                  # residual carried (4 -> 5)
        θ2, r2 = GNSSSignals._subchip_phase_split(3.25, 0, 1.0, 1.0, 1, sd)  # non-edge
        @test θ2 == 3 && 0 <= r2 < sd
    end

    @testset "very high oversampling (m > 62): EXTRAS strided interior stores" begin
        # Runs longer than the widest splat store (SW = 64) take the `EXTRAS` variant, which
        # adds SW-strided stores inside each run. Exact, incl. size-1 continuation chunks that
        # land on/inside run interiors.
        chips = rand(MersenneTwister(21), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for cps in (1 / 100, 1 / 64.7, 1 / 128, 1 / 5000)
            sn = _step_num(cps)
            n = 300_000
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, sn, _SD; backend = CL.Portable())
            @test out == _ref(chips, sn, 0, n)
        end
        sn = _step_num(1 / 100)
        eng = CL.make_fill_engine(ct, sn, _SD; backend = CL.Portable()); st = CL.fill_state(eng)
        m = 600; ref = _ref(chips, sn, 0, m); got = Int8[]
        for _ in 1:m
            b = Vector{Int8}(undef, 1); st = CL.fill_continue!(b, eng, st); append!(got, b)
        end
        @test got == ref
    end

    # ─────────────────────────────────────────────────────────────────────────────
    # #104: CPU features are detected at PRECOMPILE time and baked into `HOST_FEATURES`. If the
    # pkgimage is precompiled on an AVX-512-VBMI host and later loaded on a CPU without VBMI
    # (JULIA_CPU_TARGET=generic, shared/NFS depot, Docker/CI), the stale const would make
    # `default_backend()` return `AVX512()` and the first `gen_code!` would emit `vpermb` and
    # SIGILL. The fix re-detects features in `__init__` (into a runtime Ref) and selects the
    # backend from *those*. A true SIGILL can't be triggered safely in-process, so we test the
    # pure, CPU-independent selection function + the runtime feature Ref instead.
    @static if Sys.ARCH in (:x86_64, :i686)
        @testset "CPU feature runtime re-validation (#104)" begin
            # Pure backend selection is a function of a plain feature set: a non-VBMI set must
            # NEVER yield AVX512() (its vpermb would SIGILL), a VBMI set may (fast path kept).
            @test CL._select_backend((avx2 = false, avx512vbmi = false)) isa CL.Portable
            @test CL._select_backend((avx2 = true, avx512vbmi = false)) isa CL.AVX2
            @test !(CL._select_backend((avx2 = true, avx512vbmi = false)) isa CL.AVX512)
            @test CL._select_backend((avx2 = true, avx512vbmi = true)) isa CL.AVX512

            # __init__ populated the runtime feature Ref with the CPU actually executing.
            @test CL.RUNTIME_FEATURES[] == CL._x86_features()

            # The AVX-512 path is compiled/selectable ONLY when this build was precompiled on a
            # VBMI host (the `@static` gate on HOST_FEATURES — otherwise its <64 x i8> `vpermb`
            # can't be legalised and even compiling it aborts LLVM). That is also the only build
            # that can over-claim, so the runtime demotion (#104) lives there and is exercised on
            # the AVX-512 SDE CI job. On a non-VBMI build the guarantee is stronger and static:
            # AVX-512 is const-folded away, so gen_code! can never emit `vpermb`.
            if CL.HOST_FEATURES.avx512vbmi
                @test CL.default_backend() === CL._select_backend(CL.RUNTIME_FEATURES[])
                # Simulate a VBMI-baked pkgimage loaded on a non-VBMI CPU by pointing the runtime
                # Ref at a non-VBMI feature set: selection must demote rather than emit `vpermb`.
                saved = CL.RUNTIME_FEATURES[]
                try
                    CL.RUNTIME_FEATURES[] = (avx2 = true, avx512vbmi = false)
                    @test !(CL.default_backend() isa CL.AVX512)
                    @test CL.default_backend() isa CL.AVX2
                    CL.RUNTIME_FEATURES[] = (avx2 = false, avx512vbmi = false)
                    @test CL.default_backend() isa CL.Portable
                finally
                    CL.RUNTIME_FEATURES[] = saved
                end
            else
                @test !(CL.default_backend() isa CL.AVX512)
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Issue #105: the public entry points (`gen_code!` / `gen_code` / `code_engine`) must
# require `Unitful.Frequency` for `sampling_frequency` and `code_frequency`. A bare number
# used to be silently interpreted as Hz — e.g. `code_frequency = 1.023` (meaning MHz) became
# 1.023 Hz, passing the fs≥fc·P guard trivially and returning a constant all-ones buffer
# (0 chip transitions instead of 511). A bare frequency must raise `MethodError` again.
# ─────────────────────────────────────────────────────────────────────────────
@testset "issue #105: public API requires Unitful.Frequency" begin
    signal = GPSL1CA(); prn = 1

    # A bare Real code_frequency must NOT be silently accepted as Hz.
    @test_throws MethodError gen_code!(zeros(Int8, 4000), signal, prn, 4MHz, 1.023)
    @test_throws MethodError gen_code(4000, signal, prn, 4MHz, 1.023)
    @test_throws MethodError code_engine(signal, prn, 4MHz, 1.023)
    @test_throws MethodError code_engine(signal, prn, 4MHz, 1.023, Val(1))
    # A bare Real sampling_frequency is likewise rejected.
    @test_throws MethodError gen_code!(zeros(Int8, 4000), signal, prn, 4e6, 1.023MHz)
    @test_throws MethodError code_engine(signal, prn, 4e6, 1.023MHz)

    # The correct Unitful call resolves the intended 1.023 MHz rate: 511 chip transitions
    # over one 1 ms C/A period (vs 0 for the mis-scaled bare-number bug).
    buf = zeros(Int8, 4000)
    gen_code!(buf, signal, prn, 4MHz, 1.023MHz)
    @test sum(diff(buf) .!= 0) == 511
    @test gen_code(4000, signal, prn, 4MHz, 1.023MHz) == buf

    # The default-frequency path (code_frequency omitted, derived from the signal) still works.
    @test sum(diff(gen_code(4000, signal, prn, 4MHz)) .!= 0) == 511
end
