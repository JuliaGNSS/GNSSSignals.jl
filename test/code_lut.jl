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

    @testset "iterator (generate_code) matches array fill" begin
        chips = rand(MersenneTwister(7), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for be in _BACKENDS, cps in (0.2046, 1 / 3, 0.999)
            W = _width(be)
            n = 4W * 80
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, cps; backend = be)
            collected = Int8[]
            for code_vec in CL.generate_code(ct, cps, n; backend = be)
                append!(collected, [code_vec[j] for j in 1:W])
            end
            @test collected == out[1:length(collected)]
        end
    end

    @testset "4-way iterator (generate_code4) matches array fill" begin
        chips = rand(MersenneTwister(11), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for be in _BACKENDS, (cps, ph) in ((0.2046, 0), (0.999, 0), (1 / 3, 9))
            W = _width(be)
            n = 4W * 80
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, cps; phase = ph, backend = be)
            collected = Int8[]
            for quad in CL.generate_code4(ct, cps, n; phase = ph, backend = be)
                for v in quad, j in 1:W
                    push!(collected, v[j])
                end
            end
            @test collected == out[1:length(collected)]
        end
    end

    @testset "continuing generator (make_generator + fill_continue!)" begin
        # Concatenating consecutive fills equals one big generation.
        chips = rand(MersenneTwister(13), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        for be in _BACKENDS, cps in (0.2046, 1 / 3)
            sn = _step_num(cps)
            n = 5000
            ref = _ref(chips, sn, 0, n)
            gen = CL.make_generator(ct, cps; backend = be)
            got = Int8[]
            for m in (777, 1, 1024, 64, n - 777 - 1 - 1024 - 64)   # uneven chunk sizes
                buf = Vector{Int8}(undef, m)
                CL.fill_continue!(buf, gen)
                append!(got, buf)
            end
            @test got == ref
        end
    end

    @testset "high-oversampling run-fill path" begin
        # At high oversampling `make_generator` / `generate_code!` switch from the windowed
        # permute to a broadcast run-fill (`_runfill_*`). Its approximate samples-per-chip DDA
        # matches the exact `floor(step_num·n/2^_B)` fixed-point reference for any single fill
        # (the rounding drift only reaches ≤1 sample after ~10⁶ chips — see the drift test
        # below), so here it must be byte-identical. Sweep `m = samples/chip` across the
        # `Val`-specialised ladder AND the `m > 64` generic kernel, with non-power-of-two
        # rates (`m`/`m+1` mixed runs), several phases/lengths, and chunked continuation.
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
                        gen = CL.make_generator(ct, cps; phase = phase, backend = be)
                        got = Int8[]
                        for chunk in (1, 333, 64, 2048, n - 1 - 333 - 64 - 2048)
                            buf = Vector{Int8}(undef, chunk)
                            CL.fill_continue!(buf, gen)
                            append!(got, buf)
                        end
                        @test got == ref
                    end
                end
            end
        end
    end

    @testset "run-fill long-sequence drift + >MAXFILL segmentation" begin
        chips = rand(MersenneTwister(7), Int8[-1, 1], 1023); ct = CL.CodeTable(chips)
        # Drift bound: over a long single fill (millions of chips) the approximate DDA may
        # differ from the exact floor reference, but only by a couple of samples at boundaries.
        for cps in (1 / 8, 0.0613)            # exact and fractional run-fill rates
            sn = _step_num(cps); n = 3_000_000
            ref = _ref(chips, sn, 0, n)
            out = Vector{Int8}(undef, n)
            CL.generate_code!(out, ct, cps; backend = CL.Portable())
            @test count(!=(0), out .- ref) <= 4      # ≤ a few sub-chip-boundary shifts
        end
        # Buffers longer than `_RUNFILL_MAXFILL` exercise the internal segmentation loop (which
        # caps the per-call accumulator). The segmented output must still track the floor
        # reference to within the same small drift bound.
        n = CL._RUNFILL_MAXFILL + 12345; sn = _step_num(1 / 16)
        @test n > CL._RUNFILL_MAXFILL                       # segmentation actually triggered
        out = Vector{Int8}(undef, n)
        CL.generate_code!(out, ct, sn, CL._STEP_DEN; backend = CL.Portable())
        @test count(!=(0), out .- _ref(chips, sn, 0, n)) <= 4
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
        mods = [CL.LOC(), CL.BOC(1), CL.BOC(2), CL.TMBOC(3, [k in (0, 4, 6) for k in 0:10])]
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
                # iterator path too
                col = Int8[]
                for v in CL.generate_code(ct, cps, n; backend = CL.AVX2()), j in 1:32
                    push!(col, v[j])
                end
                @test col == ref[1:length(col)]
            end
            if Sys.ARCH === :aarch64   # NEON: same Int32-phase long-table path, W = 16
                out = Vector{Int8}(undef, n); CL.generate_code!(out, ct, cps; backend = CL.Neon())
                @test out == ref
                @test CL.default_backend(ct) isa CL.Neon   # not Portable
                col = Int8[]
                for v in CL.generate_code(ct, cps, n; backend = CL.Neon()), j in 1:16
                    push!(col, v[j])
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
# Public adapter: `CodeGeneratorLUT4` — the 4-wide, finite counterpart to
# `CodeGeneratorLUT`. It wraps `CodeLUT.generate_code4`, so its 4-tuple of
# `Vec{W,Int8}` chunks (4·W samples per step) must concatenate byte-exactly to the
# array produced by the plan `gen_code!` at the same rate.
# ─────────────────────────────────────────────────────────────────────────────
@testset "CodeGeneratorLUT4" begin
    # SIMD width of the host's default backend (what CodeGeneratorLUT4 uses internally).
    Wdef = _width(CL.default_backend())

    # Mix of baked-secondary / no-secondary signals AND non-baked-secondary signals:
    # GPSL5I (NH10) and GPSL1C_P (its 1800-chip overlay is too long to bake) both yield a
    # residual secondary that the 4-wide iterator folds in per primary period. All must
    # concatenate byte-exactly to the array `gen_code!` (which applies the same secondary).
    sigs = [GPSL1CA(), GPSL1C_D(), GalileoE1B_BOC11(), GPSL5I(), GPSL1C_P()]
    for signal in sigs
        plan = CodeReplicaLUT(signal, 1)
        P = plan.mc.subchip_factor
        fc = get_code_frequency(plan.signal)
        # A couple of rates, all satisfying fs ≥ fc·P.
        for fs in (Float64(fc) * P, Float64(fc) * P * 2.5)
            nsteps = 6
            num_samples = nsteps * 4 * Wdef
            it = CodeGeneratorLUT4(plan, fs, fc, num_samples)

            @test Base.IteratorSize(typeof(it)) == Base.HasLength()
            @test length(it) == nsteps
            @test eltype(typeof(it)) == eltype(typeof(it.inner))
            @test eltype(typeof(it)) <: NTuple{4,<:CL.SIMD.Vec{Wdef,Int8}}

            # Byte-exact: concatenated 4-tuples == first length(it)·4W samples of gen_code!.
            oracle = Vector{Int8}(undef, num_samples)
            gen_code!(oracle, plan, fs, fc)
            collected = Int8[]
            for (a, b, c, d) in it
                for v in (a, b, c, d), j in 1:Wdef
                    push!(collected, v[j])
                end
            end
            @test length(collected) == nsteps * 4 * Wdef
            @test collected == oracle[1:length(collected)]
        end
    end

    @testset "allocation-free steady iteration" begin
        # Both a baked-secondary signal (L1CA, no per-chunk secondary work) and a non-baked
        # one (L5I NH10, whose sign-fold runs every step) must iterate allocation-free.
        for signal in (GPSL1CA(), GPSL5I())
            plan = CodeReplicaLUT(signal, 1)
            P = plan.mc.subchip_factor
            fc = get_code_frequency(plan.signal)
            fs = Float64(fc) * P
            it = CodeGeneratorLUT4(plan, fs, fc, 6 * 4 * Wdef)
            # Warm up, then assert a single step allocates nothing.
            st = iterate(it)
            @test st !== nothing
            _, state = st
            step!(itr, s) = iterate(itr, s)
            step!(it, state)                       # compile
            # 0-alloc fused iteration is verified on Julia ≥ 1.11; 1.10's iteration-state
            # inference leaks a small box through the wrapped iterator, so only assert where
            # it's reliable (the inner generate_code4 is 0-alloc on all backends).
            if VERSION >= v"1.11"
                @test (@allocated step!(it, state)) == 0
            end
        end
    end

    @testset "errors" begin
        n = 4 * Wdef
        # fs < fc·subchip_factor → error (the one remaining construction error; a non-baked
        # secondary is now supported, see the main testset above).
        plan = CodeReplicaLUT(GPSL1C_D(), 1)   # BOC(1,1) → subchip_factor 2
        P = plan.mc.subchip_factor
        @test P > 1
        fc1c = get_code_frequency(plan.signal)
        @test_throws ErrorException CodeGeneratorLUT4(plan, Float64(fc1c) * P / 2, fc1c, n)
    end

    # End-to-end public adapter at HIGH oversampling, where the engine uses the broadcast
    # run-fill rather than the permute. The one-shot `gen_code!(out, plan, …)` and a *warm*
    # continuing `gen_code!(out, gen)` must produce the identical fill, and chunked
    # continuation must concatenate to one big generation — exercising the run-fill across
    # call boundaries together with the per-primary-period secondary negate (GPS L5I NH10).
    @testset "plan/generator run-fill at high oversampling (L1CA + L5I secondary)" begin
        for (signal, osr) in ((GPSL1CA(), 40), (GPSL5I(), 32))   # L5I carries the NH10 secondary
            plan = CodeReplicaLUT(signal, 1)
            fc = get_code_frequency(plan.signal)
            fs = Float64(fc) * osr * plan.mc.subchip_factor
            N = 7000
            # one-shot (builds + run-fills from phase 0)
            oneshot = Vector{Int8}(undef, N)
            gen_code!(oneshot, plan, fs, fc)
            # warm continuing generator, single fill of N
            gen = CodeGeneratorLUT(plan, fs, fc)
            warm = Vector{Int8}(undef, N)
            gen_code!(warm, gen)
            @test warm == oneshot
            # chunked continuation == one big generation (crosses run + secondary-period boundaries)
            gen2 = CodeGeneratorLUT(plan, fs, fc)
            got = Int8[]
            for chunk in (1, 999, 64, 4096, N - 1 - 999 - 64 - 4096)
                buf = Vector{Int8}(undef, chunk)
                gen_code!(buf, gen2)
                append!(got, buf)
            end
            @test got == oneshot
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# `CodeGeneratorLUT` register iteration (`for v in gen`) with a non-baked secondary.
# Folding the per-primary-period sign into each yielded `Vec{W,Int8}` chunk must reproduce
# the array `gen_code!` byte-for-byte — including chunks that straddle a secondary-period
# boundary (the per-lane sign-vector path) and the high-oversampling run-fill engine.
# ─────────────────────────────────────────────────────────────────────────────
@testset "CodeGeneratorLUT iteration with non-baked secondary" begin
    # (signal, oversampling-over-subchip): L5I NH10 (permute + run-fill), L1C_P 1800-chip
    # overlay (long TMBOC table), and L1CA (no secondary — regression guard for the fast path).
    cases = [
        (GPSL5I(),   3),    # permute path, NH10
        (GPSL5I(),   16),   # run-fill path, NH10
        (GPSL1C_P(), 2),    # long table + 1800-chip overlay
        (GPSL1CA(),  4),    # no secondary
    ]
    for (signal, osr) in cases
        plan = CodeReplicaLUT(signal, 1)
        P = plan.mc.subchip_factor
        fc = get_code_frequency(plan.signal)
        fs = Float64(fc) * P * osr
        per = plan.mc.period_subchips
        samples_per_period = round(Int, per * osr)        # = per · fs/(fc·P)
        W = GNSSSignals.CodeLUT.gen_width(CodeGeneratorLUT(plan, fs, fc).engine)
        # Span >2 secondary periods so chunks straddle boundaries (no-op for L1CA's Ls=1).
        n = ((2 * samples_per_period + 3 * W) ÷ W) * W
        oracle = Vector{Int8}(undef, n)
        gen_code!(oracle, plan, fs, fc)
        gen = CodeGeneratorLUT(plan, fs, fc)
        got = Vector{Int8}(undef, n)
        idx = 1
        for v in gen
            idx > n && break
            for j in 1:W
                idx > n && break
                got[idx] = v[j]
                idx += 1
            end
        end
        @test got == oracle
    end
end
