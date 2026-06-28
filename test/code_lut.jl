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

    @testset "fractional sub-chip start phase (rem0) — all backends == shifted reference" begin
        # A fixed-point fractional sub-chip offset `rem0 ∈ [0, 2^_B)` seeds the DDA's running
        # remainder so it tracks the SHIFTED stream `((step_num·n + rem0) >> _B + phase)`. Every
        # PERMUTE backend must match the shifted reference byte-exactly; the high-oversampling
        # broadcast RUN-FILL gets its documented ≤2-sample boundary tolerance. rem0 = 0 ⇒ the
        # original integer-phase output.
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
                            if CL._use_runfill(sn, _SD, be, n)
                                @test count(!=(0), out .- ref) <= 2   # run-fill boundary tolerance
                            else
                                @test out == ref                      # permute: byte-exact
                            end
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
                        gen = CL.make_generator(ct, sn, _SD; phase = phase, rem0 = rem0, backend = be)
                        got = Int8[]
                        for m in (777, 1, 1024, 64, n - 777 - 1 - 1024 - 64)
                            buf = Vector{Int8}(undef, m); CL.fill_continue!(buf, gen); append!(got, buf)
                        end
                        if CL._use_runfill(sn, _SD, be)
                            @test count(!=(0), got .- ref) <= 4   # accumulated run-fill drift
                        else
                            @test got == ref
                        end
                        # value-based engine (permute backends only — run-fill has no engine)
                        if !CL._use_runfill(sn, _SD, be)
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
# `code_advance` over a `CodeReplicaLUT` plan — the allocation-free, register-resident
# counterpart to filling an array with `gen_code!`. A K-way interleaved loop (built with
# `Val(K)`) holds K states `W` samples apart and advances each by K chunks per step; the
# concatenated `Vec{W,Int8}` lookups must reproduce the `gen_code!` array byte-exactly.
# ─────────────────────────────────────────────────────────────────────────────
@testset "code_engine (value-based)" begin
    # SIMD width of the host's default backend (what the engine uses internally).
    Wdef = _width(CL.default_backend())

    # Baked-secondary signals only. Excluded: GPSL5I (non-baked NH10) and GPSL1C_P
    # (its 1800-chip overlay is too long to bake) — both yield a residual secondary
    # and are covered by the error testset below.
    sigs = [GPSL1CA(), GPSL1C_D(), GalileoE1B_BOC11(), GalileoE1B()]   # E1B: CBOC, multi-level long table
    for signal in sigs
        plan = CodeReplicaLUT(signal, 1)
        P = plan.mc.subchip_factor
        fc = get_code_frequency(plan.signal)
        # A couple of rates, all satisfying fs ≥ fc·P.
        for fs in (Float64(fc) * P, Float64(fc) * P * 2.5)
            eng = code_engine(plan, fs, fc, Val(1))
            W = code_width(eng)
            @test W == Wdef
            nsteps = 6
            num_samples = nsteps * 4 * W

            # Byte-exact oracle.
            oracle = Vector{Int8}(undef, num_samples)
            gen_code!(oracle, plan, fs, fc)

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
            eng4 = code_engine(plan, fs, fc, Val(4))
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
        plan = CodeReplicaLUT(GPSL1CA(), 1)
        P = plan.mc.subchip_factor
        fc = get_code_frequency(plan.signal)
        fs = Float64(fc) * P
        eng = code_engine(plan, fs, fc, Val(1))
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
        fc = get_code_frequency(GPSL1CA())
        # Non-baked secondary (GPS L5I NH10) → error at construction.
        l5i = CodeReplicaLUT(GPSL5I(), 1)
        fc5 = get_code_frequency(GPSL5I())
        @test_throws ErrorException code_engine(l5i, Float64(fc5) * l5i.mc.subchip_factor, fc5, Val(1))
        # fs < fc·subchip_factor → error.
        plan = CodeReplicaLUT(GPSL1C_D(), 1)   # BOC(1,1) → subchip_factor 2
        P = plan.mc.subchip_factor
        @test P > 1
        fc1c = get_code_frequency(plan.signal)
        @test_throws ErrorException code_engine(plan, Float64(fc1c) * P / 2, fc1c, Val(1))
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
# Galileo E1B CBOC support in the LUT path. The subcarrier's two sqrt-power amplitudes are
# baked as an Int8 *integer approximation* (default (19,6) ≈ (√(10/11), √(1/11)), ratio √10);
# the permute/run-fill backends carry the resulting multi-level Int8 values verbatim, so the
# resampled replica reproduces the float `gen_code!` (identical signs; correlation set by the
# amplitude ratio). `cboc_amplitudes` lets the caller trade accuracy for smaller magnitudes.
# ─────────────────────────────────────────────────────────────────────────────
@testset "Galileo E1B CBOC (LUT integer approximation)" begin
    signal = GalileoE1B(); prn = 1
    fc = get_code_frequency(signal)
    P = 12                                            # 2·lcm(1,6)

    @testset "plan construction + baked composite values" begin
        plan = CodeReplicaLUT(signal, prn)           # default (19,6)
        @test plan.mc.subchip_factor == P
        @test length(plan.mc.table) == get_code_length(signal) * P
        @test sort(unique(plan.mc.table.chips)) == Int8[-25, -13, 13, 25]   # ±(19±6)
        # the documented default really is (19,6)
        explicit = CodeReplicaLUT(signal, prn; cboc_amplitudes = (19, 6))
        @test plan.mc.table.chips == explicit.mc.table.chips
    end

    @testset "custom amplitudes scale the table; sign pattern is amplitude-independent" begin
        base = CodeReplicaLUT(signal, prn)
        for (amps, vals) in (((2, 1), Int8[-3, -1, 1, 3]), ((3, 1), Int8[-4, -2, 2, 4]))
            plan = CodeReplicaLUT(signal, prn; cboc_amplitudes = amps)
            @test sort(unique(plan.mc.table.chips)) == vals
            @test sign.(plan.mc.table.chips) == sign.(base.mc.table.chips)
        end
    end

    @testset "amplitude validation" begin
        for bad in ((0, 1), (1, 0), (-1, 1), (1, -2), (100, 100))   # non-positive or a1+a2 > 127
            @test_throws ErrorException CodeReplicaLUT(signal, prn; cboc_amplitudes = bad)
        end
    end

    @testset "matches float gen_code! (signs exact; correlation tracks the amplitude ratio)" begin
        # At fs = fc·P (integer ratio P) every sample is exactly one sub-chip, so the LUT and
        # the float path align sub-chip-for-sub-chip: the sign pattern is identical for ALL
        # amplitudes, and the normalized correlation equals the cosine between (a1,a2) and the
        # float (√10, 1) — ~1 for (19,6), degrading gracefully for coarser ratios.
        fs = fc * P
        N = get_code_length(signal) * P              # one full code period
        flt = Vector{Float32}(undef, N); gen_code!(flt, signal, prn, fs, fc)
        for (amps, mincorr) in (((19, 6), 0.9999), ((3, 1), 0.9998), ((2, 1), 0.987))
            plan = CodeReplicaLUT(signal, prn; cboc_amplitudes = amps)
            lut = Vector{Int8}(undef, N); gen_code!(lut, plan, fs, fc)
            @test sign.(lut) == Int8.(sign.(flt))    # signs identical to the spec
            rho = sum(Float64.(lut) .* flt) /
                  (sqrt(sum(abs2, Float64.(lut))) * sqrt(sum(abs2, flt)))
            @test rho >= mincorr
        end
    end

    @testset "resampling all backends == fixed-point reference (multi-level table)" begin
        plan = CodeReplicaLUT(signal, prn)           # default (19,6) → ±25, ±13
        full = plan.mc.table.chips; Lf = length(full)
        for (a, b) in ((1, 1), (3, 5), (7, 9)), phase in (0, 7)
            sn = _step_num(a / b); psub = phase * P; n = 4 * 64 * 8
            ref = [full[mod(div(sn * (i - 1), _SD) + psub, Lf) + 1] for i in 1:n]
            for be in _BACKENDS                       # E1B table > typemax(Int16): AVX2/NEON use Int32 phase
                out = Vector{Int8}(undef, n)
                CL.generate_code!(out, plan.mc; code_frequency = a, sampling_frequency = b * P,
                                  phase = phase, backend = be)
                @test out == ref
            end
        end
    end

    @testset "high-oversampling run-fill: one-shot == warm == chunked (multi-level)" begin
        plan = CodeReplicaLUT(signal, prn)
        fs = fc * P * 10                              # broadcast run-fill regime
        N = 9000
        oneshot = Vector{Int8}(undef, N); gen_code!(oneshot, plan, fs, fc)
        gen = CodeGeneratorLUT(plan, fs, fc)
        warm = Vector{Int8}(undef, N); gen_code!(warm, gen)
        @test warm == oneshot
        gen2 = CodeGeneratorLUT(plan, fs, fc); got = Int8[]
        for chunk in (1, 999, 64, 4096, N - 1 - 999 - 64 - 4096)
            buf = Vector{Int8}(undef, chunk); gen_code!(buf, gen2); append!(got, buf)
        end
        @test got == oneshot
        @test all(in((-25, -13, 13, 25)), oneshot)   # run-fill preserves the multi-level values
    end

    @testset "fs < fc·subchip_factor errors" begin
        plan = CodeReplicaLUT(signal, prn)
        @test_throws ErrorException gen_code!(Vector{Int8}(undef, 64), plan, fc * (P - 1), fc)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# End-to-end fractional-phase parity against the ORIGINAL float `gen_code!`. With true
# sub-chip phase support the LUT plan reproduces the original's sub-sample phase exactly: at
# fractional `start_phase` / `start_index_shift` (incl. a negative shift) the resampled sign
# pattern is identical to `gen_code!(out, signal, prn, fs, fc, start_phase, start_index_shift)`,
# to within a handful of samples at chip transitions (differing fixed-point rounding only).
# ─────────────────────────────────────────────────────────────────────────────
@testset "fractional start-phase parity vs original gen_code!" begin
    sigs = [GPSL1CA(), GPSL1C_D(), GalileoE1B_BOC11()]   # BPSK, BOC(1,1), BOC(1,1)
    for signal in sigs
        prn = 1
        fc = get_code_frequency(signal)
        plan = CodeReplicaLUT(signal, prn)
        P = plan.mc.subchip_factor
        fs = fc * P * 2.5                         # permute regime (sub-chip oversampled), Unitful Hz
        N = 20000
        for (sp, sis) in ((0.0, 0), (0.25, 0), (0.5, 0), (3.5, 0), (0.0, -1), (1.7, 5))
            orig = Vector{Int8}(undef, N)
            gen_code!(orig, signal, prn, fs, fc, sp, sis)
            lut = Vector{Int8}(undef, N)
            gen_code!(lut, plan, fs, fc, sp, sis)
            # Signs must match (the LUT is ±1; the original may be wider). Allow a handful of
            # chip-boundary samples to differ from fixed-point rounding, as the integer-phase
            # parity comparisons do.
            @test count(sign.(orig) .!= sign.(lut)) <= 8
        end
    end
end
