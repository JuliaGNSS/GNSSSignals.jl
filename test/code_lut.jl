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

_width(be) = be isa CL.AVX512 ? 64 : be isa CL.AVX2 ? 32 : 1

# Backends to exercise: Portable always; AVX-512 / AVX2 only when the host supports them.
const _BACKENDS = (CL.Portable(),
                   (CL.HOST_FEATURES.avx512vbmi ? (CL.AVX512(),) : ())...,
                   (CL.HOST_FEATURES.avx2       ? (CL.AVX2(),)   : ())...)

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
                    be isa CL.AVX2 && length(mc) > typemax(Int16) && continue
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
            if CL.HOST_FEATURES.avx2
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
        end
        # A normal-length table passes the AVX2 length check (lengths > typemax(Int32) throw,
        # but that is not allocatable here).
        @test CL._check_avx2_length(CL.CodeTable(ones(Int8, 8)))
    end
end
