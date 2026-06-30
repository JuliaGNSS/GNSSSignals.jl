@testset "Common functions for $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    if typeof(signal) <: GalileoE1B
        @test get_code_type(signal) == Float32
    else
        @test get_code_type(signal) == Int16
    end
    @test get_codes(signal) == signal.codes
end

@testset "min_bits_for_code_length" begin
    @test min_bits_for_code_length(GPSL1CA()) == 10  # 1023 requires 10 bits
    @test min_bits_for_code_length(GPSL5I()) == 17  # 10230 * 10 = 102300 requires 17 bits
    @test min_bits_for_code_length(GalileoE1B()) == 12  # 4092 requires 12 bits
end

@testset "get_signal_name" begin
    @test get_signal_name(GPSL1CA()) == "GPS L1 C/A"
    @test get_signal_name(GPSL5I()) == "GPS L5-I"
    @test get_signal_name(GalileoE1B()) == "Galileo E1B"
end

@testset "SecondaryCode dispatch" begin
    # L5-I has a SharedSecondaryCode of length 10 (NH10).
    gpsl5i_sec = get_secondary_code(GPSL5I())
    @test gpsl5i_sec isa SharedSecondaryCode{10}
    @test GNSSSignals.secondary_code_length(gpsl5i_sec) == 10
    # NH10 = (1, 1, 1, 1, -1, -1, 1, -1, 1, -1); prn is ignored for SharedSecondaryCode.
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 0) == 1
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 4) == -1
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 9) == -1
    # Index wraps modulo length.
    @test GNSSSignals.secondary_value(gpsl5i_sec, 1, 10) == 1

    # L1 C/A and E1B have NoSecondaryCode; secondary_value returns `true`
    # (so multiplication is a true no-op preserving eltype).
    for sig in (GPSL1CA(), GalileoE1B())
        sec = get_secondary_code(sig)
        @test sec isa NoSecondaryCode
        @test GNSSSignals.secondary_code_length(sec) == 1
        @test GNSSSignals.secondary_value(sec, 1, 0) === true
        @test GNSSSignals.secondary_value(sec, 7, 42) === true
    end
end

@testset "Base.show for GNSS signals" begin
    io = IOBuffer()
    show(io, GPSL1CA())
    @test occursin("GPSL1CA", String(take!(io)))

    show(io, GPSL5I())
    @test occursin("GPSL5I", String(take!(io)))

    show(io, GalileoE1B())
    @test occursin("GalileoE1B", String(take!(io)))
end

@testset "Broadcasting GNSS signals" begin
    gpsl1ca = GPSL1CA()
    # Test that signals can be broadcast
    phases = [0.0, 1.0, 2.0]
    result = get_code.(gpsl1ca, phases, 1)
    @test length(result) == 3
    @test all(x -> x ∈ [-1, 1], result)
end

@testset "get_modulation type dispatch" begin
    @test get_modulation(GPSL1CA) == GNSSSignals.LOC()
    @test get_modulation(GPSL5I) == GNSSSignals.LOC()
    @test get_modulation(GalileoE1B) isa GNSSSignals.CBOC
end

# ─────────────────────────────────────────────────────────────────────────────
# Fixed-point scalar oracle for the embedded Int8 `gen_code!`, over a signal's baked column.
# This mirrors the kernel exactly: a sample n (1-based) reads the baked sub-chip at index
# `((step_num·(n-1) + rem0) >> _B) + phase_sub` (mod table_length), with any non-baked
# secondary applied per primary period. It IS the definition of correct LUT output (the LUT
# was validated byte-for-sign against the original generator in prior PRs), so we anchor the
# generation tests to it rather than to the now-deleted float generator. Byte-exact in the
# permute regime; the high-oversampling run-fill gets a ≤2-sample boundary tolerance.
# ─────────────────────────────────────────────────────────────────────────────
const _CL = GNSSSignals.CodeLUT

# Shared fixed-point scalar oracle (defined here so test/code_lut.jl and test/signal_lut.jl,
# both included after common.jl, can reuse it). Builds the secondary-FREE resample from the
# baked column `full` exactly as the kernel does — sample n (1-based) reads sub-chip
# `((sn·(n-1) + rem0) >> _B) + psub` (mod L) — then applies any non-baked `sec` as a
# per-primary-period range-negate, mirroring `GNSSSignals._apply_secondary!` precisely
# (period boundaries via `cld((p·per − psub)·sd, sn)`, on the integer `psub` only). Byte-exact
# vs the embedded gen_code! in the permute regime; the run-fill regime gets a ≤2-sample
# boundary tolerance at the call site. `sec` is the per-PRN residual vector (`Int8[1]` = none).
function _gen_code_scalar_oracle(full, sn::Int, sd::Int, rem0::Int, psub::Int, N::Int,
                                 sec::AbstractVector{Int8}, per::Int)
    Lf = length(full)
    ref = Vector{Int8}(undef, N)
    @inbounds for i in 1:N
        subi = ((sn * Int64(i - 1) + Int64(rem0)) >> _CL._B) + psub
        ref[i] = full[mod(subi, Lf) + 1]
    end
    Ls = length(sec)
    if any(!=(Int8(1)), sec)
        p = 0
        @inbounds while true
            T = p * per - psub
            s_p = T <= 0 ? 0 : cld(T * sd, sn)
            s_p >= N && break
            Tn = (p + 1) * per - psub
            s_next = min(Tn <= 0 ? 0 : cld(Tn * sd, sn), N)
            if sec[mod(p, Ls) + 1] == -1
                for n in (s_p + 1):s_next
                    ref[n] = -ref[n]
                end
            end
            p += 1
        end
    end
    ref
end

# Returns (oracle::Vector{Int8}, is_runfill::Bool) for `signal`/`prn` at the given rate/phase.
function _gen_code_oracle(signal, prn, fs_u, fc_u, sp, sis, N)
    lut = signal.lut
    P = lut.subchip_factor; Lf = lut.table_length; per = lut.period_subchips
    full = lut.padded[1:Lf, prn]
    sec = GNSSSignals._signal_lut_secondary(lut, prn)
    fcv = Float64(GNSSSignals._to_hz(fc_u)); fsv = Float64(GNSSSignals._to_hz(fs_u))
    sn, sd = _CL._fixed_point_step((fcv * P) / fsv)
    psub, rem0 = GNSSSignals._subchip_phase_split(sp, sis, fcv, fsv, P, sd)
    ref = _gen_code_scalar_oracle(full, sn, sd, rem0, psub, N, collect(Int8, sec), per)
    runfill = _CL._use_runfill(sn, sd, _CL.default_backend(), N)
    (ref, runfill)
end

# Assert the embedded gen_code! equals the oracle (byte-exact in permute regime; ≤2-sample
# boundary tolerance in the run-fill regime).
function _check_gen_code(signal, prn, fs_u, fc_u, sp, sis, N)
    out = Vector{Int8}(undef, N)
    gen_code!(out, signal, prn, fs_u, fc_u, sp, sis)
    ref, runfill = _gen_code_oracle(signal, prn, fs_u, fc_u, sp, sis, N)
    if runfill
        @test count(!=(0), out .- ref) <= 2
    else
        @test out == ref
    end
    out
end

@testset "gen_code! error paths" begin
    gpsl1ca = GPSL1CA()
    # Sampling frequency below code frequency · subchip_factor errors.
    code = zeros(Int8, 100)
    @test_throws ErrorException gen_code!(
        code, gpsl1ca, 1, 500e3Hz, get_code_frequency(gpsl1ca),
    )
end

# High oversampling drives the broadcast run-fill kernel; exercise it against the oracle.
@testset "High oversampling (run-fill) matches oracle" begin
    signal = GPSL1CA()
    sampling_rate = 200e6Hz   # ≈ 195× oversampling → run-fill regime
    samples = 2000
    _check_gen_code(signal, 1, sampling_rate, get_code_frequency(signal), 0.0, 0, samples)
end

@testset "Code generation $(get_signal_name(signal))" for signal in
                                                          [GalileoE1B(), GPSL1CA(), GPSL5I()]
    sampling_rate = 25e6Hz
    samples = 4000
    fc = get_code_frequency(signal)
    out = _check_gen_code(signal, 1, sampling_rate, fc, 0.0, 0, samples)
    @test eltype(out) === Int8
    # gen_code (non-bang) allocates Int8 and matches the in-place fill.
    @test out == gen_code(samples, signal, 1, sampling_rate, fc, 0.0)
    @test eltype(gen_code(samples, signal, 1, sampling_rate, fc, 0.0)) === Int8
end

@testset "Small code generation $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    sampling_rate = 25e6Hz
    samples = 100
    _check_gen_code(signal, 1, sampling_rate, get_code_frequency(signal), 3.5, 0, samples)
end

# Very short outputs (chip-boundary / tail edge cases) must still match the oracle.
@testset "Very short / odd-length code generation $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    sampling_rate = 25e6Hz
    for samples in (1, 2, 3, 5, 7, 63, 65, 127, 129)
        @test_nowarn _check_gen_code(signal, 1, sampling_rate, get_code_frequency(signal), 3.5, 0, samples)
    end
end

@testset "Code generation for different units" begin
    sampling_rate = 25MHz
    signal = GPSL1CA()
    samples = 1000
    fc = get_code_frequency(signal)
    out = _check_gen_code(signal, 1, sampling_rate, fc, 0.0, 0, samples)
    @test out == gen_code(samples, signal, 1, sampling_rate, fc, 0.0)
end

@testset "Code generation with start_phase bigger than code_length" begin
    signal = GPSL1CA()
    sampling_rate = 2.5e6Hz
    num_samples = 4000
    fc = get_code_frequency(signal)
    out = _check_gen_code(signal, 1, sampling_rate, fc, 2065.0, 0, num_samples)
    @test out == gen_code(num_samples, signal, 1, sampling_rate, fc, 2065.0)
end

@testset "Code generation $(get_signal_name(signal)) with different index" for signal in [
    GalileoE1B(),
    GPSL1CA(),
    GPSL5I(),
]
    sampling_rate = 25e6Hz
    samples = 4002
    fc = get_code_frequency(signal)
    out = _check_gen_code(signal, 1, sampling_rate, fc, 0.0, -1, samples)
    @test out == gen_code(samples, signal, 1, sampling_rate, fc, 0.0, -1)
end

@testset "Code generation with large start_phase (no overflow / BoundsError)" begin
    # Large start_phase values used to overflow the old fixed-point accumulator into negative
    # indices; the LUT DDA wraps mod table length. Still must run cleanly and match the oracle.
    signal = GPSL1CA()
    sampling_freq = 5.0e6Hz
    code_frequency = 1.0230022937236385e6Hz  # Actual value from tracking crash
    start_phase = 14113.513288791713
    start_index_shift = -2
    @test_nowarn _check_gen_code(signal, 1, sampling_freq, code_frequency, start_phase, start_index_shift, 1023)
end

@testset "Code generation with negative start_phase (no overflow / BoundsError)" begin
    signal = GPSL1CA()
    sampling_freq = 5.0e6Hz
    code_frequency = get_code_frequency(signal) + 4000Hz * get_code_center_frequency_ratio(signal)
    start_phase = -1000.0
    start_index_shift = 0
    @test_nowarn _check_gen_code(signal, 1, sampling_freq, code_frequency, start_phase, start_index_shift, 1023)
end

# Secondary-code application across small / tail-only buffers: a buffer at
# `start_phase = k·primary_length` selects secondary chip k; relative to chip 0 the whole
# buffer must be sign-flipped iff the two secondary chips differ. Keep N below one primary
# period so no secondary transition happens inside the buffer. (Equivalent coverage to the old
# tail-only-buffer regression test, against the embedded Int8 gen_code!.)
@testset "secondary applied for small / tail-only buffers" begin
    # GPS L1C-P (PerPRNSecondaryCode) at the chip-aligned 12.276 MHz rate.
    let signal = GPSL1C_P(), prn = 2, primary = get_code_length(signal),
        sr = 12.276e6Hz, cf = 1023e3Hz, sec = GNSSSignals.get_secondary_code(signal)
        s0 = GNSSSignals.secondary_value(sec, prn, 0)
        for sec_offset in 0:5, N in (12, 24, 100)
            buf0 = Vector{Int8}(undef, N); gen_code!(buf0, signal, prn, sr, cf, 0.0, 0)
            buf1 = Vector{Int8}(undef, N); gen_code!(buf1, signal, prn, sr, cf, Float64(sec_offset * primary), 0)
            sk = GNSSSignals.secondary_value(sec, prn, sec_offset)
            @test buf1 == (Int(sk) * Int(s0) > 0 ? buf0 : .-buf0)
        end
    end
    # GPS L5-I (SharedSecondaryCode = NH10).
    let signal = GPSL5I(), prn = 1, primary = get_code_length(signal),
        sr = 25e6Hz, cf = 10230e3Hz, sec = GNSSSignals.get_secondary_code(signal)
        s0 = GNSSSignals.secondary_value(sec, prn, 0)
        for sec_offset in 0:9, N in (12, 24, 100)
            buf0 = Vector{Int8}(undef, N); gen_code!(buf0, signal, prn, sr, cf, 0.0, 0)
            buf1 = Vector{Int8}(undef, N); gen_code!(buf1, signal, prn, sr, cf, Float64(sec_offset * primary), 0)
            sk = GNSSSignals.secondary_value(sec, prn, sec_offset)
            @test buf1 == (Int(sk) * Int(s0) > 0 ? buf0 : .-buf0)
        end
    end
end
