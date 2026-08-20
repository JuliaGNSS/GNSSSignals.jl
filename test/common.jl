@testset "Common functions for $(get_signal_name(signal))" for signal in [
    GalileoE1B(),
    GalileoE1C(),
    GPSL1CA(),
    GPSL5I(),
]
    if get_modulation(signal) isa GNSSSignals.CBOC
        @test get_code_type(signal) == Float32
        # Multi-level CBOC: amplitude is the RMS of the baked ±13/±25 sub-carrier table,
        # sqrt(19^2 + 6^2) for E1's (19, 6) split (both E1B and E1C, boc2_sign aside).
        tbl = Int.(signal.lut.padded[1:signal.lut.table_length, 1])
        @test get_code_amplitude(signal) ≈ sqrt(sum(abs2, tbl) / length(tbl))
        @test get_code_amplitude(signal) ≈ sqrt(19.0^2 + 6.0^2)
    else
        @test get_code_type(signal) == Int16
        @test get_code_amplitude(signal) == 1.0          # ±1 code
    end
    @test get_codes(signal) == signal.codes
end

@testset "get_code_amplitude is a positive, PRN-independent signal property" begin
    # The amplitude is the RMS of ONE baked column, so it must not be read off a PRN slot
    # the signal leaves undefined: BeiDou B2b_I defines only PRN 6..58 and stores the rest
    # as an all-zero code, which used to make the whole signal report 0.0 — a divide-by-zero
    # for the "divide the correlation by this" contract the accessor documents.
    for S in ALL_SIGNALS
        signal = S()
        amp = get_code_amplitude(signal)
        @test amp > 0.0
        lut = signal.lut
        # Every PRN that carries a code agrees with the signal-level value (only the signs
        # differ between PRNs), and the value is the RMS of that column.
        for prn in axes(lut.padded, 2)
            tbl = Int.(@view lut.padded[1:lut.table_length, prn])
            all(iszero, tbl) && continue
            @test sqrt(sum(abs2, tbl) / length(tbl)) ≈ amp
        end
    end
end

@testset "get_carrier_phase_offset" begin
    # GPS civil signals ride the QUADRATURE carrier, lagging their band's P(Y)
    # in-phase reference by 90° → −π/2 (IS-GPS-200N §3.3.1.5.1, IS-GPS-705 §3.3.1.5).
    @test get_carrier_phase_offset(GPSL1CA()) == -π / 2
    @test get_carrier_phase_offset(GPSL2CM()) == -π / 2   # nominal; CNAV can command in-phase
    @test get_carrier_phase_offset(GPSL2CL()) == -π / 2
    @test get_carrier_phase_offset(GPSL5Q()) == -π / 2
    # Galileo E5aQ pilot LEADS its E5a-I reference by 90° → +π/2 (OS SIS ICD Eq. 1).
    @test get_carrier_phase_offset(GalileoE5aQ()) == π / 2
    # BeiDou uses the same complex-envelope convention as Galileo (`s_data + j·s_pilot`,
    # RF `I·cos − Q·sin`), so its quadrature pilots also LEAD by 90° → +π/2:
    # B2aQ (BDS-SIS-ICD-B2a-1.0 Eq. 4-2/4-3) and the B1C pilot's BOC(1,1) arm, which
    # the BeiDouB1C_P replica tracks (BDS-SIS-ICD-B1C-1.0 Eq. 4-3, Table 4-2).
    @test get_carrier_phase_offset(BeiDouB2aQ()) == π / 2
    @test get_carrier_phase_offset(BeiDouB1C_P()) == π / 2
    # In-phase references and in-phase components default to 0.
    @test get_carrier_phase_offset(GPSL5I()) == 0.0
    @test get_carrier_phase_offset(GalileoE5aI()) == 0.0
    @test get_carrier_phase_offset(BeiDouB2aI()) == 0.0
    @test get_carrier_phase_offset(BeiDouB1C_D()) == 0.0
    # Single-component BeiDou signals are each their band's I component per their ICDs.
    @test get_carrier_phase_offset(BeiDouB1I()) == 0.0
    @test get_carrier_phase_offset(BeiDouB3I()) == 0.0
    @test get_carrier_phase_offset(BeiDouB2bI()) == 0.0
    # GPS L1C rides the same P(Y) in-phase carrier (IS-GPS-800J §3.2.1.6.1) → 0,
    # i.e. 90° off C/A on the same band.
    @test get_carrier_phase_offset(GPSL1C_D()) == 0.0
    @test get_carrier_phase_offset(GPSL1C_P()) == 0.0
    @test get_carrier_phase_offset(GalileoE1B()) == 0.0   # E1B/E1C anti-phase is in the code, not the carrier
    @test get_carrier_phase_offset(GalileoE1C()) == 0.0
    # C/A and L1C are physically in quadrature with each other.
    @test get_carrier_phase_offset(GPSL1C_P()) - get_carrier_phase_offset(GPSL1CA()) == π / 2
    # Type- and instance-forms agree.
    @test get_carrier_phase_offset(GPSL5Q) == get_carrier_phase_offset(GPSL5Q())
    @test get_carrier_phase_offset(GPSL1CA) == get_carrier_phase_offset(GPSL1CA())
end

@testset "min_bits_for_code_length" begin
    @test min_bits_for_code_length(GPSL1CA()) == 10  # 1023 requires 10 bits
    @test min_bits_for_code_length(GPSL5I()) == 17  # 10230 * 10 = 102300 requires 17 bits
    @test min_bits_for_code_length(GalileoE1B()) == 12  # 4092 requires 12 bits
    @test min_bits_for_code_length(GalileoE1C()) == 17  # 4092 * 25 = 102300 requires 17 bits
end

@testset "get_signal_name" begin
    @test get_signal_name(GPSL1CA()) == "GPS L1 C/A"
    @test get_signal_name(GPSL5I()) == "GPS L5-I"
    @test get_signal_name(GalileoE1B()) == "Galileo E1B"
    @test get_signal_name(GalileoE1C()) == "Galileo E1C"

    # Type-level dispatch works without constructing the signal (which would read
    # the code file from disk); the instance path forwards to it.
    for S in ALL_SIGNALS
        @test @inferred(get_signal_name(S)) isa String
    end
    @test get_signal_name(GPSL2CL) == "GPS L2CL"
    @test get_signal_name(GalileoE1B_BOC11) == "Galileo E1B (BOC(1,1) approximation)"
    @test get_signal_name(GPSL1CA()) == get_signal_name(GPSL1CA)
    @test get_signal_name(GalileoE5aQ()) == get_signal_name(GalileoE5aQ)
end

@testset "get_signal_id" begin
    @test @inferred(get_signal_id(GPSL1CA())) === :GPSL1CA
    @test @inferred(get_signal_id(GPSL5I())) === :GPSL5I
    @test @inferred(get_signal_id(GalileoE1B())) === :GalileoE1B
    @test @inferred(get_signal_id(GalileoE1C())) === :GalileoE1C
    @test @inferred(get_signal_id(BeiDouB1I())) === :BeiDouB1I
    @test @inferred(get_signal_id(BeiDouB1C_P())) === :BeiDouB1C_P
    @test @inferred(get_signal_id(BeiDouB2aI)) === :BeiDouB2aI
    # Type-level dispatch works without constructing the signal.
    @test @inferred(get_signal_id(GPSL1CA)) === :GPSL1CA
    # Finer than the band id: same PRN on two bands → two distinct signal ids.
    @test get_signal_id(GPSL1CA()) !== get_signal_id(GPSL5I())
    # BeiDou B1C data and pilot share a carrier but are distinct signals.
    @test get_signal_id(BeiDouB1C_D()) !== get_signal_id(BeiDouB1C_P())
end

@testset "get_constellation_id" begin
    @test @inferred(get_constellation_id(GPSL1CA())) === :GPS
    @test @inferred(get_constellation_id(GPSL5I())) === :GPS
    @test @inferred(get_constellation_id(GalileoE1B())) === :Galileo
    @test @inferred(get_constellation_id(GalileoE5aI())) === :Galileo
    @test @inferred(get_constellation_id(BeiDouB1I())) === :BeiDou
    @test @inferred(get_constellation_id(BeiDouB2aQ())) === :BeiDou
    # Type-level dispatch works without constructing the signal.
    @test @inferred(get_constellation_id(GalileoE1B)) === :Galileo
    # Coarser than the signal id: two bands of one constellation share an id.
    @test get_constellation_id(GPSL1CA()) === get_constellation_id(GPSL5I())
end

@testset "get_constellation_name" begin
    @test @inferred(get_constellation_name(GPSL1CA())) == "GPS"
    @test @inferred(get_constellation_name(GPSL5I())) == "GPS"
    @test @inferred(get_constellation_name(GalileoE1B())) == "Galileo"
    @test @inferred(get_constellation_name(GalileoE5aI())) == "Galileo"
    @test @inferred(get_constellation_name(BeiDouB1I())) == "BeiDou"
    @test @inferred(get_constellation_name(BeiDouB1C_D())) == "BeiDou"
    # Type-level dispatch works without constructing the signal.
    @test @inferred(get_constellation_name(GalileoE1B)) == "Galileo"
    # Coarser than the signal name: two bands of one constellation share a name.
    @test get_constellation_name(GPSL1CA()) == get_constellation_name(GPSL5I())
    # Display counterpart of the id — same granularity, String instead of Symbol.
    @test get_constellation_name(GPSL1CA()) == String(get_constellation_id(GPSL1CA()))
    @test get_constellation_name(GalileoE1B()) == String(get_constellation_id(GalileoE1B()))

    # Every constellation states its name as a literal rather than taking the derived
    # fallback: that is what makes the call allocation-free and constant-foldable, and
    # value equality alone cannot tell the two apart ("GPS" == String(:GPS)). So check
    # the method dispatch actually reaches, not just what it returns.
    derived = which(get_constellation_name, Tuple{Type{AbstractGNSSSignal}})
    for S in (GPSL1CA, GPSL5I, GalileoE1B, GalileoE5aQ, BeiDouB1I, BeiDouB1C_P)
        @test which(get_constellation_name, Tuple{Type{S}}) !== derived
        # Stated, but still exactly String() of the id — the two cannot drift apart.
        @test get_constellation_name(S) == String(get_constellation_id(S))
    end
end

# A signal of a constellation declared outside the package: it states only the
# constellation id, so the name comes from the fallback.
struct TestOnlySignal <: AbstractGNSSSignal{Matrix{Int16}} end
GNSSSignals.get_constellation_id(::Type{TestOnlySignal}) = :TestOnlyGNSS

@testset "get_constellation_name fallback for a new constellation" begin
    @test get_constellation_id(TestOnlySignal) === :TestOnlyGNSS
    @test get_constellation_name(TestOnlySignal) == "TestOnlyGNSS"
    @test get_constellation_name(TestOnlySignal()) == "TestOnlyGNSS"
    @test get_constellation_name(TestOnlySignal) ==
          String(get_constellation_id(TestOnlySignal))
    # It really is the derived fallback doing the work here.
    @test which(get_constellation_name, Tuple{Type{TestOnlySignal}}) ===
          which(get_constellation_name, Tuple{Type{AbstractGNSSSignal}})
end

@testset "every signal answers every identity accessor" begin
    # A signal added later — a new constellation's, say — goes into ALL_SIGNALS
    # and is held to the whole accessor layer here, so an omission surfaces as a
    # test failure rather than a MethodError at someone's call site.
    for S in ALL_SIGNALS
        @test get_signal_id(S) isa Symbol
        @test get_signal_name(S) isa String
        @test get_constellation_id(S) isa Symbol
        @test get_constellation_name(S) isa String
        @test get_band(S) isa Band
        @test get_band_id(S) isa Symbol
        @test get_band_name(S) isa String
        @test get_time_system(S) isa TimeSystem
        @test get_time_system_id(S) isa Symbol
        @test get_time_system_name(S) isa String
        @test get_center_frequency(S) isa Frequency
    end
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

    show(io, GalileoE1C())
    @test occursin("GalileoE1C", String(take!(io)))
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
    @test get_modulation(GalileoE1C) isa GNSSSignals.CBOC
    # E1C is CBOC(−): the BOC(6,1) component is in anti-phase.
    @test get_modulation(GalileoE1C).boc2_sign == -1
    @test get_modulation(GalileoE1B).boc2_sign == 1

    # The instance form forwards to the type form, so the modulation the code LUT
    # is baked from (type form, via build_signal_lut) can never drift from the one
    # get_code_spectrum reports (instance form).
    #
    # Equal values are not enough to pin that: a signal that states its modulation
    # twice passes as long as the two copies still agree. So also require that the
    # one generic forwarder is what dispatch actually reaches for every signal — a
    # per-signal `get_modulation(::MySignal)` would intercept it here.
    generic_forwarder = which(get_modulation, Tuple{AbstractGNSSSignal})
    for S in ALL_SIGNALS
        @test @inferred(get_modulation(S())) === get_modulation(S)
        @test which(get_modulation, Tuple{S}) === generic_forwarder
    end
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
# both included after common.jl, can reuse it). This is the INDEPENDENT definition of correct
# output: for each sample n (1-based) it forms the absolute sub-chip phase
# `A = ((sn·(n-1) + rem0) >> _B) + psub`, reads the baked chip `full[mod(A, Lf) + 1]`, and — for
# a non-baked `sec` — multiplies by the secondary chip of that sample's OWN primary period
# `fld(A, per)`. Deriving the period from the very same `A` as the chip lookup (rather than
# re-deriving period boundaries the way `_apply_secondary!` does) is deliberate: it makes the
# oracle a genuine check on the boundary math instead of a mirror of it. Byte-exact vs the
# embedded gen_code! in the permute regime; the run-fill regime gets a ≤2-sample boundary
# tolerance at the call site. `sec` is the per-PRN residual vector (`Int8[1]` = none).
function _gen_code_scalar_oracle(full, sn::Int, sd::Int, rem0::Int, psub::Int, N::Int,
                                 sec::AbstractVector{Int8}, per::Int)
    Lf = length(full); Ls = length(sec)
    has_sec = any(!=(Int8(1)), sec)
    ref = Vector{Int8}(undef, N)
    @inbounds for i in 1:N
        A = ((sn * Int64(i - 1) + Int64(rem0)) >> _CL._B) + psub
        chip = full[mod(A, Lf) + 1]
        has_sec && (chip *= sec[mod(fld(A, per), Ls) + 1])
        ref[i] = chip
    end
    ref
end

# Returns the oracle::Vector{Int8} for `signal`/`prn` at the given rate/phase.
function _gen_code_oracle(signal, prn, fs_u, fc_u, sp, sis, N)
    lut = signal.lut
    P = lut.subchip_factor; Lf = lut.table_length; per = lut.period_subchips
    full = lut.padded[1:Lf, prn]
    sec = GNSSSignals._signal_lut_secondary(lut, prn)
    fcv = Float64(GNSSSignals._to_hz(fc_u)); fsv = Float64(GNSSSignals._to_hz(fs_u))
    sn, sd = _CL._fixed_point_step((fcv * P) / fsv)
    psub, rem0 = GNSSSignals._subchip_phase_split(sp, sis, fcv, fsv, P, sd)
    _gen_code_scalar_oracle(full, sn, sd, rem0, psub, N, collect(Int8, sec), per)
end

# Assert the embedded gen_code! equals the oracle — byte-exact in BOTH regimes (the exact
# boundary fill removed the old run-fill's ≤2-sample tolerance).
function _check_gen_code(signal, prn, fs_u, fc_u, sp, sis, N)
    out = Vector{Int8}(undef, N)
    gen_code!(out, signal, prn, fs_u, fc_u, sp, sis)
    @test out == _gen_code_oracle(signal, prn, fs_u, fc_u, sp, sis, N)
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
                                                          [GalileoE1B(), GalileoE1C(), GPSL1CA(), GPSL5I()]
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
    GalileoE1C(),
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
    GalileoE1C(),
    GPSL1CA(),
    GPSL5I(),
]
    sampling_rate = 25e6Hz
    for samples in (1, 2, 3, 5, 7, 63, 65, 127, 129)
        @test_nowarn _check_gen_code(signal, 1, sampling_rate, get_code_frequency(signal), 3.5, 0, samples)
    end
end

# Regression: a fractional start phase (nonzero sub-chip residual `rem0`) together with a
# NON-baked secondary (GPS L5I's NH10, GPS L1C-P's overlay) must still flip the secondary sign
# on exactly the sample the DDA places at each primary-period boundary. The earlier code
# derived the boundary from `psub` alone (ignoring `rem0`), misplacing the flip by one sample
# at each transition. This needs enough samples to cross several secondary periods AND a
# fractional phase — the short/integer-phase cases above never exercised it. Byte-exact in the
# permute regime, so `_check_gen_code` asserts equality against the independent oracle.
@testset "Fractional phase + non-baked secondary $(get_signal_name(signal))" for signal in [
    GPSL5I(),
    GPSL1C_P(),
]
    fc = get_code_frequency(signal)
    P = signal.lut.subchip_factor
    per = signal.lut.period_subchips
    sampling_rate = fc * P * 5 // 2                       # permute regime, sub-chip oversampled
    sec = GNSSSignals._signal_lut_secondary(signal.lut, 1); Ls = length(sec)
    # Size the buffer to cross the first secondary sign transition (a couple of periods past it),
    # so the test is non-vacuous for both the 10-chip NH10 (L5I) and the 1800-chip overlay
    # (L1C-P, whose primary period — hence per-period sample count — is ~12× longer).
    s0 = sec[1]
    flip = findfirst(p -> sec[mod(p, Ls) + 1] != s0, 1:Ls)
    @test flip !== nothing                                # residual secondary must actually vary
    samples_per_period = per * 5 ÷ 2                      # per / cps, cps = 2/5
    num_samples = (flip + 2) * samples_per_period
    for (sp, sis) in ((0.3, 0), (0.7, 0), (1.7, 0), (0.0, -1), (2.4, 3))
        one = _check_gen_code(signal, 1, sampling_rate, fc, sp, sis, num_samples)
        # continuing engine (chunked) must reproduce the one-shot across period + call boundaries
        eng = code_engine(signal, 1, sampling_rate, fc; start_phase = sp, start_index_shift = sis)
        st = code_state(eng); got = Int8[]
        h = num_samples ÷ 2
        for ch in (1, h, 64, num_samples - 1 - h - 64)
            b = Vector{Int8}(undef, ch); st = gen_code!(b, eng, st); append!(got, b)
        end
        @test got == one
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
    GalileoE1C(),
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

# The composites `get_relative_power` splits: the real components of each, whose
# values must sum to 1. The BOC(1,1) stand-ins are not listed — they duplicate a
# component rather than adding one — and are covered by their own testset below.
const POWER_COMPOSITES = (
    (GPSL1C_D, GPSL1C_P),         # GPS L1: L1C, a 75/25 split
    (GPSL2CM, GPSL2CL),           # GPS L2: L2C
    (GPSL5I, GPSL5Q),             # GPS L5
    (GalileoE1B, GalileoE1C),     # Galileo E1
    (GalileoE5aI, GalileoE5aQ),   # Galileo E5a
    (BeiDouB1C_D, BeiDouB1C_P),   # BeiDou B1C, the other 75/25 split
    (BeiDouB2aI, BeiDouB2aQ),     # BeiDou B2a
    (BeiDouB1I,),                 # single-component: each its own composite
    (BeiDouB3I,),
    (BeiDouB2bI,),
)

# The constellation-and-band group each signal belongs to. Every signal the package
# defines appears exactly once, which the exhaustiveness check below enforces.
const POWER_GROUPS = (
    (GPSL1CA, GPSL1C_D, GPSL1C_P),
    (GPSL2CM, GPSL2CL),
    (GPSL5I, GPSL5Q),
    (GalileoE1B, GalileoE1C, GalileoE1B_BOC11, GalileoE1C_BOC11),
    (GalileoE5aI, GalileoE5aQ),
    # Constellation-and-band grouping keeps BeiDou apart from GPS/Galileo even where
    # the carrier is shared: B1C rides L1 and B2a rides L5, but a BeiDou value is
    # only comparable with another BeiDou value from the same band.
    (BeiDouB1C_D, BeiDouB1C_P),
    (BeiDouB2aI, BeiDouB2aQ),
    (BeiDouB1I,),
    (BeiDouB3I,),
    (BeiDouB2bI,),
)

@testset "get_relative_power" begin
    # The value is the ICD power split itself, stated as a plain constant rather
    # than derived, so each is exact.
    @test get_relative_power(GalileoE1B()) === 0.5
    @test get_relative_power(GalileoE1C()) === 0.5
    @test get_relative_power(GalileoE5aI()) === 0.5
    @test get_relative_power(GalileoE5aQ()) === 0.5
    @test get_relative_power(GPSL5I()) === 0.5      # ICD tabulates I5/Q5 equal
    @test get_relative_power(GPSL5Q()) === 0.5
    @test get_relative_power(GPSL2CM()) === 0.5     # CM/CL time-multiplexed
    @test get_relative_power(GPSL2CL()) === 0.5
    @test get_relative_power(GPSL1C_P()) === 0.75   # 75/25 split (IS-GPS-800J)
    @test get_relative_power(GPSL1C_D()) === 0.25
    @test get_relative_power(BeiDouB2aI()) === 0.5  # ICD states data:pilot 1:1
    @test get_relative_power(BeiDouB2aQ()) === 0.5
    @test get_relative_power(BeiDouB1C_P()) === 0.75  # the other 75/25 split (ICD 1:3)
    @test get_relative_power(BeiDouB1C_D()) === 0.25
    # Single-component signals: each is the only OS component on its carrier, so
    # each is its own composite.
    @test get_relative_power(BeiDouB1I()) === 1.0
    @test get_relative_power(BeiDouB3I()) === 1.0
    @test get_relative_power(BeiDouB2bI()) === 1.0

    # Being splits, the components of one composite sum to exactly 1, and the two
    # pilot/data ratios are exactly 3 — no tolerance needed for either.
    for composite in POWER_COMPOSITES
        @test sum(get_relative_power, composite) === 1.0
    end
    @test get_relative_power(GPSL1C_P()) / get_relative_power(GPSL1C_D()) === 3.0
    @test get_relative_power(BeiDouB1C_P()) / get_relative_power(BeiDouB1C_D()) === 3.0

    # The B1C value is the ICD component split, not the 29/44 the BOC(1,1) pilot
    # replica captures — same convention as the Galileo _BOC11 stand-ins below.
    @test get_relative_power(BeiDouB1C_P()) != 29 / 44

    @testset "GPS L1: C/A against L1C" begin
        # The one value that is not a bare split, and the reason the scale spans a
        # band rather than a composite: C/A is a separate signal on the L1 carrier,
        # not a third component of L1C, so it is scaled against that composite.
        # −158.5 dBW against L1C's −157.0 is 1.5 dB below it.
        @test get_relative_power(GPSL1CA()) ≈ 10.0^(-1.5 / 10)
        # Hence the pilot sits 0.25 dB above C/A and the data component 4.52 below.
        @test 10log10(get_relative_power(GPSL1C_P()) / get_relative_power(GPSL1CA())) ≈
              0.2506126339170007
        @test 10log10(get_relative_power(GPSL1C_D()) / get_relative_power(GPSL1CA())) ≈
              -4.520599913279623

        # Normalised weights over the three signals a joint L1 loop would combine.
        w = map(get_relative_power, (GPSL1CA(), GPSL1C_D(), GPSL1C_P()))
        @test sum(w ./ sum(w)) ≈ 1.0
        @test collect(w ./ sum(w)) ≈
              [0.4145013213281905, 0.14637466966795237, 0.4391240090038571]
    end

    @testset "BOC(1,1) approximations inherit their parent" begin
        # The approximation changes the correlation shape, not the power split.
        @test get_relative_power(GalileoE1B_BOC11()) === get_relative_power(GalileoE1B())
        @test get_relative_power(GalileoE1C_BOC11()) === get_relative_power(GalileoE1C())
    end

    @testset "scale is per constellation and band" begin
        # A value only means something next to another from the same group. Pinning
        # the grouping keeps that visible: it is why no cross-band ratio is offered,
        # and why GalileoE1B and GPSL5I both being 0.5 says nothing about their
        # relative strength.
        for group in POWER_GROUPS
            # A signal listed in the wrong group fails here rather than silently
            # gaining a peer it is not comparable with.
            @test allequal(map(get_constellation_id, group))
            @test allequal(map(get_band_id, group))
        end
        # Groups are disjoint and cover every signal the package defines.
        @test sort!(collect(Iterators.flatten(POWER_GROUPS)); by = nameof) ==
              sort!(collect(ALL_SIGNALS); by = nameof)
    end

    @testset "accessor layer" begin
        # No generic default, so an omitted signal is a MethodError at a caller
        # rather than a test failure — sweep every type the package defines.
        for S in ALL_SIGNALS
            @test @inferred(get_relative_power(S)) isa Float64
            @test 0 < get_relative_power(S) <= 1
            # Instance and type forms agree bit-for-bit (the instance path is a
            # plain forwarder, so a per-signal method stated on the instance
            # instead of the type would show up here).
            @test get_relative_power(S) === get_relative_power(S())
        end
    end

    @testset "folds to a compile-time constant" begin
        # Like the other type-level accessors, calling this must cost nothing. CI runs
        # the suite with coverage instrumentation (Pkg.test(coverage = true)), which
        # inserts :code_coverage_effect statements into the typed IR and allocates in
        # the measured call itself — so pin the fold as "nothing left but a constant
        # return once coverage statements are ignored", and check allocations only in
        # uninstrumented runs.
        f() = get_relative_power(GPSL1CA)
        stmts = filter(code_typed(f, ())[1][1].code) do stmt
            !(stmt isa Expr && stmt.head === :code_coverage_effect)
        end
        @test only(stmts) === Core.ReturnNode(get_relative_power(GPSL1CA))
        if Base.JLOptions().code_coverage == 0
            @test @allocated(get_relative_power(GPSL1CA)) == 0
        end
    end
end
