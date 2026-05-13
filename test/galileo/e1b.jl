@testset "Galileo E1B" begin
    galileo_e1b = GalileoE1B()
    @test @inferred(get_band(galileo_e1b)) == L1()
    @test @inferred(get_center_frequency(galileo_e1b)) == 1.57542e9Hz
    @test @inferred(get_code_length(galileo_e1b)) == 4092
    @test @inferred(get_secondary_code_length(galileo_e1b)) == 1
    @test @inferred(get_secondary_code(galileo_e1b)) isa NoSecondaryCode
    @test @inferred(get_code(galileo_e1b, 0, 1)) ≈ 1 * sqrt(10 / 11) + 1 * sqrt(1 / 11)
    @test @inferred(get_code(galileo_e1b, 0.0, 1)) ≈ 1 * sqrt(10 / 11) + 1 * sqrt(1 / 11)
    @test @inferred(get_code(galileo_e1b, 0.5, 1)) ≈ -1 * sqrt(10 / 11) + 1 * sqrt(1 / 11)
    @test @inferred(get_code(galileo_e1b, 1.0, 1)) ≈ 1 * sqrt(10 / 11) + 1 * sqrt(1 / 11)
    @test @inferred(get_code(galileo_e1b, 1.5, 1)) ≈ -1 * sqrt(10 / 11) + 1 * sqrt(1 / 11)
    @test @inferred(GNSSSignals.get_code_unsafe(galileo_e1b, 0.0, 1)) ≈
          1 * sqrt(10 / 11) + 1 * sqrt(1 / 11)
    @test @inferred(get_data_frequency(galileo_e1b)) == 250Hz
    @test @inferred(get_code_frequency(galileo_e1b)) == 1023e3Hz
    @test get_signal_name(galileo_e1b) == "Galileo E1B"

    @test GNSSSignals.get_code_factor(galileo_e1b) == 1

    @test get_code_spectrum(galileo_e1b, 0) == 0.0

    @test get_code_spectrum(galileo_e1b, 2 * get_code_frequency(galileo_e1b)) == 0.0
    @test get_code_spectrum(galileo_e1b, -2 * get_code_frequency(galileo_e1b)) == 0.0

    @test get_code_center_frequency_ratio(galileo_e1b) ≈ 1 / 1540

    # `get_code_type` for CBOC signals comes from
    # `promote_type(eltype(codes), typeof(modulation.boc1_power))`. Pin
    # it to `Float32` so a future widening of `boc1_power` to `Float64`
    # (or eltype change) is caught here rather than silently flipping
    # every `gen_code!` buffer/allocation downstream.
    @test get_code_type(galileo_e1b) === Float32
end

@testset "CBOC `multiply_with_subcarrier!` rejects integer buffers" begin
    # CBOC amplitudes are `sqrt(boc1_power)` and `sqrt(1 - boc1_power)`
    # — irrational (~0.953 and ~0.302 for E1B's 10/11). Integer
    # buffers used to throw a confusing `InexactError` deep in the
    # inner loop; the integer dispatch now throws an `ArgumentError`
    # at the call site naming the offending buffer type and pointing
    # callers at `get_code_type(signal) === Float32`. Float64 still
    # works.
    modulation = get_modulation(GalileoE1B())
    fs = 25e6Hz
    cf = 1023e3Hz

    for buf in (zeros(Int16, 100), zeros(Int32, 100), zeros(Int64, 100))
        @test_throws ArgumentError GNSSSignals.multiply_with_subcarrier!(
            buf, modulation, fs, cf, 0.0, 0,
        )
    end

    # The error message names the eltype and points at `get_code_type`.
    err = try
        GNSSSignals.multiply_with_subcarrier!(zeros(Int16, 100), modulation, fs, cf, 0.0, 0)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("Int16", err.msg)
    @test occursin("get_code_type", err.msg)

    # Float64 path: same algorithm as Float32, should run and produce
    # finite samples. The peak amplitude is
    # `sqrt(boc1_power) + sqrt(1 - boc1_power)` (≈ 1.255 for E1B), so
    # the bound is √2 rather than 1.
    buf64 = ones(Float64, 100)
    GNSSSignals.multiply_with_subcarrier!(buf64, modulation, fs, cf, 0.0, 0)
    @test all(isfinite, buf64)
    @test maximum(abs, buf64) <= sqrt(2) + 1e-6
end

@testset "Galileo E1B BOC(1,1) approximation" begin
    e1b = GalileoE1B_BOC11()
    @test @inferred(get_band(e1b)) == L1()
    @test @inferred(get_center_frequency(e1b)) == 1.57542e9Hz
    @test @inferred(get_code_length(e1b)) == 4092
    @test @inferred(get_secondary_code(e1b)) isa NoSecondaryCode
    @test @inferred(get_secondary_code_length(e1b)) == 1
    @test @inferred(get_data_frequency(e1b)) == 250Hz
    @test @inferred(get_code_frequency(e1b)) == 1023e3Hz
    @test get_signal_name(e1b) == "Galileo E1B (BOC(1,1) approximation)"

    # Modulation differs from full E1B's CBOC.
    @test @inferred(get_modulation(e1b)) == BOCsin(1, 1)
    @test GNSSSignals.get_code_factor(e1b) == 1

    # Output element type stays integer (the whole point of dropping
    # the BOC(6,1) component).
    @test get_code_type(e1b) === Int16

    # Primary code matches the full GalileoE1B's primary code chip-for-chip.
    @test get_codes(e1b) == get_codes(GalileoE1B())
end

@testset "GalileoE1B_BOC11 sample-stream matches PocketSDR reference (PRN 1, 12.276 MHz)" begin
    # Independent cross-check against the BOC(1,1) replica produced by
    # PocketSDR's `gen_code_E1B` at 2 samples per chip, upsampled by 6×
    # to 12 samples per chip. At `sampling_frequency = 12 × code_frequency`
    # the primary-chip and BOC half-cycle boundaries are both
    # sample-aligned, so the comparison is bit-exact.
    #
    # PocketSDR encodes the E1B primary code with the opposite chip-sign
    # convention from Julia (the ICD doesn't pin which bit value maps to
    # which sign); the fixture has been sign-flipped at generation time
    # to match Julia's convention so this test exercises the BOC
    # alignment specifically, not the chip-sign question.
    fixture_dir = joinpath(@__DIR__, "fixtures")
    n_samples = 12 * 4092
    ref = _load_packed_hex_fixture(
        joinpath(fixture_dir, "galileo_e1b_boc11_prn1_fs12chip.hex.gz"),
        n_samples,
    )

    sampling_rate = 12.276e6Hz
    code_rate = 1023e3Hz

    buf = zeros(Int16, n_samples)
    gen_code!(buf, GalileoE1B_BOC11(), 1, sampling_rate, code_rate, 0.0, 0)
    @test buf == ref
end
