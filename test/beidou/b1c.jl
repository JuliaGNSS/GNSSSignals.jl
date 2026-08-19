@testset "BeiDou B1C data (D)" begin
    b1c_d = BeiDouB1C_D()
    @test @inferred(get_band(b1c_d)) == L1()
    @test @inferred(get_center_frequency(b1c_d)) == 1_575_420_000Hz
    @test @inferred(get_code_length(b1c_d)) == 10230
    @test @inferred(get_secondary_code_length(b1c_d)) == 1
    @test @inferred(get_secondary_code(b1c_d)) isa NoSecondaryCode
    @test @inferred(get_data_frequency(b1c_d)) == 100Hz
    @test @inferred(get_code_frequency(b1c_d)) == 1_023_000Hz
    @test get_signal_name(b1c_d) == "BeiDou B1C data"
    # Data component is sine-phased BOC(1,1).
    @test @inferred(get_modulation(b1c_d)) == GNSSSignals.BOCsin(1, 1)
    @test GNSSSignals.get_code_factor(b1c_d) == 1
end

@testset "BeiDou B1C pilot (P)" begin
    b1c_p = BeiDouB1C_P()
    @test @inferred(get_band(b1c_p)) == L1()
    @test @inferred(get_center_frequency(b1c_p)) == 1_575_420_000Hz
    @test @inferred(get_code_length(b1c_p)) == 10230
    @test @inferred(get_secondary_code_length(b1c_p)) == 1800
    @test @inferred(get_secondary_code(b1c_p)) isa PerPRNSecondaryCode
    @test @inferred(get_data_frequency(b1c_p)) == 0Hz   # dataless pilot
    @test @inferred(get_code_frequency(b1c_p)) == 1_023_000Hz
    @test get_signal_name(b1c_p) == "BeiDou B1C pilot"
    # ICD signal is QMBOC(6,1,4/33); its real in-phase replica is BOC(1,1) (see docstring).
    @test @inferred(get_modulation(b1c_p)) == GNSSSignals.BOCsin(1, 1)
end

@testset "BeiDou B1C primary codes match the ICD" begin
    codes_d = get_codes(BeiDouB1C_D())
    codes_p = get_codes(BeiDouB1C_P())
    @test size(codes_d) == (10230, 63)
    @test size(codes_p) == (10230, 63)
    _bd_verify_octal(codes_d, ICD_B1C_DATA)
    _bd_verify_octal(codes_p, ICD_B1C_PILOT)
end

@testset "BeiDou B1C pilot secondary codes match the ICD" begin
    b1c_p = BeiDouB1C_P()
    overlay = b1c_p.overlay_codes
    @test size(overlay) == (1800, 63)
    _bd_verify_octal(overlay, ICD_B1C_PILOT_SEC)
    # And via the public secondary interface (period-1800 per-PRN overlay).
    sec = get_secondary_code(b1c_p)
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:1799] == overlay[:, 1]
    @test GNSSSignals.secondary_value(sec, 1, 1800) == overlay[1, 1]
end

@testset "BeiDou B1C gen_code! runs (BOC(1,1) LUT path)" begin
    # The BOC(1,1) bake/resample path is already covered by GPS L1C / Galileo
    # E1B; here we confirm the B1C wiring drives it to a finite, non-zero ±1
    # Int8 output (data has no secondary; pilot applies the 1800-chip overlay).
    sampling_rate = 25e6Hz
    samples = 4000
    for signal in (BeiDouB1C_D(), BeiDouB1C_P())
        buf = zeros(Int8, samples)
        gen_code!(buf, signal, 1, sampling_rate, 1_023_000Hz, 0.0, 0)
        @test all(x -> x == 1 || x == -1, buf)
        # BOC(1,1) sub-carrier is zero-mean, so the average stays well below 1.
        @test abs(sum(buf)) < samples * 0.5
    end
end
