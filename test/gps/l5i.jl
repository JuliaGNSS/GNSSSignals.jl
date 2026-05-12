@testset "Shift register" begin
    registers = 8191
    for i = 1:8191
        output_xb, registers =
            @inferred GNSSSignals.shift_register(registers, [1, 3, 4, 6, 7, 8, 12, 13])
        results = [2788, 2056, 3322, 2087, 6431]
        if (i in [266, 804, 1559, 3471, 5343])
            @test registers in results
        end
    end
    @test registers == 8191
end

@testset "GPS L5-I" begin
    gpsl5i = GPSL5I()
    @test @inferred(get_band(gpsl5i)) == L5()
    @test @inferred(get_center_frequency(gpsl5i)) == 1.17645e9Hz
    @test @inferred(get_code_length(gpsl5i)) == 10230
    @test @inferred(get_secondary_code_length(gpsl5i)) == 10
    @test @inferred(get_secondary_code(gpsl5i)) isa SharedSecondaryCode{10}
    @test @inferred(get_code(gpsl5i, 0, 1)) == 1
    @test @inferred(get_code(gpsl5i, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(gpsl5i, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl5i)) == 100Hz
    @test @inferred(get_code_frequency(gpsl5i)) == 10230e3Hz
    @test get_code.(gpsl5i, 0:10229, 1) == L5_SAT1_CODE
    @test get_signal_name(gpsl5i) == "GPS L5-I"

    @test GNSSSignals.get_code_factor(gpsl5i) == 1

    @test get_code_spectrum(gpsl5i, 0) ≈ 1.0Hz / get_code_frequency(gpsl5i)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(gpsl5i, m * get_code_frequency(gpsl5i)) == 0
        @test get_code_spectrum(gpsl5i, -m * get_code_frequency(gpsl5i)) == 0
    end
    @test sum(get_code_spectrum.(gpsl5i, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5

    @test get_code_center_frequency_ratio(gpsl5i) ≈ 1 / 115
end

@testset "Neuman sequence" begin
    gpsl5i = GPSL5I()
    code = get_code.(gpsl5i, 0:103199, 1)
    satellite_code = code[1:10230]
    neuman_hofman_code = [0, 0, 0, 0, 1, 1, 0, 1, 0, 1]
    for i = 1:10
        @test code[1+10230*(i-1):10230*i] ==
              (satellite_code .* (Int8(-1)^neuman_hofman_code[i]))
    end
    @test code[1:10230] == code[10231:20460]
end
