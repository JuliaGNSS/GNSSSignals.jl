@testset "Shift register" begin
    registers = 8191
    for i = 1:8191
        output_xb, registers = @inferred GNSSSignals.shift_register(
            registers, [1, 3, 4, 6, 7, 8, 12, 13]
        )
        results = [2788, 2056 , 3322, 2087, 6431]
        if (i in [266, 804, 1559, 3471, 5343])
            @test registers in results
        end
    end
    @test registers == 8191
end

@testset "GPS L1" begin

    gpsl5 = GPSL5()
    @test @inferred(get_center_frequency(gpsl5)) == 1.17645e9Hz
    @test @inferred(get_code_length(gpsl5)) == 10230
    @test @inferred(get_secondary_code_length(gpsl5)) == 10
    @test @inferred(get_code(gpsl5, 0, 1)) == 1
    @test @inferred(get_code(gpsl5, 0.0, 1)) == 1
    @test @inferred(get_code_unsafe(gpsl5, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl5)) == 100Hz
    @test @inferred(get_code_frequency(gpsl5)) == 10230e3Hz
    @test get_code.(gpsl5, 0:10229, 1) == L5_SAT1_CODE

end

@testset "Neuman sequence" begin
    gpsl5 = GPSL5()
    code = get_code.(gpsl5, 0:103199, 1)
    satellite_code = code[1:10230]
    neuman_hofman_code = [0,0,0,0,1,1,0,1,0,1]
    for i = 1:10
        @test code[1+10230*(i-1):10230*i] == (
            satellite_code .* (Int8(-1)^neuman_hofman_code[i])
        )
    end
    @test code[1:10230] == code[10231:20460]
end
