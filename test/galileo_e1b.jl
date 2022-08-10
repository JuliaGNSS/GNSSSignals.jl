@testset "Galileo E1B" begin

    galileo_e1b = GalileoE1B()
    @test @inferred(get_center_frequency(galileo_e1b)) == 1.57542e9Hz
    @test @inferred(get_code_length(galileo_e1b)) == 4092
    @test @inferred(get_secondary_code_length(galileo_e1b)) == 1
    @test @inferred(get_code(galileo_e1b, 0, 1)) ≈ 1 * sqrt(10/11) + 1 * sqrt(1/11)
    @test @inferred(get_code(galileo_e1b, 0.0, 1)) ≈ 1 * sqrt(10/11) + 1 * sqrt(1/11)
    @test @inferred(get_code(galileo_e1b, 0.5, 1)) ≈ -1 * sqrt(10/11) + 1 * sqrt(1/11)
    @test @inferred(get_code(galileo_e1b, 1.0, 1)) ≈ 1 * sqrt(10/11) + 1 * sqrt(1/11)
    @test @inferred(get_code(galileo_e1b, 1.5, 1)) ≈ -1 * sqrt(10/11) + 1 * sqrt(1/11)
    @test @inferred(GNSSSignals.get_code_unsafe(galileo_e1b, 0.0, 1)) ≈ 1 * sqrt(10/11) + 1 * sqrt(1/11)
    @test @inferred(get_data_frequency(galileo_e1b)) == 250Hz
    @test @inferred(get_code_frequency(galileo_e1b)) == 1023e3Hz

    @test GNSSSignals.get_code_factor(galileo_e1b) == 1
end
