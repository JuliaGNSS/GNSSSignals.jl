@testset "Galileo E1B" begin

    @test @inferred(get_center_frequency(GalileoE1B)) == 1.57542e9Hz
    @test @inferred(get_code_length(GalileoE1B)) == 4092
    @test @inferred(get_shortest_code_length(GalileoE1B)) == 4092
    @test @inferred(get_code(GalileoE1B, 0, 1)) == 1
    @test @inferred(get_code(GalileoE1B, 0.0, 1)) == 1
    @test @inferred(get_code(GalileoE1B, 0.5, 1)) == -1
    @test @inferred(get_code(GalileoE1B, 1.0, 1)) == 1
    @test @inferred(get_code(GalileoE1B, 1.5, 1)) == -1
    @test @inferred(get_code_unsafe(GalileoE1B, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(GalileoE1B)) == 250Hz
    @test @inferred(get_code_frequency(GalileoE1B)) == 1023e3Hz

end
