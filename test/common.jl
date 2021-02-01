@testset "Common" begin

@test get_code_center_frequency_ratio(GPSL1) ≈ 1/1540

@testset "Code spectra $(system)" for system = [GPSL1, GPSL5, BOCcos{GPSL1,0,15}]
    @test get_code_spectrum(system, 0) ≈ 1.0Hz/get_code_frequency(system)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(system, m*get_code_frequency(system)) == 0
        @test get_code_spectrum(system, -m*get_code_frequency(system)) == 0
    end
    @test sum(get_code_spectrum(system, -1e12:1e4:1e12))*1e4 ≈ 1 rtol = 1e-5
end

@testset "Code spectra $(system)" for system = [GalileoE1B, BOCcos{GPSL1, 15, 2.5}]
    @test get_code_spectrum(system, 0) == 0.0
    @test sum(get_code_spectrum(system, -1e12:1e4:1e12))*1e4 ≈ 1 rtol = 1e-5
end

end
