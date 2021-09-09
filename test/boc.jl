@testset "BOCcos" begin
    @testset "BOCCos($(get_system_string(system)), 0, $n)" for system in [GPSL1()], n in [1]#[GPSL1(), GPSL5(), GalileoE1B()], n in [1, 5, 10]
        boc = BOCcos(system, 0, n)
        @test @inferred(get_center_frequency(boc)) == get_center_frequency(system)
        @test @inferred(get_code_length(boc)) == get_code_length(system)
        @test @inferred(get_secondary_code_length(boc)) == get_secondary_code_length(system)
        @test @inferred(get_code(boc, 0, 1)) == get_code(system, 0, 1)
        @test @inferred(get_code(boc, 0.0, 1)) == get_code(system, 0.0, 1)
        @test @inferred(get_data_frequency(boc)) == get_data_frequency(system)
        @test @inferred(get_code_frequency(boc)) == n * get_code_frequency(system)
        @test get_code.(boc, 0:1022, 1) == get_code.(system, 0:1022, 1)
        @test GNSSSignals.get_code_unsafe.(boc, 0:1022, 1) == get_code.(system, 0:1022, 1)
        @test GNSSSignals.get_code_unsafe.(boc, 0.0:1022.0, 1) == get_code.(system, 0.0:1022.0, 1)
    end

    @testset "BOCcos(GPSL1, $m, 1) modulation" for m in [1,2,2.5,5,12.5]
        gpsl1 = GPSL1()
        rate = 1e-3
        tau = 0:rate:10
        boc_code = get_code.(BOCcos(gpsl1, m, 1), tau, 1) ./ get_code.(gpsl1, tau, 1)
        idx = findall(diff(boc_code) .== 2)
        len = mean(diff(idx))
        @test isapprox(len, 1 / (m * rate), atol = 1 / length(idx))
    end
end
