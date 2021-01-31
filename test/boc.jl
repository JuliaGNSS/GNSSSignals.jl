@testset "BOC" begin

@testset "BOC{$t, 0, $n}" for t in [GPSL1, GPSL5, GalileoE1B], n in [1, 5, 10]
    type = BOC{t,0,n}
    @test @inferred(get_center_frequency(type)) == get_center_frequency(t)
    @test @inferred(get_code_length(type)) == get_code_length(t)
    @test @inferred(get_secondary_code_length(type)) == get_secondary_code_length(t)
    @test @inferred(get_code(type, 0, 1)) == get_code(t, 0, 1)
    @test @inferred(get_code(type, 0.0, 1)) == get_code(t, 0.0, 1)
    @test @inferred(get_code_unsafe(type, 0.0, 1)) == get_code_unsafe(t, 0.0, 1)
    @test @inferred(get_data_frequency(type)) == get_data_frequency(t)
    @test @inferred(get_code_frequency(type)) == n * 1023e3Hz
    @test get_code.(type, 0:1022, 1) == get_code.(t, 0:1022, 1)
end

@testset "BOC($m,1) modulation" for m in [1,2,2.5,5,12.5]
    rate = 1e-3
    tau = 0:rate:10
    boc = get_code.(BOC{GPSL1,m,1}, tau, 1) ./ get_code.(GPSL1, tau, 1)
    idx = findall(diff(boc) .== 2)
    len = mean(diff(idx))
    @test isapprox(len, 1 / (m * rate), atol = 1 / length(idx))
end

end