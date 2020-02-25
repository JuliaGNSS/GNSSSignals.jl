@testset "Fast cis" begin
    @test @inferred(GNSSSignals.sin_fast(π / 2)) ≈ sin(π / 2)
    @test @inferred(GNSSSignals.sin_fast(-π / 2)) ≈ sin(-π / 2)
    @test @inferred(GNSSSignals.sin_fast(0.0)) ≈ sin(0.0)
    @test @inferred(GNSSSignals.sin_fast(1π)) ≈ sin(π) atol = eps(Float64)
    @test @inferred(GNSSSignals.sin_fast(-1π)) ≈ sin(-π) atol = eps(Float64)

    @test @inferred(GNSSSignals.cos_fast(π / 2)) ≈ cos(π / 2) atol = eps(Float64)
    @test @inferred(GNSSSignals.cos_fast(-π / 2)) ≈ cos(-π / 2) atol = eps(Float64)
    @test @inferred(GNSSSignals.cos_fast(0.0)) ≈ cos(0.0)
    @test @inferred(GNSSSignals.cos_fast(1π)) ≈ cos(π)
    @test @inferred(GNSSSignals.cos_fast(-1π)) ≈ cos(-π)

    x = -π:0.01:π
    @test sqrt(sum(abs2.(GNSSSignals.sin_vfast.(x) .- sin.(x))) / length(x)) < 0.0359
    @test sqrt(sum(abs2.(GNSSSignals.sin_fast.(x) .- sin.(x))) / length(x)) < 0.000597

    @test sqrt(sum(abs2.(GNSSSignals.cos_vfast.(x) .- cos.(x))) / length(x)) < 0.03589
    @test sqrt(sum(abs2.(GNSSSignals.cos_fast.(x) .- cos.(x))) / length(x)) < 0.000597

    @test sqrt(sum(abs2.(GNSSSignals.cis_vfast.(x) .- cis.(x))) / length(x)) < 0.0507
    @test sqrt(sum(abs2.(GNSSSignals.cis_fast.(x) .- cis.(x))) / length(x)) < 0.00085

    @test get_carrier_fast_unsafe.(x) == GNSSSignals.cis_fast.(x)
    @test get_carrier_vfast_unsafe.(x) == GNSSSignals.cis_vfast.(x)
end

@testset "Generate carrier" begin

    carrier = StructArray(zeros(Complex{Int16}, 2500))

    fpcarrier!(
        carrier,
        1500Hz,
        2.5e6Hz,
        0.25, # π / 2
        start_sample = 111,
        num_samples = 2390,
        bits = Val(7)
    )

    @test sqrt(mean(abs2.(carrier.re[111:2500] ./ 1 << 7 .-
                            cos.(2π * (0:2389) * 1500Hz / 2.5e6Hz .+ π / 2)))) < 6.5e-3
    @test sqrt(mean(abs2.(carrier.im[111:2500] ./ 1 << 7 .-
                            sin.(2π * (0:2389) * 1500Hz / 2.5e6Hz .+ π / 2)))) < 6.5e-3
end
