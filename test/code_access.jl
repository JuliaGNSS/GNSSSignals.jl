@testset "get_code returns get_code_type for $(signal_label(S))" for S in ALL_SIGNALS
    signal = make_signal(S)
    @test eltype(get_codes(signal)) === Int8
    @test @inferred(get_code(signal, 0.3, 1)) isa get_code_type(signal)
    @test @inferred(get_code_unsafe(signal, 0.3, 1)) isa get_code_type(signal)
end

@testset "get_code values correlate without overflowing" begin
    # `reduce(+, …)` keeps the element type, unlike `sum`; on the stored Int8 chips it would
    # wrap around at 127.
    code = get_code.(GPSL1CA(), 0:1022, 1)
    @test reduce(+, code .* code) == 1023
end

@testset "Broadcasting GNSS signals" begin
    gpsl1ca = GPSL1CA()
    # Test that signals can be broadcast
    phases = [0.0, 1.0, 2.0]
    result = get_code.(gpsl1ca, phases, 1)
    @test length(result) == 3
    @test all(x -> x ∈ [-1, 1], result)
end
