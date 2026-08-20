@testset "Galileo E5a-QP" begin
    e5a_qp = GalileoE5aQP()
    @test @inferred(get_band(e5a_qp)) == L5()
    @test @inferred(get_center_frequency(e5a_qp)) == 1.17645e9Hz
    @test @inferred(get_code_length(e5a_qp)) == 330
    @test @inferred(get_secondary_code_length(e5a_qp)) == 1
    @test @inferred(get_secondary_code(e5a_qp)) isa NoSecondaryCode
    @test @inferred(get_data_frequency(e5a_qp)) == 0Hz      # no user data
    @test @inferred(get_code_frequency(e5a_qp)) == 5115e3Hz
    @test get_signal_name(e5a_qp) == "Galileo E5a-QP"
    @test @inferred(get_modulation(e5a_qp)) == GNSSSignals.LOC()
    @test get_code_type(e5a_qp) === Int16
    # 40 codes, not 50: the ICD caps the family at 40 for cross-correlation.
    @test size(get_codes(e5a_qp)) == (330, 40)

    # Code 1 starts 0110 1011 ... (first hex symbol 6 of Annex C.9 code 1).
    @test @inferred(get_code(e5a_qp, 0, 1)) == -1
    @test @inferred(get_code(e5a_qp, 1, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(e5a_qp, 0.0, 1)) == -1
    # 330 chips = 2/31 ms; the code simply repeats, with no overlay.
    @test get_code(e5a_qp, 330, 1) == get_code(e5a_qp, 0, 1)

    # BPSK(5) spectrum, and 5.115 MHz against the 1176.45 MHz E5a carrier.
    @test GNSSSignals.get_code_factor(e5a_qp) == 1
    @test get_code_spectrum(e5a_qp, 0) ≈ 1.0Hz / get_code_frequency(e5a_qp)
    @test get_code_center_frequency_ratio(e5a_qp) ≈ 1 / 230
end

@testset "Galileo E5a-QP primary codes match the ICD electronic annex" begin
    # `GalileoE5aQP` builds each code by XOR-ing its three generative short codes
    # (OS SIS ICD v2.2 §3.4.2, Table 21) repeated up to 330 chips. The fixture
    # holds all 40 codes as listed in the electronic annex C.9 — the same codes
    # written out chip by chip — so this checks the construction end to end.
    ref = _load_packed_hex_fixture(
        joinpath(@__DIR__, "fixtures", "galileo_e5a_qp_primary.hex.gz"),
        330 * 40,
    )
    @test get_codes(GalileoE5aQP()) == reshape(ref, 330, 40)

    # Each code is built from a 15-, an 11- and a 10-chip short code, so its
    # generative codes repeat 22, 30 and 33 times within the 330 chips.
    @test 330 == 22 * 15 == 30 * 11 == 33 * 10
end

@testset "Galileo E5a-QP code generation" begin
    e5a_qp = GalileoE5aQP()
    # A 330-chip code is shorter than most buffers: generating several periods
    # exercises the wrap of a very short code.
    n = 3 * 330 + 7
    buffer = zeros(Int8, n)
    gen_code!(buffer, e5a_qp, 3, 5115e3Hz)
    @test buffer == [get_code(e5a_qp, i, 3) for i = 0:n-1]
end
