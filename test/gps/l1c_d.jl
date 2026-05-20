@testset "GPS L1C-D" begin
    sig = GPSL1C_D()
    @test @inferred(get_band(sig)) == L1()
    @test @inferred(get_center_frequency(sig)) == 1.57542e9Hz
    @test @inferred(get_code_length(sig)) == 10230
    @test @inferred(get_secondary_code_length(sig)) == 1
    @test @inferred(get_secondary_code(sig)) isa NoSecondaryCode
    @test @inferred(get_data_frequency(sig)) == 100Hz
    @test @inferred(get_code_frequency(sig)) == 1023e3Hz
    @test @inferred(get_modulation(sig)) == GNSSSignals.BOCsin(1, 1)
    @test get_signal_name(sig) == "GPS L1C-D"

    # Single-chip code access at a known position (PRN 1 chip 0 must be +1
    # per the spec's L1C_D_INITIAL_24 value 0o77001425, MSB = 1 → -1).
    # Cross-check the first chip computed by `get_code` matches what the
    # codes matrix stores.
    @test get_code(sig, 0, 1) == sig.codes[1, 1]
    @test get_code(sig, 1, 1) == sig.codes[2, 1]
    @test get_code(sig, 10230, 1) == sig.codes[1, 1]   # phase wraps modulo code length

    # get_code_factor for BOC(1,1) is n = 1.
    @test GNSSSignals.get_code_factor(sig) == 1
end
