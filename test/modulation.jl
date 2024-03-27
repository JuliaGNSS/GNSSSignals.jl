@testset "Modulation" begin
    galileo_e1b = GalileoE1B()
    boc_sin = GNSSSignals.BOCsin(1, 1)
    @test GNSSSignals.get_subcarrier_code(boc_sin, 0.0) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 0.25) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 0.75) == -1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 1.25) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_sin, 1.75) == -1.0

    boc_cos = GNSSSignals.BOCcos(1, 1)
    @test GNSSSignals.get_subcarrier_code(boc_cos, 0.0) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 0.5) == -1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 1.0) == 1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 1.5) == -1.0
    @test GNSSSignals.get_subcarrier_code(boc_cos, 2.0) == 1.0

    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCcos(0, 1)
    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCcos(1, 0)
    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCsin(0, 1)
    @test_throws ErrorException("m and n must be >= 1") GNSSSignals.BOCsin(1, 0)

    @test_throws ErrorException(
        "Power of BOC1 must be between 0 and 1 and n of both BOCs must match",
    ) GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), -1.0)
    @test_throws ErrorException(
        "Power of BOC1 must be between 0 and 1 and n of both BOCs must match",
    ) GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), 2.0)
    @test_throws ErrorException(
        "Power of BOC1 must be between 0 and 1 and n of both BOCs must match",
    ) GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 2), 1 / 11)

    cboc = GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), 1 / 2)
    @test GNSSSignals.get_subcarrier_code(cboc, 0.0) ≈ sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 0.25) ≈ sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 0.75) ≈ -sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 1.25) ≈ sqrt(1 / 2) * 2
    @test GNSSSignals.get_subcarrier_code(cboc, 1.75) ≈ -sqrt(1 / 2) * 2

    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCcos(2, 1), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCcos(2, 2), 2.3) == 4
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCsin(2, 1), 2.3) == 2
    @test GNSSSignals.get_floored_phase(GNSSSignals.BOCsin(2, 2), 2.3) == 4

    @test GNSSSignals.get_code_factor(GNSSSignals.BOCcos(2, 1)) == 1
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCcos(2, 2)) == 2
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCcos(2, 2.5)) == 2.5
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCsin(2, 1)) == 1
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCsin(2, 2)) == 2
    @test GNSSSignals.get_code_factor(GNSSSignals.BOCsin(2, 2.5)) == 2.5
    @test GNSSSignals.get_code_factor(
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 1), GNSSSignals.BOCsin(1, 1), 1 / 2),
    ) == 1
    @test GNSSSignals.get_code_factor(
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 2), GNSSSignals.BOCsin(1, 2), 1 / 2),
    ) == 2
    @test GNSSSignals.get_code_factor(
        GNSSSignals.CBOC(GNSSSignals.BOCsin(1, 2.5), GNSSSignals.BOCsin(1, 2.5), 1 / 2),
    ) == 2.5
end
