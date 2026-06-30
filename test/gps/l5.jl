# Independent GPS L5 primary-code reference generator (used for both I5 and Q5).
#
# This uses a *different* formulation from the package's `gen_l5_code` (which
# free-runs the XB register from its per-PRN initial state): here the XB
# sequence is generated once over its full 8191-chip period from the all-ones
# state, and each PRN selects a starting phase via its IS-GPS-705 "XB code
# advance" (the same advances published by GNSS-SDR's `GPS_L5i_INIT_REG` /
# `GPS_L5q_INIT_REG` and PocketSDR's `L5I_XB_adv` / `L5Q_XB_adv`). Agreement
# between the two formulations validates the initial-state tables independently
# of the implementation under test. For I5, PRN 1 is additionally anchored to
# the external `L5_SAT1_CODE` reference, pinning this generator to ground truth.
const L5I_XB_ADVANCE = [    # IS-GPS-705 I5 XB code advance, PRN 1-37
    266, 365, 804, 1138, 1509, 1559, 1756, 2084, 2170, 2303,
    2527, 2687, 2930, 3471, 3940, 4132, 4332, 4924, 5343, 5443,
    5641, 5816, 5898, 5918, 5955, 6243, 6345, 6477, 6518, 6875,
    7168, 7187, 7329, 7577, 7720, 7777, 8057,
]
const L5Q_XB_ADVANCE = [    # IS-GPS-705 Q5 XB code advance, PRN 1-37
    1701, 323, 5292, 2020, 5429, 7136, 1041, 5947, 4315, 148,
    535, 1939, 5206, 5910, 3595, 5135, 6082, 6990, 3546, 1523,
    4548, 4484, 1893, 3961, 7106, 5299, 4660, 276, 4389, 3783,
    1591, 1601, 749, 1387, 1661, 3210, 708,
]

function _reference_l5_code(advance)
    xa_indices = [9, 10, 12, 13]
    xb_indices = [1, 3, 4, 6, 7, 8, 12, 13]
    # XA sequence, period 8190 (reset to all-ones one chip before its natural
    # period, exactly as in `gen_l5_code`).
    xa = Vector{Int}(undef, 10230)
    reg = 8191
    for i = 1:10230
        out, reg = GNSSSignals.shift_register(reg, xa_indices)
        xa[i] = out
        i == 8190 && (reg = 8191)
    end
    # Full XB period (8191 chips) from the all-ones state.
    xbp = Vector{Int}(undef, 8191)
    reg = 8191
    for i = 1:8191
        out, reg = GNSSSignals.shift_register(reg, xb_indices)
        xbp[i] = out
    end
    code = Vector{Int8}(undef, 10230)
    for n = 0:10229
        code[n+1] = Int8(2 * (xa[n+1] ⊻ xbp[(advance + n) % 8191 + 1]) - 1)
    end
    code
end

@testset "Shift register" begin
    registers = 8191
    for i = 1:8191
        output_xb, registers =
            @inferred GNSSSignals.shift_register(registers, [1, 3, 4, 6, 7, 8, 12, 13])
        results = [2788, 2056, 3322, 2087, 6431]
        if (i in [266, 804, 1559, 3471, 5343])
            @test registers in results
        end
    end
    @test registers == 8191
end

@testset "GPS L5-I" begin
    gpsl5i = GPSL5I()
    @test @inferred(get_band(gpsl5i)) == L5()
    @test @inferred(get_center_frequency(gpsl5i)) == 1.17645e9Hz
    @test @inferred(get_code_length(gpsl5i)) == 10230
    @test @inferred(get_secondary_code_length(gpsl5i)) == 10
    @test @inferred(get_secondary_code(gpsl5i)) isa SharedSecondaryCode{10}
    @test @inferred(get_code(gpsl5i, 0, 1)) == 1
    @test @inferred(get_code(gpsl5i, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(gpsl5i, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl5i)) == 100Hz
    @test @inferred(get_code_frequency(gpsl5i)) == 10230e3Hz
    @test get_code.(gpsl5i, 0:10229, 1) == L5_SAT1_CODE
    @test get_signal_name(gpsl5i) == "GPS L5-I"
    @test @inferred(get_modulation(gpsl5i)) == GNSSSignals.LOC()
    @test get_code_type(gpsl5i) === Int16

    @test GNSSSignals.get_code_factor(gpsl5i) == 1

    @test get_code_spectrum(gpsl5i, 0) ≈ 1.0Hz / get_code_frequency(gpsl5i)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(gpsl5i, m * get_code_frequency(gpsl5i)) == 0
        @test get_code_spectrum(gpsl5i, -m * get_code_frequency(gpsl5i)) == 0
    end
    @test sum(get_code_spectrum.(gpsl5i, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5

    @test get_code_center_frequency_ratio(gpsl5i) ≈ 1 / 115
end

@testset "GPS L5-I primary codes match the IS-GPS-705 reference" begin
    # Cross-check every supported PRN against the independent advance-offset
    # generator (a different algorithm from `gen_l5_code`). PRN 1 is also
    # pinned to L5_SAT1_CODE above, anchoring the reference to ground truth.
    gpsl5i = GPSL5I()
    for prn = 1:37
        @test get_codes(gpsl5i)[:, prn] == _reference_l5_code(L5I_XB_ADVANCE[prn])
    end
end

@testset "Neuman sequence" begin
    gpsl5i = GPSL5I()
    code = get_code.(gpsl5i, 0:103199, 1)
    satellite_code = code[1:10230]
    neuman_hofman_code = [0, 0, 0, 0, 1, 1, 0, 1, 0, 1]
    for i = 1:10
        @test code[1+10230*(i-1):10230*i] ==
              (satellite_code .* (Int8(-1)^neuman_hofman_code[i]))
    end
    @test code[1:10230] == code[10231:20460]
end

@testset "GPS L5-Q" begin
    gpsl5q = GPSL5Q()
    @test @inferred(get_band(gpsl5q)) == L5()
    @test @inferred(get_center_frequency(gpsl5q)) == 1.17645e9Hz
    @test @inferred(get_code_length(gpsl5q)) == 10230
    @test @inferred(get_secondary_code_length(gpsl5q)) == 20
    @test @inferred(get_secondary_code(gpsl5q)) isa SharedSecondaryCode{20}
    @test @inferred(get_code(gpsl5q, 0, 1)) == 1
    @test @inferred(get_code(gpsl5q, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(gpsl5q, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl5q)) == 0Hz   # dataless pilot
    @test @inferred(get_code_frequency(gpsl5q)) == 10230e3Hz
    @test get_signal_name(gpsl5q) == "GPS L5-Q"
    @test @inferred(get_modulation(gpsl5q)) == GNSSSignals.LOC()
    @test get_code_type(gpsl5q) === Int16

    # BPSK(10) spectrum and code/center-frequency ratio (Q5 shares the L5
    # band, so the ratio is the same 1/115 as GPS L5-I).
    @test GNSSSignals.get_code_factor(gpsl5q) == 1
    @test get_code_spectrum(gpsl5q, 0) ≈ 1.0Hz / get_code_frequency(gpsl5q)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(gpsl5q, m * get_code_frequency(gpsl5q)) == 0
        @test get_code_spectrum(gpsl5q, -m * get_code_frequency(gpsl5q)) == 0
    end
    @test sum(get_code_spectrum.(gpsl5q, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5
    @test get_code_center_frequency_ratio(gpsl5q) ≈ 1 / 115
end

@testset "GPS L5-Q primary codes match the IS-GPS-705 reference" begin
    # Cross-check every supported PRN against the independent advance-offset
    # generator above. This is a different algorithm from `gen_l5_code`, so
    # agreement validates the Q5 initial-state table, not just a snapshot.
    gpsl5q = GPSL5Q()
    for prn = 1:37
        @test get_codes(gpsl5q)[:, prn] == _reference_l5_code(L5Q_XB_ADVANCE[prn])
    end
end

@testset "GPS L5-Q Neuman-Hofman secondary code (NH20)" begin
    gpsl5q = GPSL5Q()
    sec = get_secondary_code(gpsl5q)
    # NH20 = 00000100110101001110 (IS-GPS-705 §3.2.1.2), 0 -> +1, 1 -> -1.
    nh20 = Int8[1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1]
    @test [GNSSSignals.secondary_value(sec, 1, k) for k = 0:19] == nh20
    # Shared across PRNs, and wraps with period 20.
    @test [GNSSSignals.secondary_value(sec, 7, k) for k = 0:19] == nh20
    @test GNSSSignals.secondary_value(sec, 1, 20) == nh20[1]
end

@testset "GPS L5-Q Neuman sequence" begin
    gpsl5q = GPSL5Q()
    code = get_code.(gpsl5q, 0:(20 * 10230 - 1), 1)
    satellite_code = code[1:10230]
    neuman_hofman_code = [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0]
    for i = 1:20
        @test code[1+10230*(i-1):10230*i] ==
              (satellite_code .* (Int8(-1)^neuman_hofman_code[i]))
    end
end
