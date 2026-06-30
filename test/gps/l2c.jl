# Independent GPS L2 civil (L2 CM / L2 CL) reference generator.
#
# This re-derives the L2C codes straight from the IS-GPS-200N §3.3.2.4 degree-27
# modular shift-register generator and, crucially, also returns the register's
# *end state* so each PRN can be pinned to the "End Shift Register State"
# published in IS-GPS-200N Tables 3-IIa/3-IIb. Those end states are external
# ground truth (transcribed from the ICD below), independent of the per-PRN
# initial-state table under test in `src/gps/l2c_constants.jl`: stepping the
# initial state forward and landing on the published end state validates the
# initial state, the feedback polynomial, and the iteration count all at once.

# Feedback mask derived independently from the ICD polynomial (1112225171
# octal), i.e. the polynomial with its X^0 term dropped (`poly >> 1`).
const L2C_REF_FEEDBACK_MASK = UInt32(0o1112225171 >> 1)

# Generate one L2C code from `initial_state`, returning the ±1 chip vector and
# the register state that produced the final (length-th) chip — the value the
# ICD tabulates as the End Shift Register State.
function _reference_l2c_code(initial_state, code_length)
    register = UInt32(initial_state)
    code = Vector{Int8}(undef, code_length)
    end_state = register
    for i = 1:code_length
        output = register & 0x1
        code[i] = output == 0 ? Int8(1) : Int8(-1)
        end_state = register
        register = (register >> 1) ⊻ (L2C_REF_FEEDBACK_MASK * output)
    end
    return code, end_state
end

# End Shift Register States (octal) from IS-GPS-200N Tables 3-IIa/3-IIb,
# "L2 CM *" column (short-cycled period = 10230 chips), PRNs 1-63.
const L2CM_END_STATES = UInt32[
    0o552566002, 0o034445034, 0o723443711, 0o511222013, 0o463055213,
    0o667044524, 0o652322653, 0o505703344, 0o520302775, 0o244205506,
    0o236174002, 0o654305531, 0o435070571, 0o630431251, 0o234043417,
    0o535540745, 0o043056734, 0o731304103, 0o412120105, 0o365636111,
    0o143324657, 0o110766462, 0o602405203, 0o177735650, 0o630177560,
    0o653467107, 0o406576630, 0o221777100, 0o773266673, 0o100010710,
    0o431037132, 0o624127475, 0o154624012, 0o275636742, 0o644341556,
    0o514260662, 0o133501670, 0o453413162, 0o637760505, 0o612775765,
    0o136315217, 0o264252240, 0o113027466, 0o774524245, 0o161633757,
    0o603442167, 0o213146546, 0o721323277, 0o207073253, 0o130632332,
    0o606370621, 0o330610170, 0o744312067, 0o154235152, 0o525024652,
    0o535207413, 0o655375733, 0o316666241, 0o525453337, 0o114323414,
    0o755234667, 0o526032633, 0o602375063,
]

# End Shift Register States (octal) from IS-GPS-200N Tables 3-IIa/3-IIb,
# "L2 CL **" column (short-cycled period = 767250 chips), PRNs 1-63.
const L2CL_END_STATES = UInt32[
    0o267724236, 0o167516066, 0o771756405, 0o047202624, 0o052770433,
    0o761743665, 0o133015726, 0o610611511, 0o352150323, 0o051266046,
    0o305611373, 0o504676773, 0o272572634, 0o731320771, 0o631326563,
    0o231516360, 0o030367366, 0o713543613, 0o232674654, 0o641733155,
    0o730125345, 0o000316074, 0o171313614, 0o001523662, 0o023457250,
    0o330733254, 0o625055726, 0o476524061, 0o602066031, 0o012412526,
    0o705144501, 0o615373171, 0o041637664, 0o100107264, 0o634251723,
    0o257012032, 0o703702423, 0o463624741, 0o673421367, 0o703006075,
    0o746566507, 0o444022714, 0o136645570, 0o645752300, 0o656113341,
    0o015705106, 0o002757466, 0o100273370, 0o304463615, 0o054341657,
    0o333276704, 0o750231416, 0o541445326, 0o316216573, 0o007360406,
    0o112114774, 0o042303316, 0o353150521, 0o044511154, 0o244410144,
    0o562324657, 0o027501534, 0o521240373,
]

@testset "GPS L2C feedback polynomial" begin
    # IS-GPS-200N §3.3.2.4: maximal polynomial 1112225171 (octal), degree 27.
    @test GNSSSignals.GPS_L2C_CODE_POLYNOMIAL == 0o1112225171
    # The modular-LFSR feedback mask is poly >> 1, equal to the binary mask
    # used by PocketSDR (`0b100100101001001010100111100`).
    @test GNSSSignals.GPS_L2C_FEEDBACK_MASK == UInt32(0o1112225171 >> 1)
    @test GNSSSignals.GPS_L2C_FEEDBACK_MASK == 0b100100101001001010100111100
end

@testset "GPS L2 CM" begin
    gpsl2cm = GPSL2CM()
    @test @inferred(get_band(gpsl2cm)) == L2()
    @test @inferred(get_center_frequency(gpsl2cm)) == 1.2276e9Hz
    @test @inferred(get_code_length(gpsl2cm)) == 10230
    @test @inferred(get_secondary_code_length(gpsl2cm)) == 1
    @test @inferred(get_secondary_code(gpsl2cm)) isa GNSSSignals.NoSecondaryCode
    @test @inferred(get_code(gpsl2cm, 0, 1)) == 1
    @test @inferred(get_code(gpsl2cm, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(gpsl2cm, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl2cm)) == 50Hz
    @test @inferred(get_code_frequency(gpsl2cm)) == 511.5e3Hz
    @test get_signal_name(gpsl2cm) == "GPS L2CM"
    @test @inferred(get_modulation(gpsl2cm)) == GNSSSignals.LOC()
    @test get_code_type(gpsl2cm) === Int16
    @test GNSSSignals.get_code_factor(gpsl2cm) == 1

    # BPSK(0.5115) spectrum: nulls at multiples of the chip rate, unit integral.
    @test get_code_spectrum(gpsl2cm, 0) ≈ 1.0Hz / get_code_frequency(gpsl2cm)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(gpsl2cm, m * get_code_frequency(gpsl2cm)) == 0
        @test get_code_spectrum(gpsl2cm, -m * get_code_frequency(gpsl2cm)) == 0
    end
    @test sum(get_code_spectrum.(gpsl2cm, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5

    # Code-to-center-frequency ratio: 511.5 kHz / 1227.6 MHz = 1 / 2400.
    @test get_code_center_frequency_ratio(gpsl2cm) ≈ 1 / 2400
end

@testset "GPS L2 CL" begin
    gpsl2cl = GPSL2CL()
    @test @inferred(get_band(gpsl2cl)) == L2()
    @test @inferred(get_center_frequency(gpsl2cl)) == 1.2276e9Hz
    @test @inferred(get_code_length(gpsl2cl)) == 767250
    @test @inferred(get_secondary_code_length(gpsl2cl)) == 1
    @test @inferred(get_secondary_code(gpsl2cl)) isa GNSSSignals.NoSecondaryCode
    @test @inferred(get_code(gpsl2cl, 0, 1)) == 1
    @test @inferred(get_code(gpsl2cl, 0.0, 1)) == 1
    @test @inferred(GNSSSignals.get_code_unsafe(gpsl2cl, 0.0, 1)) == 1
    @test @inferred(get_data_frequency(gpsl2cl)) == 0Hz   # dataless pilot
    @test @inferred(get_code_frequency(gpsl2cl)) == 511.5e3Hz
    @test get_signal_name(gpsl2cl) == "GPS L2CL"
    @test @inferred(get_modulation(gpsl2cl)) == GNSSSignals.LOC()
    @test get_code_type(gpsl2cl) === Int16
    @test GNSSSignals.get_code_factor(gpsl2cl) == 1

    @test get_code_spectrum(gpsl2cl, 0) ≈ 1.0Hz / get_code_frequency(gpsl2cl)
    @testset "Test $(m). zero" for m = 1:10
        @test get_code_spectrum(gpsl2cl, m * get_code_frequency(gpsl2cl)) == 0
        @test get_code_spectrum(gpsl2cl, -m * get_code_frequency(gpsl2cl)) == 0
    end
    @test sum(get_code_spectrum.(gpsl2cl, -1e12:1e4:1e12)) * 1e4 ≈ 1 rtol = 1e-5

    @test get_code_center_frequency_ratio(gpsl2cl) ≈ 1 / 2400
end

@testset "GPS L2 CM codes match the IS-GPS-200N reference" begin
    # The reference generator's end state is pinned to the IS-GPS-200N End
    # Shift Register State for every PRN, then the package's stored code matrix
    # is cross-checked against that reference, column by column.
    gpsl2cm = GPSL2CM()
    for prn = 1:63
        ref_code, end_state =
            _reference_l2c_code(GNSSSignals.GPS_L2CM_INITIAL_STATES[prn], 10230)
        @test end_state == L2CM_END_STATES[prn]
        @test get_codes(gpsl2cm)[:, prn] == ref_code
    end
end

@testset "GPS L2 CL codes match the IS-GPS-200N reference" begin
    gpsl2cl = GPSL2CL()
    for prn = 1:63
        ref_code, end_state =
            _reference_l2c_code(GNSSSignals.GPS_L2CL_INITIAL_STATES[prn], 767250)
        @test end_state == L2CL_END_STATES[prn]
        @test get_codes(gpsl2cl)[:, prn] == ref_code
    end
end

@testset "GPS L2 CM and CL are distinct, balanced codes" begin
    # The CM and CL codes share a generator but differ in initial state, so the
    # 10230-chip windows must differ; CM is an exactly balanced ±1 sequence.
    gpsl2cm = GPSL2CM()
    gpsl2cl = GPSL2CL()
    for prn in (1, 32, 63)
        cm = get_codes(gpsl2cm)[:, prn]
        cl_head = get_codes(gpsl2cl)[1:10230, prn]
        @test cm != cl_head
        @test abs(sum(Int.(cm))) <= 10230 ÷ 100
    end
end
