# Cross-check the L1C code generators against the verification chunks
# published in IS-GPS-800G Table 3.2-2 (initial / final 24 chips per
# PRN) and Table 3.2-3 (initial / final 11 overlay bits per PRN). The
# constants live in src/gps/l1c_constants.jl as
# L1C_{D,P}_{INITIAL,FINAL}_24 and L1C_OVERLAY_{INIT,FINAL}_11.

# Decode a 24-chip-packed UInt32 (octal in the spec, MSB-first) into ±1 Int8.
function decode_24_chips(packed::UInt32)
    chips = Vector{Int8}(undef, 24)
    @inbounds for k = 1:24
        bit = (packed >> (24 - k)) & 0x0001
        chips[k] = bit == 0 ? Int8(1) : Int8(-1)
    end
    return chips
end

# Decode an 11-bit-packed UInt16 (octal in the spec, MSB-first) into ±1 Int8.
function decode_11_chips(packed::UInt16)
    chips = Vector{Int8}(undef, 11)
    @inbounds for k = 1:11
        bit = (packed >> (11 - k)) & 0x0001
        chips[k] = bit == 0 ? Int8(1) : Int8(-1)
    end
    return chips
end

@testset "L1C primary code generator (IS-GPS-800G §3.2.2.1.1)" begin
    # The constructors of GPSL1C_D and GPSL1C_P build the full
    # 10230 × 63 primary code matrix; we check both endpoints (first 24
    # and last 24 chips) for every PRN against the spec's published
    # values.
    @testset "L1C-D PRN $prn" for prn = 1:GNSSSignals.L1C_NUM_PRNS
        sig = GPSL1C_D()
        # Codes are stored Int16; the spec's expected chips are Int8.
        first_24 = Int8.(sig.codes[1:24, prn])
        last_24  = Int8.(sig.codes[end-23:end, prn])
        @test first_24 == decode_24_chips(GNSSSignals.L1C_D_INITIAL_24[prn])
        @test last_24  == decode_24_chips(GNSSSignals.L1C_D_FINAL_24[prn])
    end

    @testset "L1C-P PRN $prn" for prn = 1:GNSSSignals.L1C_NUM_PRNS
        sig = GPSL1C_P()
        first_24 = Int8.(sig.codes[1:24, prn])
        last_24  = Int8.(sig.codes[end-23:end, prn])
        @test first_24 == decode_24_chips(GNSSSignals.L1C_P_INITIAL_24[prn])
        @test last_24  == decode_24_chips(GNSSSignals.L1C_P_FINAL_24[prn])
    end
end

@testset "L1C overlay code generator (IS-GPS-800G §3.2.2.1.2)" begin
    sig = GPSL1C_P()
    overlay = sig.overlay_codes
    @testset "Overlay PRN $prn" for prn = 1:GNSSSignals.L1C_NUM_PRNS
        first_11 = overlay[1:11, prn]
        last_11  = overlay[end-10:end, prn]
        @test first_11 == decode_11_chips(GNSSSignals.L1C_OVERLAY_INIT_11[prn])
        @test last_11  == decode_11_chips(GNSSSignals.L1C_OVERLAY_FINAL_11[prn])
    end
end
