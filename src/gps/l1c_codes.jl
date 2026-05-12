# Internal helpers for building the GPS L1C primary and overlay codes.
# The output of these helpers is fed into the `GPSL1C_D` and `GPSL1C_P`
# constructors (in l1c_d.jl and l1c_p.jl).
#
# Algorithms follow IS-GPS-800G §3.2.2.1.1 (ranging code) and §3.2.2.1.2
# (overlay code). IRN-IS-800J-003 makes no changes here. The tests in
# test/gps/l1c_codes.jl verify every PRN's initial 24 chips,
# final 24 chips, and overlay final 11 bits against the published spec
# values.

"""
$(SIGNATURES)

Generate the length-10223 Legendre sequence used by the L1C Weil-code
construction.

Per IS-GPS-800G §3.2.2.1.1:
- `L(0) = 0`
- `L(t) = 1` if `t` is a quadratic residue mod 10223 (i.e. there exists
  an integer `x` such that `t ≡ x² (mod 10223)`), else 0.

Returns an `AbstractVector{Bool}` of length 10223.
"""
function _l1c_legendre_sequence()
    p = L1C_WEIL_LENGTH
    ls = falses(p)
    # Quadratic residues mod p: x² for x = 1..(p-1)/2. (x and -x have the
    # same square mod p, so we only need the first half.)
    for x = 1:((p - 1) ÷ 2)
        ls[mod(x * x, p) + 1] = true
    end
    # L(0) = 0 is already the default (false).
    return ls
end

"""
$(SIGNATURES)

Build the per-PRN Weil-code XOR'd with the 7-chip expansion sequence,
returning a length-10230 `Vector{Int8}` of ±1 chips for one PRN.

Used by both `GPSL1C_D` and `GPSL1C_P`; the only difference is which
Weil index `w` and insertion index `p` are provided.

`legendre` is the precomputed sequence from [`_l1c_legendre_sequence`](@ref);
caller passes it in to share the cost across all 63 PRNs.
"""
function _l1c_primary_code(legendre::AbstractVector{Bool}, weil_index::Integer, insertion_index::Integer)
    p = L1C_WEIL_LENGTH                    # 10223
    n = L1C_PRIMARY_LENGTH                 # 10230
    ins = insertion_index                  # 1-based per IS-GPS-800G

    out = Vector{Int8}(undef, n)

    # Output positions t = 0..ins-2 (1-based 1..ins-1): Weil(t) = L(t) ⊕ L((t+w) mod p)
    @inbounds for t = 0:(ins - 2)
        out[t + 1] = (legendre[t + 1] ⊻ legendre[mod(t + weil_index, p) + 1]) ?
                     Int8(-1) : Int8(1)
    end

    # Output positions t = ins-1 .. ins+5 (1-based ins..ins+6): 7-chip expansion 0110100
    # In the spec's ±1 mapping, 0 → +1 and 1 → -1.
    @inbounds for (k, bit) in enumerate(L1C_EXPANSION_SEQUENCE)
        out[ins + k - 1] = bit == 0 ? Int8(1) : Int8(-1)
    end

    # Output positions t = ins+6..n-1 (1-based ins+7..n): Weil(t-7)
    @inbounds for t = (ins + 6):(n - 1)
        s = t - 7                          # Weil source index, 0-based
        out[t + 1] = (legendre[s + 1] ⊻ legendre[mod(s + weil_index, p) + 1]) ?
                     Int8(-1) : Int8(1)
    end

    return out
end

"""
$(SIGNATURES)

Build the L1C primary-code matrix `(L1C_PRIMARY_LENGTH, L1C_NUM_PRNS)` of
`Int8` ±1 chips, given per-PRN Weil and insertion index arrays. Used by
the `GPSL1C_D` and `GPSL1C_P` constructors.
"""
function _l1c_build_primary_codes(weil_indices::AbstractVector{<:Integer},
                                  insertion_indices::AbstractVector{<:Integer})
    @assert length(weil_indices) == L1C_NUM_PRNS
    @assert length(insertion_indices) == L1C_NUM_PRNS
    legendre = _l1c_legendre_sequence()
    codes = Matrix{Int8}(undef, L1C_PRIMARY_LENGTH, L1C_NUM_PRNS)
    for prn = 1:L1C_NUM_PRNS
        codes[:, prn] = _l1c_primary_code(legendre, weil_indices[prn],
                                          insertion_indices[prn])
    end
    return codes
end

"""
$(SIGNATURES)

Build the L1C-P overlay-code matrix `(L1C_OVERLAY_LENGTH, L1C_NUM_PRNS)`
of `Int8` ±1 chips. Each PRN's overlay is the first 1800 bits of the
output of an 11-stage LFSR with per-PRN polynomial `mij` and initial
state `init_11`, both taken from IS-GPS-800G Table 3.2-3.

The LFSR output convention used here is verified against the spec's
published final-11-bits per PRN.
"""
function _l1c_build_overlay_codes()
    n = L1C_OVERLAY_LENGTH                 # 1800
    codes = Matrix{Int8}(undef, n, L1C_NUM_PRNS)
    for prn = 1:L1C_NUM_PRNS
        # mij is stored as a 12-bit value (MSB-first). The lowest-order
        # coefficient m_{i,0} is dropped (always implicit in the LFSR
        # convention), leaving 11 tap bits at positions 11..1 of the
        # 12-bit field.
        #
        # We use a left-shift register with the MSB (bit 10) as the
        # output. The 11-bit state is stored MSB-first in the upper
        # bits of a UInt16: bit 10 = first output, bit 0 = last output
        # before any feedback wraps in. This matches the spec convention
        # where the "initial 11 bits of the sequence" are the first 11
        # outputs, with the MSB of init_11 emitted first.
        taps  = (L1C_OVERLAY_MIJ[prn]      >> 1) & 0x07FF
        state =  L1C_OVERLAY_INIT_11[prn]   & 0x07FF
        for i = 1:n
            # Output the MSB of the 11-bit register.
            codes[i, prn] = (state & 0x0400) == 0 ? Int8(1) : Int8(-1)
            # Feedback bit: parity of (state AND taps). Each tap selects
            # a register position whose contribution is XOR-summed.
            fb = count_ones(state & taps) & 0x0001
            # Left-shift the 11-bit register, inject feedback at LSB,
            # mask to keep 11 bits.
            state = ((state << 1) | UInt16(fb)) & 0x07FF
        end
    end
    return codes
end
