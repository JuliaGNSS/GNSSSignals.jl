# Shared spreading-code generators for the BeiDou Open Service signals.
#
# Two constructions cover every OS ranging code:
#
#   * `_beidou_gold_code` — a Gold code formed by the modulo-2 sum of two
#     Fibonacci LFSRs. Used by B1I (two 11-stage registers), and by B3I, B2a
#     and B2b (two 13-stage registers). Verified against the ICDs' own
#     first/last-24-chip octal values (B2a, B2b) and the B3I shift-number
#     table; B1I's two polynomials are the primitive length-2047 m-sequence
#     generators from the ICD.
#
#   * `_beidou_weil_code` — a Weil code obtained by cyclically truncating the
#     modulo-2 sum of a Legendre sequence and a shifted copy of itself. Used
#     by the B1C primary codes and the B1C/B2a pilot secondary codes. Verified
#     against the ICDs' first/last-24-chip octal values.
#
# LFSR convention (matches BDS ICD figures and reproduces every published
# verification vector): the register stages are numbered 1..n left to right;
# an output/feedback tap "k" reads stage k; on each clock the register shifts
# toward higher-numbered stages (stage_i ← stage_{i-1}) and the feedback
# (modulo-2 sum of the polynomial-tap stages) enters stage 1. Chips are mapped
# 0 → -1, 1 → +1 (`2b - 1`), matching the package-wide primary-code convention
# (see `gen_l5_code`, `read_galileo_e5a_codes`).

# Parse an ICD register-state string ("s_1 s_2 … s_n", left to right) into a
# length-n bit vector indexed by stage number.
_beidou_state(bits::AbstractString) = Int8[c == '1' ? Int8(1) : Int8(0) for c in bits]

# Advance a Fibonacci LFSR one clock: return the (pre-shift) modulo-2 sum of
# the `out_taps` stages, then shift the register toward higher stages, feeding
# the modulo-2 sum of the `fb_taps` stages into stage 1. `state` is mutated.
@inline function _beidou_lfsr_step!(state::Vector{Int8}, fb_taps, out_taps)
    n = length(state)
    out = Int8(0)
    @inbounds for t in out_taps
        out ⊻= state[t]
    end
    fb = Int8(0)
    @inbounds for t in fb_taps
        fb ⊻= state[t]
    end
    @inbounds for i = n:-1:2
        state[i] = state[i-1]
    end
    @inbounds state[1] = fb
    out
end

"""
$(SIGNATURES)

Generate one BeiDou Gold-code period of length `len` as an `Int8` vector of
±1 chips.

The code is the modulo-2 sum of two Fibonacci LFSRs of length `n`. `fb1`/`fb2`
are the feedback-tap stage numbers (the polynomial exponents `i` for which the
coefficient of `Xⁱ` is 1). `out1`/`out2` are the output-tap stage numbers whose
modulo-2 sum forms each register's output — a single stage `(n,)` for a plain
maximal-length output, or several stages for a phase-selected output (BeiDou
B1I selects the G2 phase by tapping different stages). `init1`/`init2` are the
register start states (stage 1..n).

If `reset1_at > 0`, register 1 is reset to `init1` immediately after emitting
chip `reset1_at` (BeiDou B3I/B2a/B2b short-cycle register 1 to a period of 8190
chips within each 10230-chip code period). Register 2 is never reset mid-period.
"""
function _beidou_gold_code(n, len, fb1, fb2, out1, out2, init1, init2, reset1_at)
    s1 = copy(init1)
    s2 = copy(init2)
    reset1 = copy(init1)
    code = Vector{Int8}(undef, len)
    @inbounds for i = 1:len
        o1 = _beidou_lfsr_step!(s1, fb1, out1)
        o2 = _beidou_lfsr_step!(s2, fb2, out2)
        code[i] = Int8(2) * (o1 ⊻ o2) - Int8(1)
        i == reset1_at && (s1 .= reset1)
    end
    code
end

"""
$(SIGNATURES)

Build the length-`N` Legendre sequence used by the BeiDou Weil-code
construction, as a `BitVector` indexed `L[k+1]` for `k = 0 … N-1`.

Per the BDS ICDs (e.g. BDS-SIS-ICD-B1C-1.0 Eq. 5-2): `L(0) = 0`, and for
`k > 0`, `L(k) = 1` if there exists an integer `x` with `k ≡ x² (mod N)`, else
`L(k) = 0`.
"""
function _beidou_legendre(N::Integer)
    L = falses(N)
    @inbounds for x = 1:(N-1)
        L[mod(x * x, N)+1] = true
    end
    L[1] = false        # L(0) = 0
    L
end

"""
$(SIGNATURES)

Generate a truncated BeiDou Weil code of length `N0` as an `Int8` vector of ±1
chips, from a length-`N` Legendre sequence.

Per the BDS ICDs: the Weil code is `W(k; w) = L(k) ⊕ L((k + w) mod N)`, and the
ranging code is `c(n) = W((n + p − 1) mod N; w)` for `n = 0 … N0-1`, where `w`
is the phase difference and `p` the (1-based) truncation point. Chips are
mapped 0 → -1, 1 → +1.
"""
function _beidou_weil_code(N::Integer, N0::Integer, w::Integer, p::Integer,
                           L::BitVector = _beidou_legendre(N))
    code = Vector{Int8}(undef, N0)
    @inbounds for n = 0:(N0-1)
        k = mod(n + p - 1, N)
        bit = L[k+1] ⊻ L[mod(k + w, N)+1]
        code[n+1] = Int8(2) * bit - Int8(1)
    end
    code
end

# BeiDou Neuman-Hoffman NH20 secondary code (BDS-SIS-ICD-B1I-3.0 §5.2.1 /
# BDS-SIS-ICD-B3I-1.0 §5.2.1): the 20-bit sequence 0 0 0 0 0 1 0 0 1 1 0 1 0 1
# 0 0 1 1 1 0, mapped 0 → -1, 1 → +1, modulo-2 added to the ranging code at
# 1 kHz (one chip per 1 ms primary period) to form the D1 (MEO/IGSO) navigation
# tiered code.
const BEIDOU_NH20_BITS = (0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0)

# NH20 overlay is per satellite type (BDS-SIS-ICD B1I/B3I Table 4-1): the
# MEO/IGSO satellites broadcast the D1 message with the NH20 overlay, while the
# GEO satellites broadcast the faster D2 message with NO secondary code. The two
# ICDs assign the same PRN ranges — GEO: PRN 1-5 and 59-63; MEO/IGSO: PRN 6-58.
const BEIDOU_D1_PRNS = 6:58   # MEO/IGSO satellites (D1 message, NH20 overlay)

@inline _beidou_is_d1_prn(prn::Integer) = prn in BEIDOU_D1_PRNS

# Build the (20 × num_prns) NH20 overlay matrix exposed via `PerPRNSecondaryCode`:
# the NH20 sequence (±1) for the MEO/IGSO (D1) PRNs, and an all-ones (no-op)
# column for the GEO (D2) PRNs, so a GEO PRN's tiered code is just its primary
# code — matching the ICD (GEO D2 has no NH overlay). See PocketSDR `sec_code_B1I`.
function _beidou_nh20_matrix(num_prns::Integer)
    nh = ntuple(i -> Int8(2 * BEIDOU_NH20_BITS[i] - 1), Val(20))
    m = ones(Int8, 20, num_prns)          # GEO (D2) PRNs: all-ones ⇒ no overlay
    for prn = 1:num_prns
        _beidou_is_d1_prn(prn) && (@inbounds m[:, prn] .= nh)
    end
    m
end
