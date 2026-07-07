# Per-PRN G2 phase-assignment taps for the BeiDou B1I ranging code.
#
# Source: BDS-SIS-ICD-B1I-3.0 (2019-02) Table 4-1, "Phase assignment of G2
# sequence". Each satellite's G2 output is the modulo-2 sum of the G2 shift
# register stages listed here (stages numbered 1..11, left to right in ICD
# Figure 4-1); XORing that with the G1 output gives the ranging code. There
# are 63 ranging codes (PRN 1..63).
const B1I_G2_PHASE_SELECT = (
    (1, 3), (1, 4), (1, 5), (1, 6), (1, 8), (1, 9),
    (1, 10), (1, 11), (2, 7), (3, 4), (3, 5), (3, 6),
    (3, 8), (3, 9), (3, 10), (3, 11), (4, 5), (4, 6),
    (4, 8), (4, 9), (4, 10), (4, 11), (5, 6), (5, 8),
    (5, 9), (5, 10), (5, 11), (6, 8), (6, 9), (6, 10),
    (6, 11), (8, 9), (8, 10), (8, 11), (9, 10), (9, 11),
    (10, 11), (1, 2, 7), (1, 3, 4), (1, 3, 6), (1, 3, 8), (1, 3, 10),
    (1, 3, 11), (1, 4, 5), (1, 4, 9), (1, 5, 6), (1, 5, 8), (1, 5, 10),
    (1, 5, 11), (1, 6, 9), (1, 8, 9), (1, 9, 10), (1, 9, 11), (2, 3, 7),
    (2, 5, 7), (2, 7, 9), (3, 4, 5), (3, 4, 9), (3, 5, 6), (3, 5, 8),
    (3, 5, 10), (3, 5, 11), (3, 6, 9),
)
