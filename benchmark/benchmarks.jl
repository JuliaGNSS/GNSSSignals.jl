using BenchmarkTools
using GNSSSignals
using Unitful: Hz, ustrip, @u_str

# Strip a Unitful frequency to a plain Float64 Hz value (for computing N).
ustrip_hz(x) = Float64(ustrip(u"Hz", x))
# Sampling frequency as an MHz string for a benchmark label ("5 MHz", "2.5 MHz").
_mhz(fs) = (m = ustrip_hz(fs) / 1e6; (m == round(m) ? string(Int(m)) : string(m)) * " MHz")

const SUITE = BenchmarkGroup()

# `gen_code!` at the signal's nominal chip rate, phase and index, its defaults, so the same
# calls run on the compared revision, whatever its chip-rate accessor is named.

# `gen_code!` fills `Int8`.
_buf(n) = zeros(Int8, n)

# ─────────────────────────────────────────────────────────────────────────────
# 2000 sample code generation — fixed 2000-sample spot-check at a couple of representative rates. A
# small, stable sanity row for the plain (GPSL1/GPSL5) and CBOC (GalileoE1B) generators; full
# per-signal coverage lives in "1 ms code generation" below.
# ─────────────────────────────────────────────────────────────────────────────
for (name, signal, fs) in (("GPSL1", GPSL1CA(), 2e6Hz), ("GPSL5", GPSL5I(), 20e6Hz),
                           ("GalileoE1B", GalileoE1B(), 15e6Hz))
    out = _buf(2000)
    SUITE["2000 sample code generation"]["$name @ $(_mhz(fs))"] = @benchmarkable gen_code!(
        $out, $signal, $1, $fs,
    ) evals = 10 samples = 10000
end

# ─────────────────────────────────────────────────────────────────────────────
# 1 ms code generation — one gen_code! fill of an N = round(Int, fs·1e-3)-sample epoch per signal, for
# EVERY signal the package defines. `fs` is chosen per (code_frequency, subchip_factor P) family,
# so signals sharing a chip rate + modulation order share a rate and timing differences reflect the
# signal/modulation rather than a different fs; each row carries its fs in the label. fs clears the
# LUT's fs ≥ code_frequency·P requirement with headroom:
#   1.023 MHz  P=1  LOC         →  5 MHz   GPSL1CA
#   1.023 MHz  P=2  BOC(1,1)    →  5 MHz   GPSL1C_D, GalileoE1B_BOC11, GalileoE1C_BOC11,
#                                          BeiDouB1C_D, BeiDouB1C_P
#   1.023 MHz  P=12 CBOC/TMBOC  → 15 MHz   GPSL1C_P, GalileoE1B, GalileoE1C
#   0.5115 MHz P=1  LOC (L2C)   →  2.5 MHz GPSL2CM, GPSL2CL
#   2.046 MHz  P=1  LOC         → 10 MHz   BeiDouB1I
#   5.115 MHz  P=1  LOC         → 20 MHz   GalileoE5aQP, GalileoE6B, GalileoE6C
#   10.23 MHz  P=1  LOC         → 40 MHz   GPSL5I, GPSL5Q, GalileoE5aI, GalileoE5aQ,
#                                          GalileoE5bI, GalileoE5bQ, BeiDouB3I, BeiDouB2bI,
#                                          BeiDouB2aI, BeiDouB2aQ
# ─────────────────────────────────────────────────────────────────────────────
const _MS_CASES = (
    (GPSL1CA, 5e6Hz), (GPSL1C_D, 5e6Hz), (GPSL1C_P, 15e6Hz),
    (GPSL2CM, 2.5e6Hz), (GPSL2CL, 2.5e6Hz), (GPSL5I, 40e6Hz), (GPSL5Q, 40e6Hz),
    (GalileoE1B, 15e6Hz), (GalileoE1B_BOC11, 5e6Hz),
    (GalileoE1C, 15e6Hz), (GalileoE1C_BOC11, 5e6Hz),
    (GalileoE5aI, 40e6Hz), (GalileoE5aQ, 40e6Hz), (GalileoE5aQP, 20e6Hz),
    (GalileoE5bI, 40e6Hz), (GalileoE5bQ, 40e6Hz), (GalileoE6B, 20e6Hz), (GalileoE6C, 20e6Hz),
    (BeiDouB1I, 10e6Hz), (BeiDouB3I, 40e6Hz), (BeiDouB2bI, 40e6Hz),
    (BeiDouB2aI, 40e6Hz), (BeiDouB2aQ, 40e6Hz), (BeiDouB1C_D, 5e6Hz), (BeiDouB1C_P, 5e6Hz),
)
for (S, fs) in _MS_CASES
    name = String(nameof(S))
    signal = S()
    N = round(Int, ustrip_hz(fs) * 1e-3)
    out = _buf(N)
    SUITE["1 ms code generation"]["$name @ $(_mhz(fs))"] = @benchmarkable gen_code!(
        $out, $signal, $1, $fs,
    ) evals = 10 samples = 1000
end

# ─────────────────────────────────────────────────────────────────────────────
# oversampling sweep — gen_code! across oversampling ratios (fs/fc = samples per chip) at a small
# (4k) and steady-state (64k) buffer. Guards the LUT's ratio-dependent path selection: the windowed
# permute (flat in the ratio) at low oversampling, switching to the exact boundary fill once the
# baked table is oversampled ≳10× (AVX-512) / ≳4× (AVX2). One representative signal (GPS L1 C/A,
# plain BPSK): its oversampling ratio equals the table-oversampling, so 2× exercises the permute
# path and 8×–32× the boundary fill around the switch. 17×/24× are non-power-of-two ratios (mixed
# m/m+1 runs, the magic-reciprocal path); 8×/32× are power-of-two (the shift path) on the SW = 16
# and SW = 64 store widths. 64×/128× cover runs longer than the widest store (the EXTRAS strided
# interior stores) — the region where the old run-fill fell off its `Val` ladder into the slow
# generic kernel.
# ─────────────────────────────────────────────────────────────────────────────
let signal = GPSL1CA(), fc = 1023e3Hz
    for oversampling in (2, 8, 17, 24, 32, 64, 128), (slabel, n) in (("4k", 4096), ("64k", 65536))
        fs = oversampling * fc
        out = _buf(n)
        # evals = 10 (average out per-call scheduling jitter — the 4k fills run in a few
        # hundred ns, where evals = 1 made the reported minimum swing ±25% on shared runners)
        # and samples = 10000, matching the "2000 sample" group, for a stable minimum.
        SUITE["oversampling sweep"]["$(lpad(oversampling, 2, '0'))x"][slabel] =
            @benchmarkable gen_code!($out, $signal, $1, $fs, $fc, $0.0, $0) evals = 10 samples = 10000
    end
end
