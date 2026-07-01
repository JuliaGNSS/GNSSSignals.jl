# Permute steady-state profiling (AVX-512, Zen 5) — experiment branch

Scratch branch off PR #69 (`feat/code-lut-resampler`) to look for a permute-path speedup.
Baseline: `_generate!` (AVX-512 permute) steady-state ≈ **37 ps/sample** (flat in oversampling),
which is ~12× above raw AVX-512 store bandwidth — so there is theoretical headroom. All numbers
GPS L1 C/A, large N (16k–64k), `julia -t1`, `@belapsed`.

## Where the time goes (diagnostic patches, wrong output / timing only)

| variant | ps/sample | note |
|---|--:|---|
| baseline (4 streams, UInt32 rem, _B=30) | **37.4** | |
| strip the DDA advance entirely | **9.2** | → advance is **~75%** of the cost; lookup+store is only ~9 ps |
| remainder narrowed to UInt16 (_B=15) | **20.4** | halving remainder width ≈ **1.8×** → the UInt32 remainder vector is the bottleneck |
| stub the scalar carry-extract `h=c[1]` (UInt32) | **35.6** | only ~3 ps → extract/base chain is NOT the bottleneck |

Conclusion: the permute is **advance-bound**, and the advance is dominated by the
**per-lane remainder vector** `Vec{64,UInt32}` (4 zmm/stream × 4 streams) — its add/compare/
conditional-subtract throughput, not the vector→scalar extract, not the lookups.

## Why the obvious win (narrow the remainder) is blocked

`UInt16` remainder caps `_B ≤ 15` (carry `rem+frac_step` must stay < 2^16). At `_B=15` the rate
`step_num/2^15` has ~2e-5 relative error, giving ~0.25-chip drift over a 1 ms epoch — enough to
shift a chip boundary. The full test suite fails there:
`L1C-D / L1C-P sample-stream matches external reference (12.276 MHz)`.
`_B` must stay ≥ ~24 to keep byte-exactness against the external reference ⇒ **UInt32 remainder
is required**. There is no precision-preserving way to narrow it.

Other levers checked and rejected:
- More/fewer streams: the advance op *count per sample* is fixed regardless of stream count;
  the chain is throughput-bound (chain latency ≪ stride time), so more ILP doesn't help.
- Sharing the init multiply across streams (`_init_rel4`): ~2 ps only (the per-stream Int8
  narrow is inherent); not the steady-state anyway.
- Recomputing `rel` from a scalar phase each stride: needs `Vec{64,Int64}` (8 zmm) — more
  expensive than the incremental UInt32 remainder, not less.

## Bottom line
No precision-preserving speedup found for the permute steady-state; the UInt32 remainder is
both required (precision) and the throughput bottleneck. lookup+store (~9 ps) and the run-fill
path are already near-optimal. A real win would need a different resampling scheme that avoids a
full-width per-lane remainder while keeping ~2^-24 rate precision — an open research question,
not a safe incremental change.
