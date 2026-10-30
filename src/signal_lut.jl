# ─────────────────────────────────────────────────────────────────────────────
# `SignalLUT`: per-SIGNAL embedded LUT holding ALL PRNs' baked code tables as one matrix.
#
# Included before the signal structs that embed it as a `lut::SignalLUT` field
# (always populated), so the struct field types resolve; built and consumed by the LUT adapter in
# `code_lut.jl` (`build_signal_lut`, `gen_code!`). The metadata (subchip_factor P,
# table_length, period_subchips) is identical across the PRNs of a signal — only the baked
# Int8 column differs — so we store one `padded` matrix whose column `prn` is the per-PRN
# baked padded table, plus the shared metadata once. Resampled zero-copy by `gen_code!`
# (a transient view-backed `CodeLUT.CodeTable` over a matrix column). See `build_signal_lut`.
# ─────────────────────────────────────────────────────────────────────────────
struct SignalLUT
    padded::Matrix{Int8}      # (table_length + WINDOW_PAD) × num_prns; column prn = that PRN's padded baked table
    subchip_factor::Int       # P
    # Residual NON-baked secondary applied per primary period, as a (Ls × n) matrix. `Int8` with
    # a single row of 1 when none/baked. Stored per-PRN (column prn) so a long PER-PRN overlay
    # (GPS L1C-P's 1800-chip code) is applied at runtime WITHOUT baking it (which would blow the
    # matrix up ~1800×). A SHARED secondary (GPS L5I's NH10) stores one column reused for all PRNs.
    secondary::Matrix{Int8}
    table_length::Int         # L·P (length of one PRN's baked table, i.e. column length minus WINDOW_PAD)
    period_subchips::Int      # sub-chips per primary period (Lp·P), for secondary application
    # Per-sample RMS amplitude of the baked code, `sqrt(mean(code^2))` over one primary
    # table. Exactly `1` for ±1 codes (BPSK/BOC/TMBOC); for the multi-level CBOC Int8
    # approximation it is `sqrt(a1^2 + a2^2)` (≈ 19.92 for E1B's baked (19, 6)) and tracks
    # any explicit `cboc_amplitudes` override. Modulation-level and PRN-independent (PRNs
    # differ only in the primary ±1 sign), so one scalar suffices. Read via
    # `get_code_amplitude`; consumers divide correlation results by it to recover a
    # modulation-independent (unit-power) scale.
    code_amplitude::Float64
end

# Residual secondary column for PRN `prn` (the shared-secondary matrix has one column reused
# for every PRN; the per-PRN matrix has one column each).
@inline _signal_lut_secondary(lut::SignalLUT, prn::Int) =
    @view lut.secondary[:, size(lut.secondary, 2) == 1 ? 1 : prn]
