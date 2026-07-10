# Handover: P16N Cast 003 RMSE Closure

**Date**: 2026-07-10
**Status**: izm depth-registration offset root-caused and FIXED (see below).
RMSE target (< 0.05 m/s) not yet met; the remaining gap is now cleanly
attributed to Stage A editing/masking differences (the last unexplained
divergence stage).

## Session history (one line each; details in `octave_harness/REPORT.md`)

- **2026-07-05**: transform-convention fix (`f3569c4`, UL beam matrix);
  rotup2down tried and rejected; Octave differential harness built (M1-M4).
- **2026-07-06**: solver exonerated (P3: Python `compute_inverse()` vs Octave
  `getinv.m` agree to ~0.01 m/s on identical input); first-DIVERGES = Stage A;
  depth-varying izm offset quantified (P2).
- **2026-07-10**: izm offset root-caused and fixed (P4). Two mechanisms:
  (1) no caller passed `lat_deg`, so `assign_bin_depths()` used the crude
  `z = p*1.00445` fallback instead of Saunders/`p2z` — ~90 m too deep at the
  cast bottom; (2) Octave shifts its ADCP time labels by its besttlag lagdt
  (−23.24 s), which mispaired the harness comparison by exactly 16 pings
  (staggered ping-cycle aliasing — this had fooled P1's time-mismatch check).
  The physical CTD-ADCP clock offset in Python's data is only −0.5 s,
  measured by the new `estimate_ctd_adcp_lag()`.

## Current state after the P4 fix

- izm vs Octave step09: rms 47.9 m → **1.55 m** (mean +0.02 m, no depth
  correlation). Residual max 7.2 m ≈ the sound-speed bin-length scaling
  Python doesn't yet apply.
- Stage A ru rms vs Octave: 0.215 → **0.117** m/s (the 16-ping mispairing
  was a large chunk of the "unresolved" P1 residual).
- Full suite: **203 passed, 8 skipped, 2 xfailed** (the two RMSE checks).

## Current numbers (`scripts/diag_rmse_strata.py`, NEW config, 2026-07-10)

| stratum | n | u RMSE | v RMSE | r(u) |
|---|---|---|---|---|
| TOTAL | 520 | 0.0848 | 0.0595 | +0.343 |
| 0–1000 m | 120 | 0.0138 | 0.0179 | +0.973 |
| 1000–2000 m | 121 | 0.1540 | 0.0447 | +0.750 |
| 2000–3000 m | 122 | 0.0644 | 0.0916 | +0.343 |
| 3000–4500 m | 157 | 0.0464 | 0.0586 | +0.416 |

**Note the u TOTAL moved 0.0678 → 0.0848 despite the fix being verifiably
correct** (validated against p2z's check value and against Octave's own z to
1.55 m rms). The old, better-looking total was partly error cancellation:
wrong depth registration (+90 m at depth) happened to partially compensate
the Stage A velocity-structure difference in 1000–2000 m. Deep water
(3000–4500 m) — where velocity structure is weak and depth registration
dominates — improved as predicted (0.0720 → 0.0464). Do not "fix" the RMSE
by reverting the depth conversion.

**Note on running tests on this machine**: `uv run pytest` fails with a
broken entry-point shim. Use `TEST_DATA_DIR=test_data uv run python -m
pytest` instead.

## Remaining work to close the RMSE gap (priority order)

The honest first divergence is now Stage A velocities (rms ~0.12 m/s on
both-finite cells, max|diff| 14 m/s from Python-side unmasked garbage).
See REPORT.md P4 handoff for detail:

1. **Port the missing editing steps**: `loadrdi.m::outlier()` (ingestion-
   time outlier rejection; the Octave log shows it discarding ~6952 DL /
   4338 UL bins) and `edit_data.m`'s bin-1 masking ("zero blanking
   distance"). These remove exactly the garbage cells dominating Stage A.
2. **Sound-speed corrections**: velocity scaling `ss/sv`
   (`getdpthi.m:182-207`) and bin-length scaling for izm
   (`getdpthi.m:428-439`).
3. Re-measure `diff_stages.py` and `diag_rmse_strata.py` after 1–2.

## rotup2down: implemented, tested, does not help — not committed

(Unchanged from 2026-07-05; see PROGRAMMERS_NOTES.md "Validation: P16N Cast
003". NOTE: the git stash holding the implementation is **gone** — the
stash list is empty as of 2026-07-10. Reconstruct from
PROGRAMMERS_NOTES.md's line-by-line convention notes against `prepinv.m`
lines 294–418 if ever needed.)

## Files changed this session (2026-07-10)

| File | Change |
|---|---|
| `src/ladcp/ingestion/ctd.py` | `pressure_to_depth()` (Saunders p2z port), `estimate_ctd_adcp_lag()` (besttlag equivalent), `assign_bin_depths(time_offset_days=...)` |
| `tests/test_ctd_loader.py` | Unit tests: p2z check value, lat branch, time offset, synthetic lag recovery |
| `scripts/diag_rmse_strata.py` | run_pipeline passes lat_deg + measured lag; time labels shifted like loadctd.m:443 |
| `tests/integration/test_inverse_p16n_cast003.py` | Same wiring as run_pipeline |
| `octave_harness/diag_izm_root_cause.py` | New: variant test A/B/C/D/E reproducing and closing the P2 offset |
| `octave_harness/diff_stages.py` | Stage A undoes Octave's −23.24 s label shift via heading content-match |
| `octave_harness/REPORT.md` | P4 section + RESOLVED marker on the P2-era priority |
| `HANDOVER.md` | This file — rewritten |
