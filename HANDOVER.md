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

Same session, second milestone: **ported `loadrdi.m::outlier()` and
`edit_data.m` bin-1 masking** (`edit_outliers()` / `edit_mask_bins()`,
commit `25df9de`) — the Python-side unmasked garbage (14 m/s cells) that
dominated Stage A max|diff| is gone.

Third milestone: **sound-speed correction ported** (commit `d70e10f`,
`ladcp.transforms.soundspeed`): `sounds.m`/`press.m` ports plus
`getdpthi.m`'s velocity/bottom-track/izm-bin-offset scaling by ss/sv.

## Current state after the P4 fix + editing port

- izm vs Octave step09: rms 47.9 m → 1.55 m (depth fix) → **0.36 m**
  (sound-speed bin-offset scaling); max 108.8 → 2.3 m. Depth registration
  is closed.
- Stage A ru rms vs Octave: 0.215 → 0.117 (ping-pairing fix) → **0.080**
  (editing port); max|diff| 14.3 → **3.0 m/s**.
- Stage D vs Octave harness: u rms 0.093 → **0.078**, v 0.063 → **0.062**.
- Full suite: **217 passed, 8 skipped, 2 xfailed** (the two RMSE checks).

## Current numbers (`scripts/diag_rmse_strata.py`, NEW config, 2026-07-10,
after depth fix + editing port + sound-speed correction)

| stratum | n | u RMSE | v RMSE | r(u) |
|---|---|---|---|---|
| TOTAL | 520 | 0.0661 | 0.0541 | +0.626 |
| 0–1000 m | 120 | 0.0268 | 0.0349 | +0.881 |
| 1000–2000 m | 121 | 0.1006 | 0.0685 | +0.800 |
| 2000–3000 m | 122 | 0.0407 | 0.0559 | +0.716 |
| 3000–4500 m | 157 | 0.0696 | 0.0522 | +0.331 |

**Trajectory this session (u TOTAL): 0.0678 → 0.0848 (depth fix) →
0.0787 (editing port) → 0.0661 (sound-speed correction), with r(u)
0.48 → 0.63.** The intermediate rise is expected, not a
regression signal: the old 0.0678 was partly error cancellation — the
wrong depth registration (+90 m too deep at the bottom) partially
compensated the Stage A velocity-structure difference in 1000–2000 m.
Deep water (3000–4500 m), where velocity structure is weak and depth
registration dominates, improved as predicted (0.0720 → 0.0439). Do not
"fix" the RMSE by reverting the depth conversion — it is validated against
p2z's documented check value and against Octave's own z to 1.55 m rms.

**Note on running tests on this machine**: `uv run pytest` fails with a
broken entry-point shim. Use `TEST_DATA_DIR=test_data uv run python -m
pytest` instead.

## Remaining work to close the RMSE gap (priority order)

The honest first divergence is Stage A velocities, now down to rms ~0.08
m/s on both-finite cells, max|diff| ~3 m/s. See REPORT.md's updated "P4
handoff" for detail:

1. **3-beam solutions**: Octave's loadrdi computes 14422 DL / 8473 UL
   3-beam solutions where one beam is bad; Python's `beam2earth()` has no
   3-beam path — those cells either differ or inherit a bad beam. Also the
   remaining mask-policy differences (instrument-nearest bins, 7.4%
   mask_disagree).
2. Re-measure after (1); the 1000–2000 m u stratum (0.1006) is still the
   dominant archive-RMSE contributor. The deep stratum (3000–4500 m) u
   worsened 0.0439 → 0.0696 with the sound-speed correction while all
   else improved — if it persists, verify the bottom-track `sc`
   application against `getdpthi.m:188-197`.

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
| `src/ladcp/qa/editing.py` | `edit_outliers()` (loadrdi outlier() port), `edit_mask_bins()` (edit_data.m bin masking) |
| `tests/test_editing.py` | Unit tests for both new editors (spikes, UL/DL independence, bottom track, no mutation) |
| `src/ladcp/transforms/soundspeed.py` | New: `sound_speed()` (sounds.m), `depth_to_pressure()` (press.m), `apply_sound_speed_correction()` (getdpthi.m scaling) |
| `tests/test_soundspeed.py` | Unit tests incl. Octave-measured sounds.m parity value |
| `HANDOVER.md` | This file — rewritten |
