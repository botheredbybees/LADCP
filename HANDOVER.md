# Handover: P16N Cast 003 RMSE Closure

**Date**: 2026-07-11
**Status**: **VALIDATION TARGETS MET — u RMSE 0.0450, v RMSE 0.0333 vs
the archived LDEO 003.nc, both under the 0.05 tolerance, both hard test
assertions.** Stage A (ingestion→editing→weights) and super-ensemble
formation both match Octave LDEO_IX to machine precision on P16N 003;
the solver was exonerated back in P3. See REPORT.md P6 for the final
formation-parity findings and the short list of optional follow-ups
(lanarrow port, Single_Ping_Err parsing, multi-cast validation).

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

Fourth milestone: **3-beam solutions** (commit `472571e`,
`reconstruct_3beam()` + `beam2earth(allow_3beam=True)`): single-missing-
beam cells reconstructed via zero error velocity, applied to DL/UL/BT.

**2026-07-11 session — Stage A closed (see REPORT.md P5):**
- **UL→DL pairing** (`f192477`): `best_ul_shift()` ports loadrdi.m's
  bestlag merge refinement as a SEQUENCE shift (`ul_idx[k−1]`, not
  `ul_idx[k]−1` — staggered pinging makes the difference in ~20% of
  columns; heading fingerprint: 100.0% exact). Stage A velocities
  0.085 → **3.7e-6** rms.
- **Weight construction** (`2e3549c`): `build_ldeo_weights()` ports
  loadrdi.m 408-533 (median-over-beams corr, medianan(maxnan) norm,
  tilt/tiltd masking, echo penalty, non-pinging removal incl. the
  dru-twice bug). Weight rms 0.151 → **0.0020**.

## Current state (2026-07-11, end of session)

- **Archive RMSE: u 0.0450 / v 0.0333 — both validation targets MET**
  (r(u) +0.77; worst stratum 1000-2000 m u now 0.0396).
- Stage A: all fields NEAR vs Octave step09 (velocities ~4e-6, weight
  0.002).
- Formation: machine-precision match vs Octave step10 given identical
  input (formation_only.py: n_se 828=828, ru/rv 5e-8, weight 2e-17).
- Full suite: **237 passed, 8 skipped, 0 xfailed** — the u and v RMSE
  checks are hard assertions.

## Current numbers (`scripts/diag_rmse_strata.py`, NEW config, 2026-07-10,
after depth fix + editing port + sound-speed correction)

| stratum | n | u RMSE | v RMSE | r(u) |
|---|---|---|---|---|
| TOTAL | 520 | 0.0584 | 0.0499 | +0.645 |
| 0–1000 m | 120 | 0.0353 | 0.0313 | +0.882 |
| 1000–2000 m | 121 | 0.0781 | 0.0623 | +0.887 |
| 2000–3000 m | 122 | 0.0497 | 0.0420 | +0.612 |
| 3000–4500 m | 157 | 0.0608 | 0.0560 | +0.379 |

**Trajectory this session (u TOTAL): 0.0678 → 0.0848 (depth fix) →
0.0787 (editing port) → 0.0661 (sound speed) → 0.0584 (3-beam); v
0.0573 → 0.0499 (target met).** The intermediate rise is expected, not a
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

DONE — see REPORT.md P6. Optional follow-ups (not gating validation):

1. Port lanarrow (LDEO step 11) for machine-precision Stage C/D parity
   vs step12.
2. Parse Single_Ping_Err from the PD0 fixed leader so superens_std_min
   is derived, not a hardcoded 0.083833 in the P16N callers.
3. Validate against more casts (S4P processed_uv, 55+ casts) before
   declaring production-ready; then the NetCDF writer/CLI gaps from
   CLAUDE.md are the remaining feature work.

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
| `src/ladcp/transforms/beam2earth.py` | `reconstruct_3beam()`, `beam2earth(allow_3beam=...)` |
| `tests/test_transforms.py` | 3-beam reconstruction unit tests |
| `src/ladcp/ingestion/rdi.py` | `best_ul_shift()` (loadrdi UL merge bestlag, sequence semantics) |
| `src/ladcp/qa/editing.py` | + `build_ldeo_weights()` (loadrdi weight construction) |
| `octave_harness/diag_stage_a_residual.py` | New: per-row/per-phase/weight localization + UL shift scan |
| `octave_harness/formation_only.py` | New: dump-driven formation harness (step09 d -> Python formation vs step10 di) |
| `src/ladcp/solution/inverse.py` | prepinv parity: izm-row window trigger, exact expansion, half-window medianan, _stdnan, tilt weight, post-loop chain |
| `HANDOVER.md` | This file — rewritten |
