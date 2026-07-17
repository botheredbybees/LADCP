# Handover: P16N Cast 003 RMSE Closure

## NEXT: remaining exploded I7N casts, cast 018, A16N deep-cast mechanism

Three independent open threads, roughly in order of tractability:

1. **7 untested "exploded" I7N casts** (042 and 060 confirmed fixed by
   commit `5635685`; 062, 086, 099, 102, 103, 118, 119 all showed the
   identical rank-deficiency signature during classification but weren't
   individually re-run post-fix — should just work, worth a quick
   `validate_multicast.py` re-run to confirm and update
   `docs/validation/BULK_VALIDATION_REPORT.md`.
2. **I7N cast 018** — the one "exploded" cast NOT fixed by `5635685`
   (confirmed unaffected: still u RMSE 3.24, well-conditioned, scipy and
   numpy agree to full precision). This is a *different* bug — a real
   solve of bad/extreme input data, not a numerical artifact. Needs its
   own root-cause investigation (likely a data-editing gap: one bad
   SADCP or bottom-track point, or a mis-weighted super-ensemble slipping
   through). Start with `scripts/diag_rmse_strata.py` or a per-bin dump
   on cast 018 to localize which bin/constraint produces the 38.7 m/s
   value at bin index 864.
3. **A16N deep-cast (>4 km) divergence** — still open after a full
   session of diagnostics (see `test_data/2013_A16N/DOWNLOAD_NOTES.md`
   "ps.shear lead — REFUTED"). `ps.shear`/`smallfac`, solve method,
   observation-count starvation, and heave/winch contamination are all
   ruled out with direct evidence. The one unexplained concrete lead:
   merging down-cast and up-cast observations helps shallow casts a lot
   but appears to hurt deep casts relative to a (more crudely
   constrained) down-only/up-only solve — not yet isolated as
   phase-timing-related vs. instrument-source-related (DL vs UL).

Also worth doing once the above settle: re-run the full bulk validation
(`docs/validation/BULK_VALIDATION_BRIEF.md`'s harness, not the brief
itself) to get fresh I7N pass-rate numbers now that 042/060 (and likely
7 more) are fixed — the "10 exploded, 53/124 pass both" headline numbers
in `docs/validation/BULK_VALIDATION_REPORT.md` predate this session's fix.

## 2026-07-17 (fourth session): exploded-cast root cause found + fixed; A16N ps.shear lead refuted

**I7N exploded-cast bug found and fixed (commit `5635685`).** Root cause:
`_solve_lsq()` called `scipy.linalg.lstsq(A, d, check_finite=False)` with
no rank cutoff. scipy's default isn't scaled by matrix dimensions, so on
these large super-ensemble systems it failed to truncate the ~1e-13
singular value produced by a genuinely unconstrained depth bin (zero
observations, untouched by any constraint row — verified SADCP/BT/
barotropic rows only ever touch the CTD-velocity block or degenerate to
zero weight for a zero-observation column). Treated as full rank, that
free column resolved to ~1e10-1e11 m/s instead of the correct
minimum-norm ~0. Switched to `numpy.linalg.lstsq(A, d, rcond=None)`
(dimension-scaled cutoff), verified against the real dumped cast-060
matrix before changing anything. Fixes **9 of the 10** "exploded" I7N
casts (all share the identical signature: rank-deficient by exactly 1,
one all-zero column) — cast 060 now u 0.030/v 0.017 (passes both
targets), cast 042 now u 0.060/v 0.056 (sane). Cast 018 is the exception,
confirmed unaffected (different bug, see NEXT above). Full suite: 256
passed / 8 skipped, no regressions; I7N 003/010 bit-identical to pre-fix;
P16N cast 003 still meets both targets (u 0.0415, v 0.0447 — v shifted
~34% from 0.0333, a real side effect of the correct cutoff truncating
more near-zero directions on that cast too, disclosed and accepted, not
a bug). TDD: failing regression test written first
(`test_solve_lsq_zero_column_does_not_explode`), confirmed red, then
fixed, confirmed green. Built on a short-lived branch
(`fix/lstsq-rank-deficient-explosion`) and fast-forward merged.

**A16N `ps.shear` lead refuted; three more hypotheses ruled out.** Full
writeup in `test_data/2013_A16N/DOWNLOAD_NOTES.md`. Summary: the
documented lead ("port getinv.m's ps.shear rows into compute_inverse")
does not hold up — `getshear2.m` runs strictly *after* the inverse solve
in both the current reference and the exact software snapshot LDEO used
for A16N, and never writes `dr.u`/`dr.v`. `ps.shear=1` and
`ps.smallfac=[1,0]` are identical on passing and failing casts (checked
the real saved `.mat` params), so neither can be the differentiator.
Also ruled out with direct checks: solve method (matches LDEO's
minimum-norm choice), mid-column observation-count starvation (Python
retains *more* data per bin than the reference, not less), the down/up
split as a data-quality proxy (confounded — cruder constraint set, not a
fair comparison), and ship/winch heave contamination (three independent
checks — CTD descent-rate residual, ADCP tilt, W-anomaly editing
rejection rate — all come back flat or inverted relative to what the
hypothesis predicts). Root mechanism still open; see NEXT above.

**Bulk validation report finalized.** The A16N section of
`docs/validation/BULK_VALIDATION_REPORT.md` had been left as "re-run in
progress" from a prior session even though the verified re-run had
completed; a stray duplicate `BULK_VALIDATION_REPORT_A16N.md` (agent-era
data, matched the verified re-run to the digit) was folded in and
removed. I7N section (124/124, 53 pass both) was already accurate.

**Docs reorganized** (commit `36d6187`): root markdown that had
accumulated (`HANDOVER.md`, `VALIDATION_PLAN.md`,
`BULK_VALIDATION_BRIEF.md`/`REPORT.md`, `CONTINUATION_PLAN.md`) moved
under `docs/`; `README.md`/`PROGRAMMERS_NOTES.md`/
`OCEANOGRAPHERS_NOTES.md`/`USER_NOTES.md`/`CLAUDE.md` were all frozen at
a pre-2026-07-11 stale status (RMSE 0.07, xfail, gaps that are actually
implemented) and brought current.

## 2026-07-11 (third session): A16N 2013 — third cruise, two new editing ports

User downloaded UH SOEST processed products for A16N 2013 (124 ref
casts); raw located at NCEI accession 0205839 and 5 casts downloaded
(WH150 DL + WH300 UL + UH-format CTD txt — see
`test_data/2013_A16N/DOWNLOAD_NOTES.md` for everything). Commits
`4d2c5f2` (UH CTD loader, GPS-dday time base), `e6d0164`
(error-velocity elim edit + PPI edit, both LDEO-faithful with per-cast
attrs readback).

**A16N results: 003 u 0.0072/v 0.0094 (best cast yet), 010
0.0216/0.0134 — both PASS on an RDI-150BB DL, a new instrument family.
Deep casts (>4 km: 030/060/090) FAIL** with alternating-sign 0.3–0.8
m/s swings — weak-scatter mid-column underdetermination. **Lead: port
getinv.m's ps.shear=1 shear-solution constraint into compute_inverse**
(all A16N archives ran with it; P16N/I7N matched without it because
their data is dense). Full suite 255 passed / 8 skipped; P16N + I7N
targets still hard-met.

## 2026-07-11 (second session): I7N multi-cast validation unblocked

First cross-cruise validation numbers exist. Three loader/pipeline fixes
(`b4be613`, `b1b5adb`, `edd41f3` — see commit bodies for detail):
interval-based CTD time fallback (I7N cnv has no time column), coarse
CTD-ADCP lag pre-alignment (NMEA time base off by −186 s), loadctd.m
in-water ensemble trim (003 pinged on deck ~7 min), and SADCP-constraint
reconstruction from the reference NC's embedded z_sadcp/u_sadcp vars
(dominant residual: u tilt +0.21/−0.11 m/s without it).

**I7N results: 003 u 0.0488 / v 0.0669, 010 u 0.0336 / v 0.0608.** Both
u under the 0.05 target on a cruise the pipeline never saw; v slightly
over, error is mid-depth large-scale wobbles below SADCP coverage.
rot/offsetup2down measured again (I7N): rot neutral, rot+offset mixed
(010 v 0.0535 but 003 v 0.0795) — defaults stay OFF; the lanarrow port
(follow-up 1 below) is the lead candidate for closing v.

New S4P downloads ("Individual Cast Data"/"Ancillary Data" dirs): 75
LDEO .mat outputs + per-cast .txt protocols + QC PDF are new; the .nc
files duplicate processed_uv/, and there is STILL no raw S4P PD0 —
LADCP_raw.tgz was never downloaded (internet outage; 7 dead .crdownload
files). I7N raw casts 020 (truncated CTD) and 030..110 also await
re-download. Full suite: 247 passed, 8 skipped.

---

# Previous handover (2026-07-11, first session)

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
