# Continuation Plan — Octave differential harness: run LDEO_IX itself, diff stage-by-stage

**For:** a Sonnet (or cheaper) session. **Date:** 2026-07-05.

**Why this instead of more correction-porting:** a month of validation work has
followed the same loop — port a step, measure end-to-end RMSE, get a new mystery.
Every bug actually found (`medianan`'s `round(n/2)` selection, the up/down
beam-matrix convention, `np.round` half-to-even vs MATLAB half-away-from-zero) is a
semantic hairline crack that line-by-line reading misses and numeric diffing catches
instantly. And every faithful port of a step LDEO definitely runs (`rotup2down`,
`offsetup2down`) made agreement *worse* — which mathematically means the inputs
reaching those steps in our pipeline differ from the inputs that reached them in
LDEO's. The bug is upstream, and the end-to-end RMSE scalar cannot localize it.

**The fix: stop treating LDEO as a black box.** The complete LDEO_IX M-code is in
`docs/legacy/` (flat tree: `process_cast.m`, `loadrdi.m`, `prepinv.m`, `getinv.m`,
`default.m`, plus `LDEO_IX_Software.tar`). The raw inputs for P16N cast 003 (DL/UL
PD0 + CTD cnv) are in `test_data/2015_P16N/`. Run their pipeline under GNU Octave,
dump the `d`/`p`/`di`/`dr` structs after every processing step, and diff each stage
against our Python pipeline. The first stage where the arrays diverge materially IS
the bug. This also regenerates the intermediate bin×ensemble data whose absence
killed the original "solver-only harness" plan (VALIDATION_PLAN.md Phase 1.5 §3).

This supersedes the previous CONTINUATION_PLAN (suspect-by-suspect investigation of
the rotup2down/offsetup2down underperformance) — the harness answers those questions
as a by-product.

## A gift: LDEO recorded their own configuration

`test_data/2015_P16N/003.nc` global attributes contain the **complete `p` struct**
from LDEO's actual run of this cast. Do not guess any parameter — read them from the
file (`netCDF4.Dataset(...)` attributes, or the `LOG_Inverse_log` attribute for the
step-by-step log). Key facts already extracted:

- Software: `Version IX_13beta` (check what version the in-tree code is — look for a
  version string in `default.m`/`misc*.m`; note any mismatch in your report, don't
  try to "fix" it).
- `ladcpdo/ladcpup`: `../data/raw/003DL000.000` / `003UL000.000` — we have both.
- `ctd` and `nav` BOTH point to `../data/CTD/2Hz/003.2Hz`: an ASCII table, 2 Hz,
  `ctd_header_lines: 0`, `ctd_fields_per_line: 11`, time in field 1 (Julian,
  `ctd_time_base: 0`), pressure field 2, temperature field 3, salinity field 4,
  lat field 10, lon field 11, bad values `-999`.
  **We do not have this file — you must generate it** from
  `test_data/2015_P16N/003_01.cnv` (check whether the cnv has NMEA lat/lon columns;
  if not, `0221195_lonlat.txt` has the cast position — constant lat/lon per sample is
  acceptable for a first pass, note it as a caveat). Match the field layout exactly.
  The cnv is 24 Hz-ish; decimate/average to 2 Hz. Compare your generated file's time
  span against the recorded `ctd_starttime: 2457124.7505` / `ctd_endtime:
  2457124.8796` to verify the time base is right before running the pipeline.
- `drot: 12.318441` (magnetic deviation — should be auto-computed by `magdev.m` from
  lat/lon/date, verify it comes out ≈12.32; if magdev errors under Octave, hardcode
  this recorded value in the cast params with a comment).
- No SADCP input for the Octave run (LDEO used one, but `INPUT_SADCP_profile_avail`
  handling can be deferred — see milestones; our own pipeline runs with
  `sadcp_003.npz`, so final-solution comparisons should note the SADCP difference).

## Hard rules

- `docs/legacy/` and `test_data/` are READ-ONLY. Copy the M-code into a new
  `octave_harness/ldeo_ix/` directory and do all editing there. Keep the edits
  minimal and record them (`git diff` of the copy vs original, or a PATCHES.md).
- Commit the harness (driver scripts, converter, patches, README), but gitignore the
  bulky outputs: add `octave_harness/dumps/` and `octave_harness/work/` to
  `.gitignore`.
- All `save` calls that Python must read: use `save('-v6', ...)` (scipy.io.loadmat
  cannot read Octave's default text format or v7.3/HDF5).
- Do not modify `src/ladcp/` in this session. This is an instrumentation session;
  fixing whatever the diff reveals is the NEXT session's job (with the diff report
  in hand).
- Test commands: `TEST_DATA_DIR=test_data uv run pytest` (or
  `./.venv/Scripts/python.exe -m pytest` if uv errors with "trampoline failed to
  canonicalize" — known issue, workaround confirmed).
- Timebox any single obstacle to ~3 focused attempts, then take the documented
  fallback or stop and report. The milestones below are ordered so that partial
  completion is still valuable.

## Getting Octave (Windows host, no Octave installed)

Preferred: Docker (installed on this machine). `docker run --rm -v
"<abs path to LADCP>:/work" -w /work docker.io/gnuoctave/octave:9.2.0 octave-cli ...`
(or latest available tag; `octave-cli` = no GUI, no graphics toolkit needed if
plotting is stubbed). Alternative if Docker is painful: `choco install
octave.portable`. Put the choice and version in the report.

## Milestones (each independently valuable — stop at a milestone boundary if stuck, don't half-finish the next)

**M1 — LDEO ingestion dump (highest value per effort).** Get `loadrdi.m` alone
running on `003DL000.000` and `003UL000.000` under Octave; `save -v6` the resulting
`d`, `p` structs. Then write `octave_harness/diff_ingestion.py`: load the dump, load
the same files with our `load_rdi()`, and report field-by-field max/RMS differences
(velocities, headings, pitch/roll, time, per-beam data; mind MATLAB 1-based bin
indices and any unit differences). **This single milestone answers "does our
ingestion match theirs exactly?" — which has never actually been verified numerically.**
Watch for: `fread` with `'char'` behaves differently in Octave (returns char not
double in some versions — cast explicitly if loadrdi errors); if loadrdi is
unrunnable after timeboxing, export our Python-loaded struct to .mat and note that
M1 is skipped, ingestion equivalence unverified — do NOT silently treat it as passed.

**M2 — full process_cast run with per-step dumps.** Build
`octave_harness/set_cast_params.m` from the recorded p-struct attributes (file paths
pointing at a `work/` dir containing the generated `003.2Hz` and copies/links of the
PD0s), stub plotting, and run `process_cast` through step 14. Instrument it: after
each `end_processing_step`, `save('-v6', sprintf('dumps/step%02d.mat',
pcs.cur_step), 'd', 'p')` plus `di`/`dr` when they exist (wrap in try/catch so a
missing var doesn't kill the run). Plotting under octave-cli: first try `setenv
GNUTERM dumb` / `graphics_toolkit gnuplot`; if figures still error, create
`octave_harness/stubs/` with no-op `figure.m, plot.m, subplot.m, hold.m, axis.m,
title.m, xlabel.m, ylabel.m, text.m, legend.m, colorbar.m, pcolor.m, contourf.m,
streamer.m, orient.m, print.m, pause.m` etc. (add to path AHEAD of the real ones;
add stubs lazily as errors name functions — don't pre-write fifty). Do NOT stub
`set`/`get` unless forced (they're used non-graphically too); if forced, make them
pass through for struct inputs.
Sanity gate for M2 completion: the run's printed log should broadly match the
recorded `LOG_Inverse_log` (e.g. step 10 prints "mean heading offset from compasses
= -89.978 deg" or close; magnetic deviation 12.32). Quote the matching lines in the
report.

**M3 — stage-diff report.** `octave_harness/diff_stages.py`: for each dump, compare
against the corresponding stage of our Python pipeline (reuse the pipeline assembly
from `scripts/diag_rmse_strata.py`; factor it into an importable helper if needed —
that's a scripts/ change, allowed). Map conventions explicitly: MATLAB 1-based
`d.izd`/`d.izu` vs our 0-based; `d.izm` sign; our combined-array bin ordering
(UL reversed on top of DL — see PROGRAMMERS_NOTES.md "Combined DL+UL array layout");
degrees vs radians. Output a table: stage | field | max abs diff | RMS diff |
%-finite mismatch, and a one-line verdict per stage: MATCH (diff within float
noise), NEAR (small, explainable), DIVERGES (material). **The first DIVERGES row is
the product of this whole effort.**
Also: the step-12 `di` dump is the true solver input — save it prominently; it
enables the original solver-only harness (feed it to `compute_inverse()` directly
in a later session).

**M4 (only if M1-M3 done cleanly) — final-solution comparison.** Compare the Octave
run's final `dr` u/v profile against both `003.nc` (LDEO's archived result; should
match closely — if it doesn't, our in-tree code version differs from IX_13beta in a
way that matters, which is itself a critical finding) and our Python result. Note
the SADCP asymmetry (their archived run had SADCP; this Octave run may not).

## Report format

Per milestone: done/skipped + the fallback taken. For M1/M3: the diff tables, and
an explicit statement of the FIRST stage/field that materially diverges. Quote the
sanity-gate log lines (M2). List every edit made to the copied M-code. End with:
what the next session should investigate, based on where the first divergence is
(ingestion / transform / editing / superensemble / solver). Numbers only from
commands you ran.
