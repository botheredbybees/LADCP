# Handover: P16N Cast 003 RMSE Closure

**Date**: 2026-07-05
**Status**: Transform-convention bug found and fixed; medianan fix retained; rotup2down
tried and rejected (worsens RMSE). RMSE target (< 0.05 m/s) not yet met — remaining gap
not yet root-caused.

## Current numbers (`scripts/diag_rmse_strata.py`, convention fix applied, no rotup2down)

| stratum | n | u RMSE | v RMSE | r(u) |
|---|---|---|---|---|
| TOTAL | 520 | 0.0678 | 0.0573 | +0.483 |
| 0–1000 m | 120 | 0.0142 | 0.0195 | +0.970 |
| 1000–2000 m | 121 | 0.1034 | 0.0352 | +0.748 |
| 2000–3000 m | 122 | 0.0463 | 0.0735 | +0.576 |
| 3000–4500 m | 157 | 0.0720 | 0.0737 | +0.014 |

Full test suite (`TEST_DATA_DIR=test_data uv run python -m pytest`): 190 passed, 8
skipped, 2 xfailed (the two RMSE checks in
`tests/integration/test_inverse_p16n_cast003.py`).

**Note on running tests on this machine**: `uv run pytest` fails with `uv trampoline
failed to canonicalize script path` (broken entry-point shim). Use
`uv run python -m pytest` instead — that works.

## What was previously wrong, and what actually fixed it

The ~87–90° DL/UL heading disagreement seen throughout this cast was **not** a
compass fault on either instrument. It was a Python `beam2earth` transform-convention
bug: the uplooker needs a different beam→instrument sign convention than the
downlooker (see `loadrdi.m::b2earth`'s `beams_up`-dependent matrix), and our
implementation wasn't applying it. See `VALIDATION_PLAN.md` for the full evidence
trail that led to this reframe (the previously-proposed "hardcode +87° onto UL
heading" fix, `HANDOVER.md`'s old "Option A", was retired — it would have baked a
cast-specific mounting angle into the pipeline).

Fixed in commit `f3569c4`. `scripts/diag_ul_dl_rotation.py` (the E1 diagnostic)
confirms the residual UL→DL rotation after the fix is noise-dominated around 0°
(not a coherent 87–90° bias): circ mean −8°, downcast +5°, upcast −35°, correlation
with heading/cast-phase near zero, model fit rho ~0.1–0.2 (i.e., the "rotation" model
barely explains the data — it's noise, not a systematic offset).

The separately-applied `_ref_medianan()` fix (MATLAB `medianan` vs `np.nanmedian` for
4-bin reference subtraction) is unrelated and still correctly in place.

## rotup2down: implemented, tested, does not help — not committed

A faithful port of `prepinv.m`'s `rotup2down=1` (harmonize per-ensemble DL/UL heading
fluctuation, after removing the cast-mean offset) is written, unit-tested, and was
wired into the integration test and `diag_rmse_strata.py`. Measured effect: TOTAL u
RMSE goes from 0.0678 to 0.0755 (worse), and the previously-excellent 0–1000 m
stratum degrades most (u RMSE 0.0142 → 0.0327). A sign-bug was ruled out by swapping
the DL/UL rotation directions and re-measuring (0.0754 — same degradation either way;
a real physical-misalignment correction would instead swing the other direction when
flipped). Conclusion: the per-ensemble heading residual here (mean 0.7°, std 7.3°) is
most likely per-instrument compass jitter, not a real mechanical flex between the
frames — applying the correction injects that jitter into u/v.

**The implementation is in `git stash` (not the working tree, not committed).** Run
`git stash list` / `git stash show -p` to recover it if a future cast or instrument
pair shows a genuine per-ensemble heading disagreement worth correcting. See
`PROGRAMMERS_NOTES.md`'s "Validation: P16N Cast 003" section for full detail and the
line-by-line convention verification against `prepinv.m` lines 294–418.

## Remaining work to close the RMSE gap

The 1000–2000 m stratum (u RMSE 0.10) is now the largest single contributor to the
TOTAL RMSE gap, and its cause is **not yet identified** — it survived both the
medianan fix and the transform-convention fix, and rotup2down does not close it.

Recommended next step: `VALIDATION_PLAN.md` Phase 2 — stop debugging via the
single end-to-end RMSE number (it conflates ingestion/transform/solver stages) and
build the solver-only harness (feed LDEO's own S4P `001/002/003.nc` GPS/CTD/SADCP/BT
inputs directly into `prepare_superensembles()`/`compute_inverse()`, compare to
LDEO's own output in the same file) to isolate whether the remaining gap is in the
solver or upstream of it.

## Files changed this session (2026-07-05)

| File | Change |
|---|---|
| `PROGRAMMERS_NOTES.md` | Updated validation section with current RMSE numbers, corrected the "compass offset" root-cause framing, documented rotup2down attempt and rejection |
| `HANDOVER.md` | This file — full rewrite to reflect current state |
