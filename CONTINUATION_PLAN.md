# Continuation Plan — investigate why the paired rotup2down+offsetup2down loop still underperforms

**For:** a Sonnet (or cheaper) session. **Date:** 2026-07-05.

**Situation:** Phase 2 §1a of `VALIDATION_PLAN.md` is done in this same session (by
the orchestrating Fable session directly, not handed off) — `offsetup2down()` is
implemented, TDD'd (4 unit tests), and wired into `scripts/diag_rmse_strata.py`
alongside the previously-stashed `rotup2down()` (both are now committed as tested
library functions, not stashed). Read `VALIDATION_PLAN.md`'s "Phase 2 result"
section and `PROGRAMMERS_NOTES.md`'s matching section in full before doing
anything — they have the complete numbers and reasoning. Short version:

| config | TOTAL u RMSE | 0–1000m u RMSE |
|---|---|---|
| convention fix only (baseline, no correction) | 0.0678 | 0.0142 |
| + rotup2down only | 0.0755 | 0.0327 |
| + rotup2down + offsetup2down (iterative) | 0.0865 | 0.0207 |

Pairing the two corrections (which is what LDEO actually does by default on every
cast — confirmed from `default.m` and the processing logs) claws back about half of
`rotup2down`-alone's damage, but the full pairing is still worse everywhere than
applying **neither** correction. This means: the pairing hypothesis was right (it
explains why testing rotup2down alone regressed), but something else is still
wrong — porting both corrections faithfully didn't close the gap.

Neither correction is wired into the production pipeline
(`tests/integration/test_inverse_p16n_cast003.py`'s fixture uses convention-fix-only;
this is deliberate, not an oversight — do not wire either in without new evidence).

## Hard rules

- Do not wire `rotup2down`/`offsetup2down` into the production integration fixture
  unless you have RMSE numbers proving the combination beats the current baseline
  (0.0678 TOTAL / 0.0142 at 0-1000m) — not just beats each other.
- Do not touch `docs/legacy/` or `test_data/` (read-only).
- Commands: `TEST_DATA_DIR=test_data uv run pytest` (full suite, or
  `./.venv/Scripts/python.exe -m pytest` directly if `uv run` errors with "uv
  trampoline failed to canonicalize script path" — a real issue hit this session,
  workaround confirmed working),
  `uv run python scripts/diag_rmse_strata.py` (stratified RMSE, 4 configs),
  `uv run ruff check src/ tests/ scripts/`.
- Every number in your report must come from a command you ran.

## Investigate, in priority order (from VALIDATION_PLAN.md's "remaining suspects")

1. **The first-guess simplification.** The current loop uses one preliminary solve
   (rotup2down only, no outlier trimming) as the first guess fed to
   `offsetup2down`. LDEO's actual step 11 ("REMOVE SUPER-ENSEMBLE OUTLIERS") runs
   `lanarrow` — an iterative ~1%-per-round outlier trim — before its first guess is
   considered final. Try: (a) iterate the rotup2down→solve→offsetup2down loop 2-3
   times (feed each iteration's output back as the next first guess) and see if
   RMSE converges toward the baseline or diverges further; (b) approximate outlier
   trimming with a simple percentile clip on `SuperEnsemble.weight` or residuals
   before the first solve, and re-measure. Report numbers either way — a
   convergence trend (even if not fully closing the gap) is meaningfully different
   from no change.

2. **Re-run the E1 diagnostic with both corrections applied.** `scripts/diag_ul_dl_rotation.py`
   (from commit `a56ec48`) measures the residual UL→DL rotation angle; Phase 1
   found it noise-dominated (~0°) with the convention fix alone. Run it again on
   ensembles that have gone through `rotup2down` + `offsetup2down` — if the
   residual rotation becomes coherent (not noise) after these corrections, that
   points to the corrections interacting badly with a residual transform defect
   rather than the corrections being independently fine.

3. **If (1) and (2) don't explain it**, the next real lever is Phase 3 §2 of
   `VALIDATION_PLAN.md`: get a second raw+reference cast (S4P raw PD0, not yet
   downloaded — `test_data/sources.md` has the accession pattern) to check whether
   "pairing helps directionally but doesn't close the gap" replicates on a
   different cast, or is specific to P16N cast 003's particular UL/DL geometry.

## Report format

For whichever of (1)/(2) you run: exact RMSE numbers (stratified, matching the
table format above) before/after, and E1 residual-rotation numbers if you ran that
diagnostic. State plainly whether the investigation moved the needle or not — a
clean negative result is as valuable as a positive one here (this validation
effort has already spent real time chasing false leads; don't oversell an
inconclusive result). Recommend the next concrete step.
