# Continuation Plan — Phase 2: port offsetup2down + the iterative first-guess loop

**For:** a Sonnet (or cheaper) session. **Date:** 2026-07-05.

**Situation:** Phase 1 is done. `rotup2down` (heading-fluctuation harmonization) is
implemented, unit-tested, and verified correct in isolation — but tested alone it
worsened RMSE, and was correctly NOT committed (preserved in `git stash@{0}`, "WIP on
main: a56ec48 ..."). A prior orchestrating session (Fable) then re-read the LDEO
processing logs and `process_cast.m`/`prepinv.m` and found *why*: see
`VALIDATION_PLAN.md`'s "Phase 1.5" section (added 2026-07-05) before doing anything
else. Summary of what changed:

1. LDEO applies `rotup2down` **and** `offsetup2down` by default on every cast
   (`docs/legacy/default.m:223` and the line above it) — these are paired, not
   alternatives. `offsetup2down` is a **velocity offset** correction between UL/DL
   (distinct from `rotup2down`'s **heading rotation**), and it needs a first-guess
   ocean-velocity profile as input.
2. LDEO's real pipeline is **iterative**, not single-pass: form super-ensembles with
   `rotup2down` only (`process_cast.m` step 10) → preliminary solve to get a first-guess
   profile, trimming ~1% outliers (step 11, `lanarrow`) → **re-form** super-ensembles
   using that first guess, this time applying `offsetup2down` and an incremental
   `rotup2down` adjustment (step 12) → final solve (step 14).
3. Testing `rotup2down` alone, once, before any solve, is testing an incomplete
   procedure out of its real context — that's the most likely reason it made RMSE
   worse. It is not evidence the heading-rotation port itself is wrong (`_compoff`
   independently checks out: −90.84° vs LDEO's own logged −89.978° for P16N cast 003).
4. The originally-planned "solver-only harness" (feed LDEO's own `.nc`/`.mat` inputs
   directly into `prepare_superensembles()`/`compute_inverse()`) **will not work** —
   confirmed via `scipy.io.loadmat` that `test_data/2018_S4P/003.mat` contains only
   LDEO's final per-depth output profiles and final per-super-ensemble nav/CTD series,
   no bin×ensemble matrix. Do not spend time trying to parse these files for solver
   inputs; that data isn't there.

## Hard rules

- NEVER hardcode a compass angle or a velocity-offset constant. Everything must be
  computed per-cast from the data, as LDEO does.
- Do not touch `docs/legacy/` or `test_data/` (read-only).
- Do not attempt the "solver-only harness" as literally described anywhere upstream —
  it's retired; see point 4 above.
- Commands: `TEST_DATA_DIR=test_data uv run pytest` (full suite),
  `uv run pytest tests/test_inverse.py` (fast unit),
  `uv run python scripts/diag_rmse_strata.py` (stratified RMSE), `uv run ruff check src/ tests/`.
- Every number in your report must come from a command you ran.

## Steps (in order, verify each before the next)

1. **Read context** (30-45 min): `VALIDATION_PLAN.md` in full (especially "Phase 1.5"),
   then `docs/legacy/process_cast.m` lines ~340-405 (the step 10-12 sequence) and
   `docs/legacy/prepinv.m` lines ~34, 162-230 (the `offsetup2down` block itself — note
   it needs `dr`, a first-guess profile struct with fields `z`, `u`, `v` at minimum).
   Recover the stashed `rotup2down` port: `git stash apply stash@{0}` (apply, not pop,
   until this phase is proven out) — read what it already gives you
   (`rotup2down()` in `src/ladcp/solution/inverse.py`, its unit tests in
   `tests/test_inverse.py`).

2. **Design the iteration** (this is the actual new work): figure out how to slot a
   preliminary `compute_inverse()` call into the current pipeline to produce a
   first-guess `u(z)`, `v(z)` profile, matching `dr`'s role in `prepinv.m`. You do not
   need to replicate LDEO's outlier-trimming (`lanarrow`) exactly — a single
   preliminary inverse solve without `offsetup2down` (i.e. exactly what
   `prepare_superensembles()`/`compute_inverse()` already produce today, using only
   `rotup2down`) is a reasonable first-guess source; note this as a documented
   simplification, not a silent shortcut.

3. **Port `offsetup2down`** as a new function alongside `rotup2down` in
   `src/ladcp/solution/inverse.py`: faithful to `prepinv.m:177-208` — interpolate the
   first-guess profile onto each ensemble's bin depths, subtract from raw DL/UL
   velocities, take the median residual per instrument (`medianan` — already ported,
   reuse it), split the UL−DL difference in half with opposite sign onto UL/DL (and
   bottom track, per the MATLAB block), scaled by `p.offsetup2down` (use 1, matching
   `default.m`). Unit-test it the same way `rotup2down` was tested (synthetic
   dual-instrument ensemble, known offset in → known correction out).

4. **Wire the full loop** into the P16N cast 003 integration path (mirror
   `test_inverse_p16n_cast003.py`'s existing structure): first solve (no offset
   correction) → derive `dr` → `rotup2down` (full) + `offsetup2down` using `dr` →
   final solve. Confirm this matches the *intent* of process_cast.m steps 10-12, not
   necessarily every incidental detail (e.g. LDEO's "rotated earlier, use difference"
   optimization is a performance shortcut, not behaviorally required — recomputing the
   full rotation each time should be numerically equivalent; note if it isn't).

5. **Baseline first**: before wiring the loop in, run
   `uv run python scripts/diag_rmse_strata.py` at current HEAD (convention fix only, no
   rotup2down, no offsetup2down) and record stratified u/v RMSE — this is your
   baseline, since the last recorded baseline (in `PROGRAMMERS_NOTES.md`) predates
   today's session.

6. **Measure with the full loop**: full test suite + `diag_rmse_strata.py` again.
   Compare to step 5. Also compute hoff/hrot diagnostics if useful for a sanity check
   against the LDEO log's own numbers (P16N cast 003: mean heading offset from
   compasses ≈ −89.978°, from pitch/roll ≈ −91.15°; S4P cast 003: ≈ −88.85° / −85.30° —
   these are in the `LOG_Inverse_log` global attribute of the respective `.nc` files,
   readable via `netCDF4.Dataset(...).LOG_Inverse_log`).

7. **Decide by numbers**: if the full loop improves (or doesn't harm) RMSE relative to
   step 5's baseline and all non-xfail tests pass → commit as
   `feat: port prepinv.m offsetup2down + iterative first-guess loop`. If it still
   worsens things, do not commit the pipeline wiring — stash, document exact numbers
   in `HANDOVER.md`/`PROGRAMMERS_NOTES.md`, and report; that's still a valid, useful
   result (it would mean the remaining gap is upstream of both corrections — worth
   knowing).

8. **Documentation**: update `PROGRAMMERS_NOTES.md` (current stratified RMSE, what's
   applied: convention fix + rotup2down + offsetup2down [note which were kept]) and
   `HANDOVER.md`. Update `VALIDATION_PLAN.md`'s Phase 1.5/Phase 2 status. If stratified
   RMSE now meets the Phase 2 §4 criteria, un-xfail the RMSE tests; otherwise update
   their reason strings to current numbers.

## Report format

Baseline (step 5) vs final (step 6) stratified RMSE; whether offsetup2down/the
iterative loop helped and by how much; test counts; commits made; explicit
recommendation for what's next (if the gap persists: Phase 3 test-data expansion from
`VALIDATION_PLAN.md`, particularly getting S4P raw PD0 so a second cast can confirm
this isn't P16N-specific).
