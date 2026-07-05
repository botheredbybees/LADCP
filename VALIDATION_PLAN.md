# Plan: RMSE closure and validation strategy

**Date:** 2026-07-05 (Claude Fable 5, after reading git history, HANDOVER.md,
PROGRAMMERS_NOTES.md, the LDEO cast-003 processing log, prepinv.m, and the test-data
inventory). Supersedes the "Remaining Work" section of `HANDOVER.md` (2026-06-27).

**Update 2026-07-05 (later same day):** Phase 1's `rotup2down` implementation is done,
tested, and correctly ported (see "Phase 1 result" below) — but it was tested in
isolation and made RMSE worse. Two follow-up findings from re-reading the LDEO logs and
`prepinv.m`/`process_cast.m` change what Phase 2 should actually do; see "Phase 1.5" and
the corrected Phase 2 §1 below. **Read those before starting Phase 2 work.**

## Phase 1 result (2026-07-05): rotup2down alone does not close the gap — as expected

`rotup2down` (prepinv.m's per-ensemble DL/UL heading-fluctuation harmonization) was
implemented as a faithful port (`_compoff` + `rotup2down` in `src/ladcp/solution/inverse.py`,
unit-tested), verified against P16N cast 003, and found to worsen u RMSE almost everywhere
(0–1000 m RMSE more than doubled, 0.0142→0.0327). A sign-bug was ruled out by swapping the
DL/UL rotation directions (0.0754 vs 0.0755 — nearly identical). **Not committed**; the
implementation is preserved in `git stash` (`stash@{0}`, "WIP on main: a56ec48 ...") rather
than deleted, recoverable via `git stash list` / `git stash show -p stash@{0}`.

This result is *expected* given Phase 1.5 below — do not read it as "rotup2down is broken."
It was tested as a standalone one-shot correction on raw data, which is not how LDEO
actually invokes it (see next section). The `_compoff` offset computation itself is
independently verified correct: run without proper UL-index time-alignment it still gives
hoff = −90.84° for P16N cast 003 vs LDEO's own logged −89.978° — well within expected
tolerance for the alignment shortcut used in the check.

## Phase 1.5 — Why rotup2down alone can't reproduce LDEO output (new, 2026-07-05)

Two things confirmed by reading the LDEO processing logs (embedded as the
`LOG_Inverse_log` global attribute in every reference `.nc`) and the legacy MATLAB source,
that were not visible from `prepinv.m` alone:

**1. LDEO applies TWO up/down harmonization corrections by default, not one.**
`docs/legacy/default.m:223` sets `p.offsetup2down = 1` (alongside `p.rotup2down = 1`,
line ~219) — both defaults, used on every cast unless a cruise's `cruise_params.m`
overrides them (checked: neither `test_data/2018_S4P/set_cast_params.m` nor
`test_data/ancillary/set_cast_params.m` does). `offsetup2down` is a **velocity offset**
correction between UL and DL (`prepinv.m:177-208`), separate from `rotup2down`'s
**heading rotation**. It shifts UL/DL velocities by half the median UL−DL residual
velocity (computed after subtracting a first-guess ocean velocity profile `dr`), split
with opposite sign the same way `rotup2down` splits its heading residual.

**2. offsetup2down requires a first-guess solve — LDEO's pipeline is iterative, not
single-pass.** `process_cast.m` steps 10-12 show the real sequence:
- Step 10 "FORM SUPER ENSEMBLES": `prepinv(d,p)` — applies `rotup2down` only (no `dr` yet).
- Step 11 "REMOVE SUPER-ENSEMBLE OUTLIERS": runs a preliminary solve (`lanarrow`) to get a
  first-guess profile, iteratively trimming ~1% outliers.
- Step 12 "RE-FORM SUPER ENSEMBLES": `prepinv(d,p,dr)` — re-invoked **with** the first-guess
  profile `dr`, this time applying `offsetup2down` (logged as "adjusted for velocity offset
  in up and down looking ADCP") and re-applying `rotup2down` ("rotated earlier, use
  difference" — it does NOT redo the full rotation, only an incremental adjustment).

Confirmed present in the log for **both** reference casts (`test_data/2018_S4P/003.nc`
and `test_data/2015_P16N/003.nc`) — this is standard production behavior, not a
cast-specific quirk.

**Consequence:** the Python port's current `rotup2down` call — a single invocation on raw
data before `prepare_superensembles()`, with no first-guess subtraction, no paired
`offsetup2down`, and no iteration — implements a different, incomplete procedure from what
produced the reference output. That plausibly explains why it made RMSE worse rather than
better: it's not "rotup2down is wrong," it's "half of a two-part, iterative correction,
applied out of context."

**3. The archived `.nc`/`.mat` reference files are LDEO's FINAL output, not an intermediate
checkpoint — this breaks the literal Phase 2 §1 plan below.** Checked
`test_data/2018_S4P/003.mat` (scipy.io.loadmat): its `dr`/`da`/`p`/`ps` structs contain only
final per-depth-bin profiles (`u`, `v`, `u_do`, `u_up`, `v_do`, `v_up` on the 344-level `z`
grid) and final per-super-ensemble nav/CTD series (the 550-length `tim` grid) — the same
content as the `.nc`, just also in `.mat` form. **There is no `d.ru`/`d.rv`/`d.weight`/
`d.izu`/`d.izu` matrix (bin × ensemble) anywhere in the archive** — the actual boundary
type between "prep" and "solve" in our own pipeline. LDEO's log references a MATLAB
`checkpoints/003_1` save file that would have this, but it isn't in `test_data/` and isn't
retrievable without going back to Thurnherr/LDEO or re-running their MATLAB stack.
**"Feed LDEO's own inputs directly into `prepare_superensembles()`/`compute_inverse()`"
(Phase 2 §1, original wording) cannot be done with data currently on disk.**

## The reframe: the "87° compass offset" is probably not a compass problem

Evidence assembled this session:

1. **The LDEO processing log for this exact cast** (`test_data/2015_P16N/003.txt`) shows
   Thurnherr processed 003 with NO compass correction — only magnetic deviation 12.318°.
   If the UL compass were genuinely 87° wrong, the reference output would be garbage too.
2. **The 87° difference is in the RAW heading readings** (DL 115.9°/UL 29.7° downcast),
   before any Python transform — so it's a property of the instruments, and the most
   ordinary explanation is a **physical mounting-azimuth difference**: the two instruments
   are simply bolted to the rosette rotated ~87° apart, and BOTH compasses read correctly.
3. In that world, each instrument's beam2earth using its OWN heading yields Earth-frame
   velocities that already agree — no offset correction needed, which is exactly why the
   LDEO log shows none, and why `prepinv.m`'s `rotup2down` only harmonizes the small
   per-ensemble fluctuation after removing the mean offset (`hoff`).
4. Therefore: if Python's UL Earth velocities come out ~87°-rotated from DL's (the
   observed "UL u measures v_ocean"), **the defect is in our UL transform path** — the
   UL's heading is being lost, misapplied, or applied with mirrored sense — not in the
   data. Prime suspect: the "negate pitch" shortcut (`beam2earth(..., -rdi_ul.pitch, ...)`)
   is not the full inverted-instrument convention (roll/beam-permutation handling), and an
   upside-down compass's rotation composes with opposite sense if the frame handedness
   isn't flipped consistently.

**Consequence: retire HANDOVER Option A (add +87° to UL heading).** It hardcodes a
cast-specific physical mounting angle into the pipeline; the next rosette (including
Nuyina's own) will have a different angle and the fudge becomes a new mystery. Option B
(rotup2down) alone won't close the gap either, per the correct reading of prepinv.m — it
corrects fluctuations, not the constant, because the constant needs no correction when the
transform is right.

## Phase 0 — Preserve completed work (half a day)

- `src/ladcp/solution/inverse.py` (`_ref_medianan` fix) and `PROGRAMMERS_NOTES.md` are
  **uncommitted** — this is validated root-cause work sitting in the working tree. Commit
  first, before anything else touches inverse.py.
- Commit `HANDOVER.md`, the `scripts/diag_*.py` diagnostics (they're the debugging
  toolkit, worth keeping), and clean the strays: `$HOME/`, `calude.bat`,
  `"The progress ledger still shows the SADC.txt"` (a mis-saved prompt), `opencode.json`
  (decide: track or ignore).

## Phase 1 — Discriminate, then fix the UL transform (the core work)

**E1 (the decisive experiment, ~an hour):** During ensembles where UL and DL bins sample
the SAME depth at the same time (mid-cast overlap), compare Python's Earth-frame (u,v) per
instrument, each using its own compass. Fit the rotation angle UL→DL per ensemble; plot it
against package heading and cast phase.
- ~0° → transform fine; the RMSE gap is elsewhere (go to Phase 2 with the solver harness).
- Constant ≈ 87° → UL heading effectively unused/cancelled in our transform.
- Varies with heading (e.g. ∝ 2·hdg) → mirrored-frame composition (the inverted-instrument
  convention bug).

**E2 (reference behaviour):** Same comparison inside LDEO's own products —
`test_data/2018_S4P/001.nc` embeds per-instrument DL/UL profiles. Confirm LDEO's UL and DL
agree in overlap (expected), establishing the target behaviour with zero new code.

**C1 (convention audit):** Line-diff our `beam2earth()` inverted handling against the two
authorities: `docs/legacy/ADCPtools/janus5beam2earth.m` (explicit uplooking treatment —
CLAUDE.md already names it authoritative) and `loadrdi.m`'s UL path (incl. the sysconfig
up/down bit, which our PD0 parser should read and assert against — we currently trust the
caller to know which file is the UL). Also parse the fixed-leader heading-alignment/bias
(EA/EB) words for both instruments and confirm they're zero (our parser currently ignores
them; loadrdi accounts for instrument config).

**Fix + verify:** Implement the correct inverted-instrument convention (not a per-cast
constant), re-run E1 (expect ~0°), then the P16N integration test. Only AFTER the
transform is right, implement `rotup2down` fluctuation smoothing as a faithful
prepinv.m port (it's small once the constant is gone) — note prepinv also offers a
velocity-derived rotation (`hrotvel = angle(uu/ud)`), which doubles as a permanent
built-in E1 diagnostic worth exposing in `qa/`.

## Phase 2 — Abandon single-scalar, single-cast validation

The medianan and heading episodes cost weeks because one end-to-end RMSE number conflates
five pipeline stages. Change the validation architecture:

1. **~~Solver-only harness~~ (retired as originally stated — see Phase 1.5 finding 3):**
   the archived `.nc`/`.mat` files do not contain superensemble-matrix inputs, so they
   cannot drive `prepare_superensembles()`/`compute_inverse()` directly. Revised step 1,
   in priority order:
   a. **Port `offsetup2down` + the iterative first-guess loop** (steps 10-12 of
      `process_cast.m`: form super-ensembles with `rotup2down` only → preliminary solve →
      re-form with `offsetup2down` + incremental `rotup2down` using the first-guess `dr`).
      This is not optional/nice-to-have — `default.m` turns both corrections on for every
      cast, so no real cast will validate without it. The stashed `rotup2down` port
      (`stash@{0}`) is reusable as the heading half of this; recover it first
      (`git stash apply stash@{0}` — apply, not pop, until this is proven out).
   b. **Re-test end-to-end RMSE only after (a) is wired in** — this is the fair test of
      whether the transform layer (Phase 1) is actually clean, since Phase 1's isolated
      rotup2down test was confounded by the missing offsetup2down/iteration structure.
   c. If a true stage-isolated solver check is still wanted after that, it requires either
      obtaining LDEO's own intermediate checkpoint (contact Thurnherr/LDEO — not in this
      repo) or accepting that our own `prepare_superensembles()` output (post our full
      ingestion+transform+edit+prepinv-equivalent stages) is the only per-ensemble matrix
      available, which reintroduces upstream stages rather than isolating the solver.
2. **Transform-only check:** E1/E2 above, promoted to a permanent integration test.
3. **Ingestion-only checks:** already exist (PD0 header tests).
4. **Acceptance criteria:** replace flat `RMSE < 0.05` with (a) RMSE within LDEO's own
   per-cast error estimate (the reference files carry `uerr`-style uncertainty), and
   (b) depth-stratified thresholds (0–1000 m is already at 0.025; the deep gap is the
   signal to isolate, not average away). Keep the xfail tests until criteria met, but
   re-point them at the stratified thresholds.

## Phase 3 — More test data (priority order, mostly already local)

1. **Nuyina cast 202324050_004** (`test_data/Nuyina/` — DL+UL PD0 + SBE hex + XMLCON,
   already on disk). The SBE hex decoder was completed in the last commits, so this cast
   is now processable end-to-end. It has NO reference output — its role is *robustness*
   (does the pipeline run on our own instrument config?) and it's the actual mission
   target. Make it a smoke test.
2. **S4P raw PD0** — the known-but-not-downloaded archive (NCEI GOSHIP-LADCP collection;
   accession pattern as per `test_data/sources.md`, which already maps the sources).
   Even 3–5 stations give raw+processed pairs from the SAME processor as our references —
   turns validation statistical instead of forensic, and directly tests whether the P16N
   87° mounting angle is cast-specific (it will differ per cruise — the proof the fudge
   was wrong).
3. **`processed_noedit` variants** (CLIVAR archives carry both, per sources.md) — lets us
   separate editing-stage differences from solver differences when a cast disagrees.
4. Correction to CLAUDE.md while there: `test_data/cruise_data/` holds processed
   tarballs only — the raw I7N PD0 claimed in CLAUDE.md is not present.

## Phase 4 — Only then: features

CLI wiring, NetCDF writer polish, remaining `ladcp2cdf` parity — all deferred until
Phases 1–2 are green. Rationale: every feature built on an unvalidated transform layer
inherits the doubt.

## Assumptions explicitly retired / retained

| Assumption | Verdict |
|---|---|
| "UL compass is 87° off; correct it" | **Retired** — evidence points to mounting azimuth + a Python transform-convention bug; fix the convention, never hardcode the angle |
| "rotup2down is the missing step for the constant offset" | **Retired as stated** — it handles fluctuations only; still worth porting after the transform fix |
| "rotup2down alone should close (part of) the RMSE gap" | **Retired 2026-07-05** — LDEO always pairs it with `offsetup2down` inside an iterative first-guess loop (`default.m` sets both to 1); tested alone, out of that context, it made RMSE worse. Not a rotup2down bug — port `offsetup2down` + the loop before re-testing |
| "S4P/P16N `.nc`/`.mat` files can drive a solver-only harness directly" | **Retired 2026-07-05** — confirmed via `scipy.io.loadmat` that they hold only final per-depth profiles and final per-super-ensemble nav/CTD series, no bin×ensemble matrix; LDEO's own intermediate checkpoint isn't archived here |
| "negate UL pitch is the complete inverted-instrument treatment" | **Under test (E1/C1)** — prime suspect |
| "validate end-to-end RMSE<0.05 on one cast" | **Retired** — stage-wise harness + stratified, uncertainty-aware criteria |
| "P16N cast 003 is the only raw+reference pair available" | **Retired** — S4P raw is retrievable; Nuyina cast already local |
| Validation-first principle (reproduce LDEO before new features) | **Retained** — this plan is that principle, applied with better instrumentation |
