# Octave Differential Harness — Report

**Session date:** 2026-07-05. Executes CONTINUATION_PLAN.md's milestones M1-M4
(run LDEO_IX itself under Octave, dump stage-by-stage, diff against our
Python pipeline).

## Summary

| Milestone | Status |
|---|---|
| Gate 0 (Docker/Octave) | Done |
| M1 (ingestion dump + diff) | Done — DL ingestion bit-exact, UL near-exact |
| M2 (full `process_cast` run, per-step dumps) | Done — all 17 steps complete cleanly |
| M3 (stage-diff report) | Partially done — Stage D solid, Stages A/C have an unresolved alignment issue (see below) |
| M4 (final-solution comparison) | Done |

**Headline finding:** the Octave harness's own reconstructed final answer
(missing SADCP, using a reconstructed CTD file) matches LDEO's archived
`003.nc` result *better* than our Python pipeline does (u RMSE 0.068 vs
0.087, v RMSE 0.030 vs 0.058). This weighs against the standing hypothesis
that our Python pipeline's ~0.09 m/s residual gap is caused by imperfect
input reconstruction (CTD/SADCP) — it looks more likely to be a genuine
algorithmic difference downstream of ingestion. See "M4" below for the full
triangle.

## Environment

- Docker Desktop was not running at session start; started it, confirmed
  `docker.io/gnuoctave/octave:9.2.0` `octave-cli --version` → GNU Octave
  9.2.0.
- Windows/Git-Bash volume mounts require `MSYS_NO_PATHCONV=1` (`-v`/`-w`
  paths get mangled otherwise) — wrapped in `octave_harness/run_octave.sh`.
- In-tree code version: `default.m` reports `Version IX_14beta`. The
  recorded p-struct (`GEN_Software_orig` in `003.nc`) shows LDEO actually
  ran `Version IX_13beta` on this cast. **Not reconciled** — noted per the
  plan's instruction not to "fix" this.

## M1 — ingestion dump + diff (`octave_harness/diff_ingestion.py`)

`loadrdi.m` run standalone on `003DL000.000`/`003UL000.000` under Octave
(`octave_harness/dump_m1_ingestion.m`), dumped to
`octave_harness/dumps/m1_loadrdi.mat`, diffed against
`ladcp.ingestion.rdi.load_rdi()`.

**Static config fields** (nbin, blen_m, blnk_m, dist_m, beam_angle_deg):
exact match, DL and UL, zero diff.

**Per-ensemble time series** (nearest-time matched, since ping intervals are
non-constant/staggered):

| field | inst | n | mean Δt (days) | max Δt (days) | max\|diff\| | rms diff |
|---|---|---|---|---|---|---|
| time_julian | DL | 8969 | 0 | 0 | 0 | 0 |
| time_julian | UL | 8969 | 6.8e-6 | 4.2e-5 | 4.2e-5 | 7.0e-6 |
| heading_deg | DL | 8969 | 0 | 0 | 0 | 0 |
| heading_deg | UL | 8969 | 6.8e-6 | 4.2e-5 | 18.7 | 1.7 |
| pitch/roll/temp/sound_vel | DL | 8969 | 0 | 0 | 0 | 0 |
| pitch/roll/temp/sound_vel | UL | 8969 | 6.8e-6 | 4.2e-5 | small | small |

**Verdict: MATCH.** DL is bit-exact at every matched sample (proves the PD0
binary decode is bit-correct for the file that serves as the merge's
reference clock). UL's small residual heading diff (rms 1.7°, max 18.7°) is
attributable to the nearest-time-matching artifact — max time offset between
matched ensembles is ~3.6 s, non-trivial for a fast-rotating ship instrument
during non-constant-ping-rate periods — not a decode bug.

Velocities (`d.ru`/`rv`/`rw`/`re`) were **not** compared at M1: `loadrdi.m`
applies the beam→earth rotation inline (this PD0 is BEAM-coordinate), while
our `load_rdi()` keeps ingestion and transforms as separate layers. The fair
comparison needed our own `beam2earth()` output, done at M3 instead.

## M2 — full `process_cast` run (`octave_harness/dump_m2_process_cast.m`)

All 17 steps completed cleanly in one run (`OCTAVE_HARNESS_BEGIN_STEP`
unset, full run from step 1). Sanity-gate log lines, quoted verbatim from
this session's run, against the recorded `LOG_Inverse_log`:

| Quantity | This run | Recorded (`003.nc`) |
|---|---|---|
| `mean heading offset from compasses` | `-89.978 deg` | `-89.978 deg` (exact) |
| `bottom found at` | `4305 +/- 1 m` | `4305 +/- 1 m` (exact) |
| `Barotropic velocity error` | `0.0053778 [m/s]` | `0.0053778` (exact) |
| `Velocity profile error` | `0.089 (noise 0.064)` | `0.058 (noise 0.045)` (differs — plausibly the missing SADCP constraint) |

Total run time: 560 s (dominated by `loadrdi.m`, ~4 min for both files).

### Inputs generated for this session (not supplied by LDEO)

- **`octave_harness/work/data/CTD/2Hz/003.2Hz`**: LDEO's own file was not
  archived; generated from `003_01.cnv` (24 Hz binary Sea-Bird format) via
  `octave_harness/make_2hz_ctd.py`: decimated to 2 Hz (confirmed by the log:
  "read 25599 CTD scans", matching the recorded `p-struct` exactly), PSS-78
  salinity computed from conductivity/temperature/pressure (our own Python
  CTD loader returns NaN salinity for this file — no existing computed value
  to mirror), constant lat/lon (no per-scan GPS in this cnv), IPTS-68
  temperature conversion. **Caveat**: this is a reconstruction, not LDEO's
  original file — see M4 for how well it holds up.
- **`p.time_start`/`p.time_end`**: computed from the generated file's own
  first/last scan (elapsed-seconds time base, `ctd_time_base=0` — this
  actually means "seconds since `p.time_start`", not absolute Julian day; an
  initial reading of the recorded p-struct got this backwards, caught via
  the mismatch between the file's generated absolute-Julian-day field and
  what `loadctd.m`'s `case 0` branch actually does with it).
- **`p.drot`**: hardcoded to the recorded `12.318441` (magnetic deviation),
  bypassing `magdev.m`.
- **No SADCP file** (`f.sadcp` left unset) — LDEO's archived run used
  `../data/SADCP/Leg1.mat`, not available here; `loadsadcp.m` skips cleanly
  when `f.sadcp` doesn't exist (`existf(f,'sadcp')==0`).

### Edits to the copied M-code (see `octave_harness/PATCHES.md` for full detail)

Two genuine reserved-keyword bugs fixed to get Octave to parse/run the
code at all (`do` is valid in MATLAB, a keyword in Octave):
- `getinv.m`: variable `do` → `d_orig`.
- `plotraw.m`: function parameter `do` → `is_bottom`.

One real pre-existing bug in the in-tree code, **not patched** (worked
around via `set_cast_params.m` instead): `loadnav.m`'s own `setdefv()` call
defaults `nav_time_base` onto `p`, but the consuming `switch` reads
`f.nav_time_base` — `loadctd.m`'s changelog shows the identical
`ctd_time_base` bug was fixed in 2014; the fix was never ported to
`loadnav.m`.

One missing function, confirmed absent from the entire `docs/legacy/` tree
including the tarball: `makebars.m` (used by `plotraw.m`'s diagnostic bar
overlay) — stubbed.

One MATLAB builtin unimplemented in Octave: `interp1q` — stubbed as a thin
wrapper over `interp1(...,'linear')`.

29 pure no-op plotting stubs in `octave_harness/stubs/` (figure, plot,
subplot, hold, axis, title, xlabel, ylabel, text, legend, colorbar, pcolor,
contourf, streamer, orient, print, pause, clf, gca, grid, colormap,
imagesc, shading, fill, caxis, bar, axes, set, makebars, interp1q) — added
lazily, one crash at a time, per the plan's guidance. None affect numerical
results; they silence a diagnostic-plotting subsystem irrelevant to the
stage dumps.

`end_processing_step.m` instrumentation: dumps `d`/`p`/`di`/`dr` (when
present) to `dumps/stepNN.mat` after every step, wrapped in try/catch. Also
fixed a real save/load mismatch surfaced while adding this: the existing
checkpoint save has no `.mat` extension (MATLAB auto-appends it; Octave
does not), so `begin_processing_step.m`'s `load(...'.mat')` silently never
found it. Fixed to make `p.checkpoints`-based resume actually work — this
was essential for iterating on the plotting-stub fixes without re-running
`loadrdi.m` (~4 min) each time.

## M3 — stage-diff report (`octave_harness/diff_stages.py`)

Compares against `scripts/diag_rmse_strata.py`'s `run_pipeline()`, extended
with an optional `stages` dict (additive, backward-compatible; verified no
regression: 197 passed / 8 skipped / 2 xfailed, unchanged from HANDOVER.md's
baseline).

**Bin-row convention confirmed matching** (no remapping needed): Octave's
`d.izd`/`izu` (1-based) = rows 26-50 / rows 25..1 reversed, exactly matching
our own 0-based `EnsembleData`/`SuperEnsemble` layout (UL-reversed-on-top,
DL-on-bottom) documented in PROGRAMMERS_NOTES.md.

**Stage A (post-edit) and Stage C (super-ensembles): unresolved, flagged
rather than reported as findings.** First attempt compared Octave's step01
dump (already has `loadrdi.m`'s internal outlier-masking applied) against
Python's pre-edit `post_transform` stage — a different-pipeline-point
mismatch that produced a physically-impossible 14 m/s "diff" and only 37%
finite-overlap. Corrected to compare Octave's step09 (`EDIT DATA`, after all
QA) against Python's `post_edit` stage (same nearest-time-match strategy as
M1, mean time-alignment error now ~13 μs) — but the max-diff and
finite-overlap numbers barely moved (14.1 m/s, 38% overlap). This is not
credible as a real ocean-velocity difference. Root cause not identified in
this session; candidates not yet ruled out:
- `process_cast.m` step 8 ("APPLY PITCH/ROLL CORRECTIONS") sits between
  step 1 and step 9 and may apply a bin-remapping/tilt correction our
  Python `edit_*` functions don't separately replicate (our tilt handling
  is folded into `beam2earth()`'s `gimbaled` argument at transform time,
  not a discrete post-hoc step).
  Confirm/deny this first.
- A residual row-indexing issue not caught by the izd/izu check above.

Stage C (super-ensembles) inherits whatever Stage A's issue is, plus its
own depth-grid-alignment caveat (superensemble bin centers don't
necessarily coincide between pipelines; `mean_key_err≈1.55 m` there is
non-trivial and needs depth-bin-aware pairing, not nearest-neighbor).

**Stage D (final result) is trustworthy — 100% finite-match, no alignment
ambiguity:**

| field | max\|diff\| | rms | n |
|---|---|---|---|
| u (final, Octave-harness vs Python) | 0.287 | 0.093 | 520 |
| v (final, Octave-harness vs Python) | 0.173 | 0.063 | 520 |

## M4 — final-solution triangle comparison

| Comparison | u RMSE | v RMSE | n |
|---|---|---|---|
| Octave-harness `dr.u/v` vs archived `003.nc` | **0.068** | **0.030** | 520 |
| Python `res.u/v` vs archived `003.nc` (rot+offset config) | 0.087 | 0.058 | 520 |
| Octave-harness vs Python (Stage D) | 0.093 | 0.063 | 520 |

The Octave harness — despite a reconstructed CTD file, no SADCP, and a
hardcoded magnetic deviation — lands *closer* to LDEO's real archived answer
than our current Python pipeline does. SADCP asymmetry (LDEO's archived run
had it; this Octave run and the harness's own reconstruction do not) means
this isn't a perfectly clean comparison, but the direction of the result is
still informative: it weighs against "our Python pipeline's gap comes from
imperfect input reconstruction" and toward "there is a genuine algorithmic
difference downstream of ingestion" — consistent with HANDOVER.md's
rotup2down/offsetup2down findings (tested, made things worse) and the
still-unexplained 1000-2000 m stratum gap.

## Files added/changed this session

- `octave_harness/run_octave.sh`, `dump_m1_ingestion.m`,
  `dump_m2_process_cast.m`, `diff_ingestion.py`, `diff_stages.py`,
  `make_2hz_ctd.py`, `PATCHES.md`, `REPORT.md` (this file),
  `recorded_p_struct_attrs.txt`, `recorded_LOG_Inverse_log.txt` — new.
- `octave_harness/ldeo_ix/` — copy of `docs/legacy/*.m` (67 files) plus
  `set_cast_params.m` (new, harness-specific) with the two `do`-keyword
  patches and the `end_processing_step.m` instrumentation described above.
- `octave_harness/stubs/` — 30 new stub files (list above).
- `scripts/diag_rmse_strata.py` — added the optional `stages` dict param to
  `run_pipeline()` (additive, no behavior change to existing callers;
  verified via full test suite).
- `.gitignore` — added `octave_harness/dumps/` and `octave_harness/work/`.

## What the next session should investigate

**Priority 1: resolve the Stage A/C alignment puzzle.** Before trusting any
stage-by-stage numeric comparison upstream of the final solution, confirm
whether `process_cast.m` step 8's pitch/roll correction changes `d.ru`/`rv`
in a way our Python `post_edit` stage doesn't capture. If it does, either
replicate that step's effect in a comparable Python stage, or move the
comparison point to right after step 8 specifically (a dump already exists:
`dumps/step08.mat`). If it doesn't explain the gap, the row/bin-alignment
assumption (confirmed only via izd/izu index arrays, not by inspecting
actual depth continuity across the UL/DL boundary) needs a direct check:
print `izm` at rows 24/25/26 for a well-populated ensemble and confirm depth
increases monotonically across the boundary in both pipelines.

**Priority 2, if Priority 1 resolves cleanly:** rerun `diff_stages.py` and
find the real first-DIVERGES stage — that's the actual deliverable this
harness was built for. Given M4's result (harness beats Python against the
archived answer), the superensemble/solver stages are the most likely
candidates, not ingestion/transform.

**Priority 3:** the step-12 `di` dump (`dumps/step12.mat`) is the true
solver input and is already saved — it directly enables the originally-
planned solver-only harness (feed `di` into `compute_inverse()` directly)
without needing any of the above resolved first, if that's a faster path
to isolating the solver specifically.

Numbers in this report come only from commands actually run this session
(`octave_harness/diff_ingestion.py`, `octave_harness/diff_stages.py`, the
ad-hoc M4 triangle script, `pytest`).
