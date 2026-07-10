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

Priorities 1-3 above (Stage A/C alignment puzzle, rerun `diff_stages.py`,
solver-only harness) were all executed this session (2026-07-06) — see the
P3, P1, and P2 sections below for the results. **Do not re-open those
questions; read their verdicts first.** In one line: the solver is
exonerated (P3), the step-8 pitch/roll hypothesis is dead (P1a), and the
first genuine DIVERGES stage is Stage A `ru`, driven by a depth-varying
`izm` (bin-depth-assignment) offset between the two pipelines (P1b/P2).

**[RESOLVED 2026-07-10 — see P4 above. The izm offset was root-caused
(missing Saunders pressure→depth conversion + Octave-side time-label
shift) and fixed; do not re-open. The remaining historical text below is
kept for the evidence trail.]**

**Priority 1 for the next session: find where the depth-varying `izm`
offset originates.** P2 quantified it precisely — mean +34.2 m, median
+31.0 m, rms 47.9 m, max 108.8 m, correlation of per-column mean diff with
per-column mean depth = −0.80 (sign convention: diff = Octave − Python with
`izm` negative-down, so the positive mean means Octave assigns bins
*shallower* on average — Octave is slightly deeper near the cast start
(column-mean diff −20.8 m) and increasingly shallower at depth (up to
~+90–105 m); the gap grows with cast depth rather than being a fixed
vertical shift, and is row-independent, ruling out a bin-index/row error)
— but did not isolate *why* it happens. Candidate mechanisms, in the order they touch the depth
calculation:
- Python's `assign_bin_depths()` in `src/ladcp/ingestion/ctd.py:350`
  vs. Octave's `docs/legacy/getdpthi.m` — compare the two functions
  line-by-line for a sound-speed-corrected bin-range difference (a
  depth-dependent, row-independent offset growing to ~2.5% of depth at
  4300 m is the classic signature of a sound-speed correction applied in
  one pipeline but not the other, or applied with a different reference
  profile).
- Pressure→depth conversion (different gravity/EOS constants or formula).
- CTD-time alignment feeding a time-varying pressure into the bin-depth
  calculation (would explain why the offset grows with elapsed cast time,
  not just nominal depth).

**Technique to reuse (one stage earlier than `solver_only.py`):** the
P3 solver-only comparison worked by feeding an Octave dump (`step12.mat`'s
`di`) directly into the Python function that would normally consume
Python's own upstream output (`compute_inverse()`), isolating that one
function from everything upstream of it. Apply the same trick one stage
earlier:
1. Load Octave's step09 dump (`octave_harness/work/dumps/step09.mat`, struct `d`, has `d.izm`
   already computed) and Octave's step12 dump (`octave_harness/work/dumps/step12.mat`, struct
   `di`) the same way `octave_harness/diff_stages.py` and
   `octave_harness/solver_only.py` already do (see their `_load()` /
   `scipy.io.loadmat` usage for the loader pattern).
2. Feed Octave's `d.izm`/`d.ru`/`d.rv`/`d.rw` (step09) into Python's
   `prepare_superensembles()` (`src/ladcp/solution/inverse.py:346`) in
   place of the Python `post_edit` stage's own `ens.izm`, and compare the
   resulting super-ensemble `z`/`ru`/`rv` against Octave's own `di` (step12)
   super-ensembles. If the Stage C gap collapses when Octave's izm is used
   as input, that confirms the depth-registration difference (not the
   super-ensemble averaging logic) is what's corrupting Stage C/D.
3. In parallel (cheaper, more direct), call `assign_bin_depths()` and
   `getdpthi.m` on identical CTD-time/pressure inputs and diff their
   outputs directly — this isolates the depth-assignment function itself
   from any downstream masking/pairing artifacts and is likely the fastest
   way to find the exact formula/constant difference.

**Housekeeping (done in the final-review fix wave, 2026-07-06):**
- `octave_harness/diff_stages.py`'s stale "~15 m" izm NOTE (P1's two-column
  estimate) has been updated to reference P2's full-array numbers (rms
  47.9 m, growing in magnitude with cast depth) so the script's own output
  doesn't understate the finding.
- The izm summary statistics quoted in P2 (mean/median/min/max/
  depth-correlation) are now reproduced by `diff_stages.py` itself (printed
  immediately after the Stage A izm row), not ad-hoc scratch code — see the
  `izm diff (oct - py): ...` line in its output.

Numbers in this report come only from commands actually run this session
(`octave_harness/diff_ingestion.py`, `octave_harness/diff_stages.py`, the
ad-hoc M4 triangle script, `octave_harness/solver_only.py`,
`octave_harness/diag_stage_a.py`, `pytest`).

## P3 — solver-only comparison (session 2026-07-06)

Script: `octave_harness/solver_only.py` (feeds the step-12 `di` dump into
`compute_inverse()`; both solvers see identical input, no SADCP).

Sanity block: `n_se=828`, `n_bins=50`, `izd (0-based): 25..49`, `izu: 24..0`,
`z` range `[-4294.3, -11.1] m` (negative-down, no assertion failure), `izm`
range `[-4498.5, 194.8] m`, median `dt` 10.9 s, 25 ensembles with finite
bottom-track velocity. All match the brief's expectations except `dt`:
the brief estimated "order of 10²-10³ s" but the actual median is 10.9 s
(one to two orders below) — consistent with 828 super-ensembles over a
~2.5-hour cast, and since Octave's own `dr` (step14) was computed from
this same `di.dt`, both solvers see the identical value, so the comparison
below is unaffected.

| Comparison | u RMSE | v RMSE | n |
|---|---|---|---|
| Python (velerr=0.05) vs Octave dr        | 0.0087 | 0.0098 | 520 |
| Python (velerr=presolve) vs Octave dr    | 0.0082 | 0.0095 | 520 |
| Python (velerr=0.05) vs archived 003.nc  | 0.0673 | 0.0287 | 520 |
| Python (velerr=presolve) vs archived 003.nc | 0.0677 | 0.0287 | 520 |

Presolve-informed velerr = 0.0933 (vs Python default 0.05).

**getinv.m's use of the presolve dr:** the incoming presolve `dr` (from
`step12.mat`, i.e. the lanarrow/shear-based presolve, not the final inverse
solution) is used in exactly two places in `getinv.m`. (1) Lines 123-130:
if `dr` exists and `dr.uerr` has any finite values, `ps.velerr` is
overridden to `medianan(dr.uerr)` (line 126) — this is the only coupling
this harness's two-config run exercises, and it's what produced the
0.0933 value above. (2) Lines 236-278, gated by `ps.dragfac>0`: if the ship
is under way and drag correction is enabled, `ut=interp1(dr.z,dr.u,-di.z)`
and `vt=interp1(dr.z,dr.v,-di.z)` (lines 240-241) are used to build
`shipdragvel = shipvel - (ut+i*vt)`, which feeds a wire-drag-corrected CTD
velocity estimate (`ctdvel`, built through lines 244-278) that is later
added as an explicit constraint row via `laindrag()` (line 518, weight
`ps.dragfac`) and also recorded as `uctd_drag`/`vctd_drag` diagnostics
(lines 574-575). For this cast, `recorded_p_struct_attrs.txt` confirms
`dragfac = 0` (also set as the `getinv.m:51` default), so pathway (2) was
**not** exercised in the run these dumps came from — the only real
presolve-dr coupling active here is the velerr override in (1). This means
the "identical input" claim in this comparison is not absolute (`ps.velerr`
is presolve-derived in the Octave run and can be matched or left at
Python's default), but the coupling actually in effect is small, single-
valued, and reproduced explicitly by the harness's second config.

**Verdict:** Python-vs-Octave-dr RMSE (0.0082-0.0098 for both velerr
configs) is roughly an order of magnitude below the full-pipeline Stage D
gap (u 0.093 / v 0.063) and far below the Python-vs-archive gap (u 0.067-
0.068 / v 0.029). **Solver exonerated.** Given identical `di` input, the
Python `compute_inverse()` and Octave `getinv.m` solvers agree closely
(u/v RMSE ~0.01 m/s); the full-pipeline gap against LDEO's archived answer
must therefore arise upstream of the solver — in super-ensemble formation
(`prepinv.m`/`prepare_superensembles`), the rotup2down/offsetup2down
alignment, or the editing stages — making P1/P2 (Stage A/C alignment and
the first-diverging stage) the critical path for closing the remaining
gap. The presolve-informed velerr config is only marginally closer to
Octave's `dr` (0.0082/0.0095 vs 0.0087/0.0098) and essentially identical
against the archive (0.0677/0.0287 vs 0.0673/0.0287) — consistent with
`dragfac=0` meaning the presolve `dr` coupling exercised here is weak.

## P1 — Stage A/C alignment puzzle (session 2026-07-06)

**Step-8 hypothesis:** RESOLVED-NO-OP. p.tiltcor = 0.0 (scalar);
step07 vs step08 d.ru/rv/rw identical (0 cells differing across all three
arrays, shape (50, 7682) each). The REPORT.md 2026-07-05 candidate "step 8
applies a correction our Python pipeline does not capture" is refuted; the
Stage A residual must come from comparison alignment or genuine
editing/masking differences (see below).

### P1b — row/mask alignment diagnostics (`diag_stage_a.py` part 2)

Comparison: Octave step09 `d.ru` (50, 7682) vs Python `post_edit` `ens.u`,
columns matched by nearest ensemble time (n=7682, mean time error 13.4 ms,
max 63.9 ms — well under the ~1 s ping interval, so column matching is
sound; note the task brief's expectation of ~13 us was off by 1000x, which
changes nothing at these ping rates).

**izm boundary check (best-populated column, octave col 497 / python col
1369):** depth is monotonic across the UL/DL boundary (numpy rows 24|25)
in BOTH pipelines — the combined-array assembly order is correct on both
sides. However, the izm *values* at the same rows and matched time differ
systematically:

    rows 22..27  octave: [-594.9 -602.8 -610.7 -627.1 -635.1 -643.0]
                 python: [-580.3 -588.3 -596.3 -612.2 -620.2 -628.2]

Octave registers every bin 14.4–14.9 m deeper than Python (bin length
~8 m, so ~1.8 bin lengths). The same offset (~15.0 m) appears at the
max-diff neighborhood near the start of the cast (cols 90–94, rows 19–23:
octave izm −34.5…−71.0 vs python −19.5…−55.4) — i.e. a roughly constant
**~15 m depth-registration offset between the two pipelines**, present in
both the UL and DL blocks.

**Row-shift scan (mean per-row rms of octave row r vs python row r+s,
both-finite cells only, rows with >100 overlapping cells):**

| shift s | UL rows 0-24 | DL rows 25-49 |
|---|---|---|
| −2 | 0.3818 | 0.2066 |
| −1 | 0.3079 | 0.2016 |
|  0 | 0.2733 | 0.1944 |
| +1 | 0.2687 | **0.1941** |
| +2 | **0.2536** | 0.1952 |

UL block flipped: rms 0.9635 (n=8689) — much worse; **UL orientation flip
ruled out**. No shift collapses the residual (best improvement is UL
s=+2, 0.2733→0.2536, ~7%; DL is flat to <1%), so **a fixed row
misalignment in the comparison is ruled out** as the cause of the Stage A
rms. The mild monotonic preference for positive s in the UL block is
consistent with the ~1.8-bin depth offset above (python row r+2 sits
nearest octave row r in depth) where near-surface shear is strong, but it
is a small effect, not the headline.

**Mask breakdown:** both-finite 38.0%, octave-only 5.9%, python-only
5.8%, both-nan 50.4%. The low both-finite fraction is dominated by cells
neither pipeline keeps (both-nan 50.4% — far bins beyond range), not by
mask disagreement (11.7% of cells total). The disagreement is strongly
row-structured:

- Rows 24–25 (the instrument-nearest bin of UL and DL): Octave keeps only
  1709/1647 of 7682 columns vs Python's 7551/7202 — Octave's editing
  masks the first bin of each looker far more aggressively.
- Outer UL rows: Octave keeps more than Python (row 19: 4587 vs 3719;
  row 20: 6267 vs 5167) — Python masks more far-range UL cells,
  directionally consistent with Python registering bins ~15 m shallower
  so its surface/sidelobe editor cuts more (plausible link, not proven
  here).

**Max-diff cell:** 14.110 m/s at row 21, matched col 92 (near-surface,
early cast). The neighborhoods show Octave with plausible velocities
(−0.147…−0.415 m/s) while Python contains −14.42 and +2.19 m/s garbage
cells its editor did not mask. **The headline "14 m/s divergence" is a
Python-side unmasked outlier, not a real velocity difference** — the
value is physically impossible and absent from Octave's array because
Octave's editing chain (loadrdi `outlier()` + `edit_data.m`) removed it.

**Verdict (decision-table rows that matched):**

1. *Named for the headline symptoms:* the 14 m/s max|diff| is an
   editing/masking-policy difference (Python retains near-surface garbage
   Octave masks); the 38% finite-overlap is dominated by both-nan (50.4%)
   plus row-structured mask disagreement (rows 24–25 and outer UL rows).
   Neither is a real 14 m/s velocity divergence.
2. *New genuine finding:* a ~15 m (~1.8 bin) depth-registration offset in
   `izm` between the pipelines, uniform across UL/DL at both inspected
   columns. This does not change `ru` at a given (bin, time) but shifts
   which depth every sample is later binned to — directly relevant to the
   superensemble (Stage C) and final-profile comparisons. Task 5 should
   add a `d.izm` vs `ens.izm` rms/offset row to Stage A.
3. *Ruled out:* UL-block flip (rms 0.9635 vs ≤0.27 unflipped), fixed row
   shift as the cause (no s collapses the residual), combined-array
   assembly error (izm monotonic across the boundary in both pipelines),
   and column/time mismatch (13.4 ms mean vs ~1 s ping interval).
4. *Still unresolved:* the residual ~0.19 (DL) / 0.27 (UL) m/s rms on
   both-finite cells at s=0 — not explained by any alignment artifact
   tested here. Single most informative next measurement: full-array
   `d.izm − ens.izm` statistics over matched columns to confirm the
   ~15 m offset is constant and trace its source (CTD-time alignment,
   distance-to-first-bin / blank handling, or surface offset); in
   parallel Task 5 should report mask-disagreement % as its own Stage A
   metric and keep rms restricted to both-finite cells.

(Test suite not run this session: no tracked pipeline code under `src/`
or `scripts/` was touched — only the diagnostic script
`diag_stage_a.py` and this report.)

## P4 — izm root cause found and fixed (session 2026-07-10)

**The depth-varying izm offset is root-caused and closed.** Two mechanisms,
both in `loadctd.m`, both previously absent from Python
(`octave_harness/diag_izm_root_cause.py` reproduces this end-to-end):

1. **Pressure→depth formula (dominant, the depth-correlated component).**
   `loadctd.m::p2z` is Saunders & Fofonoff (1976) with the cast latitude
   (`p.poss(1)` = −15). Python's `assign_bin_depths()` had the same formula
   in its `lat_deg` branch, but **no caller passed `lat_deg`**, so every run
   used the fallback `z = p*1.00445` — ~89 m too deep at this cast's
   ~4400 dbar bottom. Variant test: baseline reproduces P2 exactly (mean
   +34.2 m, rms 48.0, corr −0.802); switching to Saunders alone gives mean
   0.03 m, corr 0.004.
2. **Time registration (the ~±25 m depth-uncorrelated residual).** The
   remaining variant-B residual correlates 0.95 with descent rate — a pure
   time offset (fitted τ = +21 s). This is Octave's `loadctd.m` label shift
   (its besttlag lagdt, measured −23.24 s from the dumps by heading
   content-match between step03 and step09), NOT a physical correction
   Python was missing: the physical CTD(cnv)-ADCP clock offset in Python's
   data is only **−0.5 s** (measured by the new `estimate_ctd_adcp_lag()`,
   corr 0.96). Octave's ctdtimoff (−23.08 s on CTD time) and lagdt
   (−23.5 s on ADCP time) mostly cancel physically — both pipelines were
   already near-correctly registered; their *labels* differ by ~23 s.

**Comparison-methodology consequence (fixed in `diff_stages.py`):** this
cast's staggered ping pattern (1.33/1.58 s) repeats every 2.91 s, and
23.24 s ≈ 16 pings almost exactly — so nearest-time matching had been
pairing ensembles **16 pings apart with only ~13 ms apparent error**. P1's
"column/time mismatch ruled out (13.4 ms)" was fooled by this. Stage A now
measures and undoes the label shift via heading content-match.

**Pipeline fixes (commits `5726afa`, `035a475`, `fc60f20`):**
`pressure_to_depth()` (Saunders, matches p2z's 9712.654 m check value),
`estimate_ctd_adcp_lag()` (besttlag whole-series equivalent, resamples the
quantized 24 Hz cnv time base onto a ≥0.5 s grid), `assign_bin_depths()`
gains `time_offset_days`; `run_pipeline()` + the integration fixture pass
`lat_deg` and the measured lag and shift `ens.time_jul` like
`loadctd.m:443`.

**Measured effect (diff_stages.py rerun):**

| metric | before | after |
|---|---|---|
| izm rms / mean / corr(depth) | 47.9 m / +34.2 m / −0.80 | **1.55 m / +0.02 m / +0.05** |
| Stage A ru rms | 0.2146 | **0.1172** |
| Stage A rv rms | 0.2087 | 0.1250 |
| Stage A rw rms | 0.4653 | 0.1867 |
| Stage C ru/rv rms | 0.1118 / 0.1076 | 0.0980 / 0.1013 |
| Stage D u/v rms | 0.0933 / 0.0630 | 0.1031 / 0.0755 |

izm's residual (max 7.2 m, row-dependent) is consistent with the
sound-speed bin-length scaling `getdpthi.m:428-439` applies and Python does
not (yet) — a few metres at far bins.

**End-to-end RMSE vs archived 003.nc moved mixed** (NEW config: u TOTAL
0.0678 → 0.0848, v 0.0573 → 0.0595; 3000–4500 m u improves 0.0720 → 0.0464,
1000–2000 m u worsens 0.1034 → 0.1540): with physically-correct depths the
profile re-registers vertically against the still-open Stage A velocity
divergence — the old total partially benefited from error cancellation
between wrong depth registration and wrong velocity structure. The honest
open problem is now solely the Stage A editing/masking + residual velocity
difference (rms ~0.12 on both-finite cells; max-diff cells are Python-side
unmasked garbage, e.g. 14 m/s near-surface values Octave's
loadrdi-`outlier()`/`edit_data.m` chain removes).

**P4 continued (same session): editing port executed.** Item 1 below was
done immediately (commit `25df9de`): `edit_outliers()` (loadrdi
`outlier()` port — per-5-min-block, two-sweep 4σ/3σ, DL/UL independent,
NaNs velocities not just weights, bottom-track/hbot rejection included)
and `edit_mask_bins()` (edit_data.m bin masking; pipeline masks bin 0 of
each instrument because both have zero blanking distance). Measured
effect: Stage A ru rms 0.117 → **0.080** (session cumulative 0.215 →
0.080), max|diff| 14.3 → **3.0 m/s** (the Python-side garbage cells are
gone), mask_disagree 9.4 → 7.4%; Stage C ru rms 0.098 → 0.094; archive
RMSE (NEW config) u 0.0848 → **0.0787**, v 0.0595 → **0.0552**.

**What the next session should investigate (P4 handoff, updated):**
1. ~~Port `loadrdi.m::outlier()` + bin-1 masking~~ — DONE, see above.
2. Sound-speed corrections Python lacks: velocity scaling `ss/sv`
   (`getdpthi.m:182-207`, needs a `sounds.m` port + CTD temp at the
   instrument) and bin-length scaling for izm (`getdpthi.m:428-439`,
   explains the remaining izm ±7 m row-dependent residual).
3. The residual Stage A rms (~0.08 u/v on both-finite cells, max|diff|
   ~3 m/s) — remaining candidates: 3-beam solutions (Octave computes
   14422 DL / 8473 UL 3-beam solutions where one beam is bad; Python's
   beam2earth() has no 3-beam path, so those cells either differ or get
   garbage from a bad beam), and the remaining mask policy differences
   (rows 24-25 instrument-nearest-bin masking, 7.4% mask_disagree).
4. Re-measure Stage C/D and archive RMSE after (2)+(3); the 1000-2000 m
   u stratum (0.1414) is still the dominant archive-RMSE contributor.

## P2 — stage-diff rerun with corrected methodology (session 2026-07-06)

Methodology changes vs 2026-07-05: **masking-policy branch applied** (P1's
decision-table row 1) — `_row_diff` keeps its existing both-finite
restriction on rms/max|diff|, and now also prints `mask_disagree=`, the %
of matched cells where the two pipelines disagree on finite-vs-masked, so
masking differences are visible separately from velocity differences. No
row/column remapping was applied — P1 ruled that out (UL flip, fixed row
shifts). A new Stage A row was added diffing `d9.izm` vs Python's post-edit
`ens.izm` (Task 4's/P1's named next measurement, same time-key matching as
the velocity fields), to quantify the ~15 m depth-registration offset P1
found at two sampled columns.
Stage C now pairs super-ensembles only within dz_se/2 = **2.51 m**
(dz_se = 5.01 m, the median super-ensemble depth spacing) — **648/828**
pairs kept for all four Stage C fields (ru/rv/rw/weight all share the same
z-based pairing, so the kept-count is identical across them).

Full rerun output (`uv run python octave_harness/diff_stages.py`):

```
--- Stage A: post-edit (masking + sidelobe + large-vel + w-outlier QA) ---
ru (east vel, all bins)     time      max|diff|=     14.11  rms=    0.2146  %finite-match= 38.0  mean_key_err=1.55e-07  mask_disagree= 11.6%  [DIVERGES]
rv (north vel, all bins)    time      max|diff|=     11.56  rms=    0.2087  %finite-match= 38.0  mean_key_err=1.55e-07  mask_disagree= 11.6%  [DIVERGES]
rw (vert vel, all bins)     time      max|diff|=     5.153  rms=    0.4653  %finite-match= 38.0  mean_key_err=1.55e-07  mask_disagree= 11.6%  [DIVERGES]
weight (post-edit)          time      max|diff|=    0.9559  rms=    0.1679  %finite-match= 41.5  mean_key_err=1.55e-07  mask_disagree= 51.1%  [DIVERGES]
    NOTE: izm units are metres of depth, not m/s -- DIVERGES here means the ~15 m registration offset (P1), not a velocity gap.
izm (depth registration, all bins)       max|diff|=     108.8  rms=     47.91  %finite-match=100.0  mean_key_err=1.55e-07  mask_disagree=  0.0%  [DIVERGES]

--- Stage C: super-ensembles (the solver's actual input) ---
    median super-ensemble depth spacing dz_se = 5.01 m
    [ru (super-ens east vel)] depth-tolerance 2.51: kept 648/828 pairs
ru (super-ens east vel)     depth     max|diff|=    0.9057  rms=    0.1118  %finite-match= 38.4  mean_key_err=1.09  mask_disagree=  9.5%  [DIVERGES]
    [rv (super-ens north vel)] depth-tolerance 2.51: kept 648/828 pairs
rv (super-ens north vel)    depth     max|diff|=    0.7186  rms=    0.1076  %finite-match= 38.4  mean_key_err=1.09  mask_disagree=  9.5%  [DIVERGES]
    [rw (super-ens vert vel)] depth-tolerance 2.51: kept 648/828 pairs
rw (super-ens vert vel)     depth     max|diff|=     2.339  rms=     1.228  %finite-match= 38.4  mean_key_err=1.09  mask_disagree=  9.5%  [DIVERGES]
    [weight (super-ens)] depth-tolerance 2.51: kept 648/828 pairs
weight (super-ens)          depth     max|diff|=     0.901  rms=    0.1839  %finite-match= 41.3  mean_key_err=1.09  mask_disagree= 55.9%  [DIVERGES]

--- Stage D: final inverse solution (u/v profile) ---
u (final east vel)          depth     max|diff|=     0.287  rms=    0.0933  %finite-match=100.0  mean_key_err=2.5  mask_disagree=  0.0%  [DIVERGES]
v (final north vel)         depth     max|diff|=     0.173  rms=   0.06298  %finite-match=100.0  mean_key_err=2.5  mask_disagree=  0.0%  [DIVERGES]

--- Summary ---
stage/field                    max|diff|         rms    %match     verdict
ru (east vel, all bins)            14.11      0.2146     38.0%    DIVERGES
rv (north vel, all bins)           11.56      0.2087     38.0%    DIVERGES
rw (vert vel, all bins)            5.153      0.4653     38.0%    DIVERGES
weight (post-edit)                0.9559      0.1679     41.5%    DIVERGES
izm (depth registration, all bins)       108.8       47.91    100.0%    DIVERGES
ru (super-ens east vel)           0.9057      0.1118     38.4%    DIVERGES
rv (super-ens north vel)          0.7186      0.1076     38.4%    DIVERGES
rw (super-ens vert vel)            2.339       1.228     38.4%    DIVERGES
weight (super-ens)                 0.901      0.1839     41.3%    DIVERGES
u (final east vel)                 0.287      0.0933    100.0%    DIVERGES
v (final north vel)                0.173     0.06298    100.0%    DIVERGES

FIRST DIVERGES: ru (east vel, all bins)
FIRST DIVERGES (velocity fields only, excl. izm depth-registration): ru (east vel, all bins)
```

**Regression sentinel check:** Stage D u rms = 0.0933, v rms = 0.06298 —
matches the pre-refactor 0.093/0.063 (P1/2026-07-05 baseline). The
`_row_diff` refactor did not disturb Stage A/D behavior.

**izm full-array result (Task 4's named next measurement):** over all
50 rows × 7682 time-matched columns (100% finite-match — `izm` is a
coordinate field, not independently masked, so mask_disagree=0% here),
`d9.izm − ens_pe.izm` has **mean +34.2 m, median +31.0 m, rms 47.9 m,
max|diff| 108.8 m, per-cell min diff −24.2 m** (this per-cell min is now
reproduced directly by `diff_stages.py`'s `izm diff (oct - py): ...` line,
printed immediately after the Stage A izm row; an earlier draft of this
paragraph quoted "−20.8 m" for the min, which was actually the min of
*per-column-mean* diff, not the per-cell min — see below). This is larger
and more variable than P1's two-column estimate (~14.4–15.0 m) — the
2026-07-06 follow-up investigation (per-row and per-column breakdown, not
committed to `diff_stages.py` since it is a one-off diagnostic breakdown
beyond Fix 3's mean/median/min/max/corr scope) found the offset is
**uniform across all 50 rows** (row-independent, ruling out a row/bin-index
error) but **varies strongly with time/cast-depth**: column-mean diff
ranges from about **−20.8 m** near the start of the cast to ~+90–105 m near
the cast's maximum depth (correlation of per-column mean diff with
per-column mean depth = −0.80; the per-cell correlation now printed by
`diff_stages.py`, −0.800, coincides because the offset is row-uniform;
sign convention: `izm` is negative-down, so this correlation means the
offset grows in magnitude as the cast gets deeper). It is not a clean proportional (percentage-of-depth)
scaling either (diff/|depth| ratios range roughly −0.03 to +0.08 across the
depth range, not a single constant fraction). **Net effect on the P1
"~15 m offset" characterization:** confirmed as a real, non-trivial,
depth/time-varying depth-registration difference between the two
pipelines — larger overall (mean ~34 m) than the two sampled columns
suggested, growing with cast depth rather than being a fixed vertical
shift. Root cause (CTD-time alignment, distance-to-first-bin/blank
handling, or a sound-speed/pressure-integration difference that
compounds with depth) is still not isolated — flagged for a future
session, not fixed here per this task's comparison-only scope.

**FIRST DIVERGES: ru (east vel, all bins), Stage A.** The `izm`
depth-registration row is technically the largest-magnitude DIVERGES
entry in the table (rms 47.9), but its units are metres of depth, not
m/s, and it is excluded from the velocity-only tally the harness now
prints separately: **FIRST DIVERGES (velocity fields only) is also `ru`
(east vel, all bins), Stage A** (rms 0.215, max|diff| 14.11). The
`mask_disagree` columns show this isn't primarily a masking-policy
artifact either: only 11.6% of Stage A velocity cells disagree on
finite-vs-masked (P1's "38% finite-overlap" is dominated by both-nan
cells neither pipeline keeps, as P1 established), and the max|diff| cell
is itself a both-finite cell (P1 traced it to a Python-side unmasked
near-surface outlier) — so the residual 0.19–0.27 m/s Stage A rms on
clean, both-finite cells (P1's "still unresolved" item 4) persists after
this session's fixes and is the honest first velocity divergence.

Read together with the P3 solver-only result (**solver exonerated**: given
identical `di` input, Python-vs-Octave-dr RMSE is 0.008–0.010 m/s, roughly
10× below the full-pipeline Stage D gap and far below the Python-vs-archive
gap), the Python pipeline's gap against LDEO is attributed to: **Stage A
edit/masking policy plus the newly-quantified depth-registration
difference, both acting before the solver ever runs.** The solver
(`compute_inverse()`/`getinv.m`) is not the source of the gap — P3 showed
it reproduces Octave's answer closely given identical input. The gap
instead accumulates upstream: Stage A already shows a genuine ~0.2 m/s rms
velocity divergence on both-finite cells (not explained by masking
disagreement or comparison alignment — P1 ruled out row shifts/UL flip),
and this session's izm measurement shows the two pipelines don't even
agree on which physical depth a given bin/time cell represents (offset
growing from near 0 to ~100 m over the cast). Because a depth-registration
difference of that size changes which raw cells land in which
super-ensemble depth bin, it plausibly explains why Stage C's rms (0.11
u / 0.11 v, now measured on depth-tolerance-filtered pairs so it isn't a
comparison artifact) and Stage D's rms (0.093 u / 0.063 v) don't shrink
back toward P3's solver-only ~0.01: the solver is being fed a
super-ensemble input that itself already differs from Octave's, because
the two pipelines assign raw ADCP bins to different depths before
averaging. The single highest-leverage next step is isolating *why* `izm`
diverges (CTD time-base alignment, distance-to-first-bin/blank handling,
or sound-speed/pressure integration) — that is the one remaining
unexplained mechanism standing between the solver-exoneration result and
a fully closed gap.
