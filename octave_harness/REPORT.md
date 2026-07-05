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
