# Octave Harness Follow-up (P3 → P1 → P2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Execute the three "next session" priorities from `octave_harness/REPORT.md` — a solver-only differential harness (its P3), resolution of the Stage A/C alignment puzzle (its P1), and a rerun of the stage diff to name the true first-DIVERGES stage (its P2).

**Architecture:** Everything here is *diagnostic* work living in `octave_harness/`. No Octave/Docker runs are needed — all 17 per-step `.mat` dumps from the 2026-07-05 session already exist in `octave_harness/work/dumps/`. Each task loads those dumps with `scipy.io.loadmat`, runs pure-Python comparisons against `src/ladcp/`, and records findings in `octave_harness/REPORT.md`.

**Task order rationale:** REPORT.md numbers these Priority 1/2/3, but this plan runs them **P3 first**: the solver-only harness is fully independent (needs nothing resolved first), and its result is the highest-information experiment — it directly tests the prime suspect (the solver) identified by the M4 triangle. P1 (alignment puzzle) comes second; P2 is gated on P1's outcome and runs last.

**Tech Stack:** Python 3.11 via `uv`, numpy, scipy.io, netCDF4. Git Bash for commands.

## Global Constraints

- **Do NOT modify anything under `src/ladcp/` or `octave_harness/ldeo_ix/`** in this plan. This is a diagnostic session; the only editable code is new/existing scripts directly under `octave_harness/` plus `REPORT.md`/`HANDOVER.md`. (`scripts/diag_rmse_strata.py` is imported, never edited.)
- **Do NOT modify anything under `docs/legacy/`** (read-only reference, per repo CLAUDE.md).
- The full test suite must remain at its baseline after every commit: **197 passed, 8 skipped, 2 xfailed**.
- Test command (from LADCP repo root, Git Bash): `TEST_DATA_DIR=test_data uv run python -m pytest`. **`uv run pytest` (without `python -m`) is broken on this machine** ("uv trampoline" error, see HANDOVER.md) — always use `uv run python -m pytest`.
- Run all commands from the LADCP repo root: `C:\Users\peter_sha\Documents\sourcecode\Nuyina\LADCP`.
- Dump files (`octave_harness/dumps/`, `octave_harness/work/`) are gitignored — never `git add -f` them.
- If a needed `.mat` dump were somehow missing, re-generating it requires Docker Desktop running and `./octave_harness/run_octave.sh octave_harness/dump_m2_process_cast.m` (~10 min). Do not do this unless a load actually fails — all 17 exist as of 2026-07-05.
- **Report numbers honestly.** Several steps below have genuinely unknown outcomes; the deliverable is the measured number and its interpretation, not a "pass". Never tweak a comparison until it looks good — every methodology choice must be written down in REPORT.md.

## Background facts (verified 2026-07-06, cite these instead of re-deriving)

- Step dumps: `octave_harness/work/dumps/step01.mat` … `step17.mat`. Each contains MATLAB structs `d`, `p`, and (from step 10 on) `di`, (from step 11 on) `dr`. Load with `sio.loadmat(path, struct_as_record=False, squeeze_me=True)`.
- `process_cast.m` step names: 7 = FIND SURFACE & SEA BED, 8 = APPLY PITCH/ROLL CORRECTIONS, 9 = EDIT DATA, 10 = FORM SUPER ENSEMBLES, 11 = REMOVE SUPER-ENSEMBLE OUTLIERS (lanarrow presolve → first `dr`), 12 = RE-FORM SUPER ENSEMBLES, 14 = CALCULATE INVERSE SOLUTION.
- Step 8 applies `uvwrot` **only if `length(p.tiltcor) > 1`** (`process_cast.m:317-327`). The default is scalar `0` (`ldeo_ix/default.m:309`), `set_cast_params.m` does not set it, and the recorded LDEO p-struct has `tiltcor = 0` (`recorded_p_struct_attrs.txt:232`). Expectation: step 8 was a no-op.
- `prepinv.m:451-475` builds `di` with fields `ru, rv, rw, ruvs, weight, izm, z, dt, time_jul, bvel (= d.bvel', so beams/components along dim 0), bvels, hbot, slat, slon, izd, izu` — everything Python's `SuperEnsemble` dataclass (`src/ladcp/solution/inverse.py:38-59`) needs. `di.dt` is already in **seconds** with the same edge-mirroring as Python (`prepinv.m:684-686` vs `inverse.py:480-486`).
- `di.izd`/`di.izu` are **1-based** MATLAB row indices (this cast: izd = 26..50, izu = 25..1 reversed). Python's are 0-based. Subtract 1.
- Octave step-14 call: `[p,dr,ps,de]=getinv(di,p,ps,dr,1)` (`process_cast.m:434`). `getinv.m` **uses the presolve `dr` it is handed**: line 125-126 sets `ps.velerr = medianan(dr.uerr)`, and lines ~240-241 interpolate `dr.u/dr.v` onto `-di.z` (purpose to be confirmed in Task 2 Step 6). Python's `compute_inverse` has no presolve input; its `InverseParams.velerr` defaults to `0.05` (`inverse.py:881`).
- Python pipeline entry: `run_pipeline(data_dir, legacy=False, rot=True, offset=True, stages=dict())` in `scripts/diag_rmse_strata.py:45` — `rot=True, offset=True` is the config matching the Octave run (see `diff_stages.py` module docstring). `stages` is populated with keys `post_transform`, `post_edit`, `superensembles`, `result`.
- Data dir: `test_data/2015_P16N/` (contains `003DL000.000`, `003UL000.000`, `003_01.cnv`, archived `003.nc`).
- Known baseline numbers from REPORT.md (2026-07-05 session):
  - Stage D (Octave harness vs Python, nearest-depth): u rms 0.093, v rms 0.063, n=520.
  - Octave harness vs archived `003.nc`: u RMSE 0.068, v RMSE 0.030.
  - Python vs archived `003.nc` (rot+offset): u RMSE 0.087, v RMSE 0.058.
  - Stage A (Octave step09 vs Python post_edit): max|diff| ≈ 14.1 m/s, ~38% finite-overlap — the unresolved puzzle.
- Git state at plan time: `octave_harness/` untracked; `.gitignore` and `scripts/diag_rmse_strata.py` modified — all from the 2026-07-05 session, intentionally kept, not yet committed.

---

### Task 1: Baseline — verify tests, commit the 2026-07-05 harness work

**Files:**
- No file edits. Commits existing work: `octave_harness/` (scripts, `ldeo_ix/`, `stubs/`, docs), `.gitignore`, `scripts/diag_rmse_strata.py`.

**Interfaces:**
- Produces: a clean git baseline so every later task has a readable diff.

- [ ] **Step 1: Run the full test suite to confirm the baseline**

Run:
```bash
cd /c/Users/peter_sha/Documents/sourcecode/Nuyina/LADCP
TEST_DATA_DIR=test_data uv run python -m pytest -q
```
Expected: `197 passed, 8 skipped, 2 xfailed` (warnings OK). If the counts differ, STOP and report — do not commit on a changed baseline.

- [ ] **Step 2: Confirm the bulky outputs are ignored**

Run:
```bash
git check-ignore -v octave_harness/dumps/m1_loadrdi.mat octave_harness/work/dumps/step01.mat
git status --short
```
Expected: both paths matched by `.gitignore` rules (`octave_harness/dumps/`, `octave_harness/work/`); `git status` shows `octave_harness/` untracked plus the two modified files, and **no** `.mat` files listed individually.

- [ ] **Step 3: Commit**

```bash
git add octave_harness .gitignore scripts/diag_rmse_strata.py
git commit -m "feat: Octave differential harness (M1-M4 session 2026-07-05)

Ingestion diff, full 17-step process_cast run under Octave 9.2 (Docker),
per-step .mat dumps, stage diff vs Python pipeline. Findings in
octave_harness/REPORT.md; M-code patches documented in PATCHES.md."
```
Expected: commit succeeds; `git status --short` afterwards shows only this plan file (if not yet committed) — commit it too if present:
```bash
git add docs/superpowers/plans/2026-07-06-octave-harness-followup.md
git commit -m "docs: plan for solver-only harness + stage A/C alignment follow-up"
```

---

### Task 2: Solver-only harness (REPORT.md Priority 3 — run FIRST)

Feed the exact `di` struct that Octave's `getinv.m` received (the step-12 dump) into Python's `compute_inverse()`. Both solvers then see **identical input**, so any output difference is attributable to the solver implementations alone.

**Files:**
- Create: `octave_harness/solver_only.py`
- Modify: `octave_harness/REPORT.md` (append results section)

**Interfaces:**
- Consumes: `SuperEnsemble`, `InverseParams`, `compute_inverse` from `src/ladcp/solution/inverse.py` (signatures in Background facts); dumps `step12.mat` (input `di` + presolve `dr`), `step14.mat` (Octave's inverse-solution `dr`).
- Produces: printed comparison table + REPORT.md section "P3: solver-only comparison". Task 6 reads that section.

- [ ] **Step 1: Write the script**

Create `octave_harness/solver_only.py`:

```python
"""Solver-only differential harness (REPORT.md Priority 3).

Feeds the exact solver input LDEO_IX's getinv.m received -- the `di`
struct dumped after process_cast.m step 12 (RE-FORM SUPER ENSEMBLES) --
into our Python compute_inverse(). Both solvers see IDENTICAL input, so
output differences isolate the solver implementations (getinv.m vs
compute_inverse), with one known coupling: getinv.m tunes ps.velerr from
the lanarrow presolve (getinv.m:125-126), so we run Python both with the
default velerr and with the presolve-informed value.

Run from the LADCP repo root:  uv run python octave_harness/solver_only.py
"""
from pathlib import Path
import sys

import netCDF4
import numpy as np
import scipy.io as sio

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from ladcp.solution.inverse import (  # noqa: E402
    InverseParams,
    SuperEnsemble,
    compute_inverse,
)

DUMPS = REPO / "octave_harness" / "work" / "dumps"
DATA_DIR = REPO / "test_data" / "2015_P16N"


def _load(step: int):
    return sio.loadmat(
        DUMPS / f"step{step:02d}.mat", struct_as_record=False, squeeze_me=True
    )


def _f(x):
    return np.asarray(x, dtype=float)


def build_se(di) -> SuperEnsemble:
    """Map the Octave di struct (prepinv.m:451-475) onto SuperEnsemble."""
    bvel = np.atleast_2d(_f(di.bvel))
    bvels = np.atleast_2d(_f(di.bvels))
    n_se = _f(di.z).size
    # prepinv stores d.bvel' -- components along dim 0 -- but verify, and
    # accept either orientation. Rows are u, v, w[, err]; keep first 3.
    if bvel.shape[0] == n_se and bvel.shape[1] in (3, 4):
        bvel, bvels = bvel.T, bvels.T
    assert bvel.shape[0] in (3, 4) and bvel.shape[1] == n_se, bvel.shape
    izd = np.asarray(di.izd, dtype=int).ravel() - 1  # MATLAB 1-based
    izu = np.asarray(di.izu, dtype=int).ravel() - 1
    return SuperEnsemble(
        ru=_f(di.ru), rv=_f(di.rv), rw=_f(di.rw),
        ruvs=_f(di.ruvs), weight=_f(di.weight),
        izm=_f(di.izm), z=_f(di.z), dt=_f(di.dt),
        time_jul=_f(di.time_jul),
        bvel=bvel[:3, :], bvels=bvels[:3, :],
        hbot=_f(di.hbot), slat=_f(di.slat), slon=_f(di.slon),
        izd=izd, izu=izu,
    )


def sanity(se: SuperEnsemble) -> None:
    print(f"n_se={se.z.size}  n_bins={se.ru.shape[0]}  n_cols={se.ru.shape[1]}")
    print(f"izd (0-based): {se.izd.min()}..{se.izd.max()}   "
          f"izu: {se.izu[0]}..{se.izu[-1]}")
    print(f"z range: [{np.nanmin(se.z):.1f}, {np.nanmax(se.z):.1f}] m "
          f"(expect <= 0, negative-down)")
    print(f"izm range: [{np.nanmin(se.izm):.1f}, {np.nanmax(se.izm):.1f}] m")
    print(f"dt: median {np.nanmedian(se.dt):.1f} s")
    print(f"bvel finite ensembles: {np.isfinite(se.bvel[0]).sum()}")
    assert se.ru.shape == se.izm.shape == se.weight.shape
    assert se.ru.shape[1] == se.z.size == se.dt.size
    assert np.nanmax(se.z) <= 0, "di.z not negative-down -- investigate before proceeding"
    assert se.izd.min() >= 0 and se.izd.max() < se.ru.shape[0]


def compare(tag, z_a, u_a, v_a, z_b, u_b, v_b):
    """Interp b onto a's depth grid (both positive-down), print RMSE."""
    ub = np.interp(z_a, z_b, u_b, left=np.nan, right=np.nan)
    vb = np.interp(z_a, z_b, v_b, left=np.nan, right=np.nan)
    ok = np.isfinite(u_a) & np.isfinite(ub)
    ru = float(np.sqrt(np.mean((u_a[ok] - ub[ok]) ** 2)))
    rv = float(np.sqrt(np.mean((v_a[ok] - vb[ok]) ** 2)))
    mu = float(np.max(np.abs(u_a[ok] - ub[ok])))
    print(f"{tag:<52} u_rmse={ru:.4f}  v_rmse={rv:.4f}  u_max={mu:.4f}  n={ok.sum()}")
    return ru, rv


def main() -> None:
    step12 = _load(12)
    di, dr_presolve = step12["di"], step12["dr"]
    dr_oct = _load(14)["dr"]  # Octave's inverse solution (step 14)

    se = build_se(di)
    sanity(se)

    ds = netCDF4.Dataset(DATA_DIR / "003.nc")
    u_ship, v_ship = float(ds.uship), float(ds.vship)
    ref_z = np.array(ds.variables["z"][:])
    ref_u = np.array(ds.variables["u"][:])
    ref_v = np.array(ds.variables["v"][:])
    ds.close()

    velerr_presolve = float(np.nanmedian(_f(dr_presolve.uerr)))
    print(f"\npresolve-informed velerr (median dr.uerr, getinv.m:125): "
          f"{velerr_presolve:.4f}  (Python default: 0.05)")

    oct_z, oct_u, oct_v = _f(dr_oct.z), _f(dr_oct.u), _f(dr_oct.v)

    print("\n--- solver-only comparison (identical di input, no SADCP) ---")
    for tag, params in [
        ("velerr=default(0.05)", InverseParams()),
        (f"velerr=presolve({velerr_presolve:.4f})",
         InverseParams(velerr=velerr_presolve)),
    ]:
        res = compute_inverse(se, params=params, u_ship=u_ship, v_ship=v_ship)
        compare(f"Python[{tag}] vs Octave dr (same di)",
                oct_z, oct_u, oct_v, res.z, res.u, res.v)
        compare(f"Python[{tag}] vs archived 003.nc",
                ref_z, ref_u, ref_v, res.z, res.u, res.v)

    print("\ncontext: full-pipeline Stage D gap (REPORT.md) was "
          "u_rmse 0.093 / v_rmse 0.063;")
    print("Octave-vs-archive was u 0.068 / v 0.030; Python-vs-archive u 0.087 / v 0.058.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it — sanity block first**

Run:
```bash
cd /c/Users/peter_sha/Documents/sourcecode/Nuyina/LADCP
uv run python octave_harness/solver_only.py
```
Expected sanity output (values known from the 2026-07-05 session): `n_bins=50`, `izd (0-based): 25..49`, `izu: 24..0`, z range negative (roughly −4300..0 m), median dt on the order of 10²–10³ s, and no assertion failures. If the `z <= 0` assertion fails, print the actual sign convention, negate `z`/`izm` consistently, note the deviation in REPORT.md, and re-run.

- [ ] **Step 3: Record the comparison numbers**

The four `compare()` lines print. There is no pre-known "correct" value — record all of them verbatim. Interpretation guide (write the matching conclusion into REPORT.md in Step 5):

| Observation | Conclusion |
|---|---|
| Python-vs-Octave-dr rms ≪ 0.093 (e.g. < 0.02) for the better velerr config | **Solver exonerated.** The full-pipeline gap arises upstream: super-ensemble formation, rotup2down/offsetup2down, or editing. P1/P2 become the critical path. |
| Python-vs-Octave-dr rms ≈ 0.09 (same order as Stage D) | **Solver implicated.** The divergence lives in `getinv.m` vs `compute_inverse` (constraint weighting, presolve coupling, lsq details). |
| Presolve-velerr config much closer than default | The presolve coupling matters; note that Python has no equivalent of `getinv.m`'s dr-based tuning. |

- [ ] **Step 4: Check what else getinv.m does with the presolve dr**

Read `octave_harness/ldeo_ix/getinv.m` lines 110–260 (already known: 125-126 sets `ps.velerr`; ~240-241 computes `ut=interp1(dr.z,dr.u,-di.z)`). Determine what `ut`/`vt` feed (weight down-weighting? outlier rejection? a prior?) and write one paragraph in REPORT.md describing every use of the incoming `dr`, with line numbers. This bounds how "identical-input" the comparison truly is.

- [ ] **Step 5: Append results to REPORT.md**

Append to `octave_harness/REPORT.md` a section:

```markdown
## P3 — solver-only comparison (session 2026-07-06)

Script: `octave_harness/solver_only.py` (feeds the step-12 `di` dump into
`compute_inverse()`; both solvers see identical input, no SADCP).

| Comparison | u RMSE | v RMSE | n |
|---|---|---|---|
| Python (velerr=0.05) vs Octave dr        | <measured> | <measured> | <n> |
| Python (velerr=presolve) vs Octave dr    | <measured> | <measured> | <n> |
| Python (velerr=0.05) vs archived 003.nc  | <measured> | <measured> | <n> |
| Python (velerr=presolve) vs archived 003.nc | <measured> | <measured> | <n> |

Presolve-informed velerr = <measured> (vs Python default 0.05).

**getinv.m's use of the presolve dr:** <paragraph from Step 4, with line numbers>

**Verdict:** <one of the three interpretation-guide conclusions, stated with
its supporting numbers>
```
Replace every `<measured>` with the actual printed number — leaving a placeholder is a task failure.

- [ ] **Step 6: Verify no regressions, commit**

Run:
```bash
TEST_DATA_DIR=test_data uv run python -m pytest -q
```
Expected: `197 passed, 8 skipped, 2 xfailed`.

```bash
git add octave_harness/solver_only.py octave_harness/REPORT.md
git commit -m "feat: solver-only differential harness (P3) -- feed Octave di dump into compute_inverse"
```

---

### Task 3: Stage A puzzle, part 1 — kill or confirm the step-8 hypothesis (Priority 1a)

REPORT.md's leading hypothesis for the Stage A mess is that step 8 (APPLY PITCH/ROLL CORRECTIONS) changes `d.ru/rv` in a way Python doesn't replicate. The Background facts strongly predict step 8 was a **no-op** (`tiltcor` is scalar 0). Prove it either way with a bit-diff.

**Files:**
- Create: `octave_harness/diag_stage_a.py` (part 1 of 2; Task 4 extends this file)
- Modify: `octave_harness/REPORT.md`

**Interfaces:**
- Consumes: dumps `step07.mat`, `step08.mat`.
- Produces: `check_step8_noop()` in `octave_harness/diag_stage_a.py`; a REPORT.md subsection "P1 findings" that Task 4 appends to. Task 4 adds `main()` wiring — keep part 1 callable standalone via `python -c` as shown below.

- [ ] **Step 1: Write the check**

Create `octave_harness/diag_stage_a.py`:

```python
"""Stage A/C alignment diagnostics (REPORT.md Priority 1).

Part 1 (Task 3): is process_cast.m step 8 (APPLY PITCH/ROLL CORRECTIONS)
a no-op on this cast? process_cast.m:317 gates it on length(p.tiltcor)>1;
default.m:309 defaults tiltcor to scalar 0 and the recorded LDEO p-struct
agrees (recorded_p_struct_attrs.txt:232) -- so we expect step07 == step08
bit-for-bit.

Part 2 (Task 4): row/column alignment diagnostics for the 14 m/s Stage A
residual -- see the functions below check_step8_noop().

Run from the LADCP repo root:  uv run python octave_harness/diag_stage_a.py
"""
from pathlib import Path
import sys

import numpy as np
import scipy.io as sio

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "octave_harness"))

DUMPS = REPO / "octave_harness" / "work" / "dumps"
DATA_DIR = REPO / "test_data" / "2015_P16N"


def _load(step: int):
    return sio.loadmat(
        DUMPS / f"step{step:02d}.mat", struct_as_record=False, squeeze_me=True
    )


def _f(x):
    return np.asarray(x, dtype=float)


def check_step8_noop() -> bool:
    s7, s8 = _load(7), _load(8)
    d7, d8, p8 = s7["d"], s8["d"], s8["p"]
    print(f"p.tiltcor = {p8.tiltcor!r}  (scalar => step 8 gate is False)")
    all_identical = True
    for name in ("ru", "rv", "rw"):
        a, b = _f(getattr(d7, name)), _f(getattr(d8, name))
        same = (np.isnan(a) & np.isnan(b)) | (a == b)
        ident = bool(same.all())
        all_identical &= ident
        print(f"d.{name}: shapes {a.shape}=={b.shape}  identical={ident}  "
              f"n_cells_differing={(~same).sum()}")
    print(f"\nSTEP-8 NO-OP: {all_identical}")
    return all_identical


if __name__ == "__main__":
    check_step8_noop()
```

- [ ] **Step 2: Run it**

Run:
```bash
cd /c/Users/peter_sha/Documents/sourcecode/Nuyina/LADCP
uv run python octave_harness/diag_stage_a.py
```
Expected: `p.tiltcor = 0` (scalar) and `STEP-8 NO-OP: True`. If instead cells differ: report which rows (UL block 0-24 vs DL block 25-49) and how many columns are affected, and record in REPORT.md that step 8 is real on this cast — Task 4 must then also compare against Python with this correction in mind.

- [ ] **Step 3: Record in REPORT.md**

Append to `octave_harness/REPORT.md` (start the P1 section; Task 4 extends it):

```markdown
## P1 — Stage A/C alignment puzzle (session 2026-07-06)

**Step-8 hypothesis:** <RESOLVED-NO-OP | CONFIRMED-ACTIVE>. p.tiltcor = <value>;
step07 vs step08 d.ru/rv/rw <identical | differ in N cells, detail...>.
<If no-op:> The REPORT.md 2026-07-05 candidate "step 8 applies a correction
our Python pipeline doesn't capture" is refuted; the Stage A residual must
come from comparison alignment or genuine editing/masking differences
(see below).
```

- [ ] **Step 4: Commit**

```bash
git add octave_harness/diag_stage_a.py octave_harness/REPORT.md
git commit -m "diag: step-8 pitch/roll correction is a no-op on P16N 003 (P1a)"
```
(Adjust the message if Step 2 found step 8 active.)

---

### Task 4: Stage A puzzle, part 2 — row/mask alignment diagnostics (Priority 1b)

The Stage A comparison (Octave step09 `d.ru` vs Python `post_edit` u) showed max|diff| ≈ 14 m/s at only ~38% finite-overlap — physically impossible as an ocean-velocity difference. Localize it: row misalignment, mask disagreement, or genuine data difference.

**Files:**
- Modify: `octave_harness/diag_stage_a.py` (append part 2)
- Modify: `octave_harness/REPORT.md`

**Interfaces:**
- Consumes: dump `step09.mat`; `run_pipeline(DATA_DIR, legacy=False, rot=True, offset=True, stages=...)` from `scripts/diag_rmse_strata.py` (stages key `post_edit` → `EnsembleData` with `.u/.v/.izm/.time_jul`, rows 0-24 = UL reversed, 25-49 = DL); `_nearest_match(key_a, key_b)` from `octave_harness/diff_stages.py` (returns `(indices_into_b, abs_key_err)`).
- Produces: diagnostic output + REPORT.md P1 conclusions consumed by Task 5.

- [ ] **Step 1: Append part 2 to `octave_harness/diag_stage_a.py`**

Add below `check_step8_noop()` (and replace the `__main__` block):

```python
def _load_stage_a():
    """Octave step09 d + Python post_edit ens, columns matched by time."""
    from diag_rmse_strata import run_pipeline
    from diff_stages import _nearest_match

    d9 = _load(9)["d"]
    stages: dict = {}
    run_pipeline(DATA_DIR, legacy=False, rot=True, offset=True, stages=stages)
    ens = stages["post_edit"]
    idx, key_err = _nearest_match(_f(d9.time_jul), ens.time_jul)
    print(f"column match: n={idx.size}  mean_time_err={key_err.mean()*86400:.3g} s  "
          f"max={key_err.max()*86400:.3g} s")
    return d9, ens, idx


def izm_boundary_check(d9, ens, idx) -> None:
    """REPORT.md's requested direct check: depth continuity across the
    UL/DL row boundary (numpy rows 24|25) in BOTH pipelines."""
    oct_ru = _f(d9.ru)
    col_o = int(np.argmax(np.isfinite(oct_ru).sum(axis=0)))  # best-populated
    col_p = int(idx[col_o])
    oct_izm = _f(d9.izm)[22:28, col_o]
    py_izm = ens.izm[22:28, col_p]
    print(f"\nizm at rows 22..27, octave col {col_o} / python col {col_p}:")
    print(f"  octave: {np.array2string(oct_izm, precision=1)}")
    print(f"  python: {np.array2string(py_izm, precision=1)}")
    for name, v in (("octave", oct_izm), ("python", py_izm)):
        dv = np.diff(v[np.isfinite(v)])
        mono = bool((dv < 0).all() or (dv > 0).all())
        print(f"  {name}: monotonic across UL/DL boundary = {mono}")


def row_shift_scan(d9, ens, idx) -> None:
    """If a fixed row shift (or UL-block flip) collapses the residual, the
    Stage A 'divergence' is a comparison-indexing artifact."""
    oct_ru = _f(d9.ru)
    py_u = ens.u[:, idx]
    n_rows = oct_ru.shape[0]
    for label, rows in (("UL rows 0-24", range(0, 25)),
                        ("DL rows 25-49", range(25, n_rows))):
        print(f"\n{label}: rms(octave row r vs python row r+s)")
        for s in (-2, -1, 0, 1, 2):
            vals = []
            for r in rows:
                if 0 <= r + s < n_rows:
                    dd = oct_ru[r] - py_u[r + s]
                    fin = np.isfinite(dd)
                    if fin.sum() > 100:
                        vals.append(float(np.sqrt(np.mean(dd[fin] ** 2))))
            print(f"  s={s:+d}: mean_row_rms={np.mean(vals):.4f}  ({len(vals)} rows)")
    # UL-block orientation check: was the UL block flipped?
    ul_flip = oct_ru[:25][::-1]
    dd = ul_flip - py_u[:25]
    fin = np.isfinite(dd)
    print(f"\nUL block FLIPPED: rms={np.sqrt(np.mean(dd[fin] ** 2)):.4f} "
          f"(n={fin.sum()})")


def mask_structure(d9, ens, idx) -> None:
    """Where does the 38% finite-overlap loss come from, and where is the
    14 m/s cell?"""
    oct_ru = _f(d9.ru)
    py_u = ens.u[:, idx]
    fo, fp = np.isfinite(oct_ru), np.isfinite(py_u)
    both, oct_only, py_only = fo & fp, fo & ~fp, ~fo & fp
    print(f"\ncells: both-finite {both.mean():.1%}  octave-only {oct_only.mean():.1%}  "
          f"python-only {py_only.mean():.1%}  both-nan {(~fo & ~fp).mean():.1%}")
    print("per-row finite counts (row, octave, python, both):")
    for r in range(oct_ru.shape[0]):
        print(f"  {r:2d}  {fo[r].sum():5d}  {fp[r].sum():5d}  {both[r].sum():5d}")
    diff = np.abs(oct_ru - py_u)
    diff[~both] = -1.0
    r, c = np.unravel_index(int(np.argmax(diff)), diff.shape)
    print(f"\nmax|diff|={diff[r, c]:.3f} m/s at row {r}, matched col {c}")
    r0, r1 = max(0, r - 2), min(oct_ru.shape[0], r + 3)
    c0, c1 = max(0, c - 2), min(oct_ru.shape[1], c + 3)
    print("octave neighborhood:")
    print(np.array2string(oct_ru[r0:r1, c0:c1], precision=3))
    print("python neighborhood:")
    print(np.array2string(py_u[r0:r1, c0:c1], precision=3))
    print("octave izm neighborhood:")
    print(np.array2string(_f(d9.izm)[r0:r1, c0:c1], precision=1))
    py_izm = ens.izm[:, idx]  # column-matched, same as py_u
    print("python izm neighborhood:")
    print(np.array2string(py_izm[r0:r1, c0:c1], precision=1))


if __name__ == "__main__":
    check_step8_noop()
    d9, ens, idx = _load_stage_a()
    izm_boundary_check(d9, ens, idx)
    row_shift_scan(d9, ens, idx)
    mask_structure(d9, ens, idx)
```

- [ ] **Step 2: Run it** (takes a few minutes — `run_pipeline` loads both PD0 files)

Run:
```bash
uv run python octave_harness/diag_stage_a.py
```
Expected mechanics: column match mean_time_err ≈ 0.013 ms-order (13 μs was measured last session); both izm series monotonic is the *hoped* outcome but genuinely unknown. Record everything printed.

- [ ] **Step 3: Interpret with this decision table**

| Symptom in output | Diagnosis | Action (this plan) |
|---|---|---|
| `s=0` clearly minimal in both blocks, flipped-UL rms worse | Rows are aligned; artifact is elsewhere | Continue down this table |
| Some `s≠0` or the flipped-UL variant clearly minimal | Row/orientation misalignment in the *comparison* | Note the correct mapping; Task 5 applies it to `diff_stages.py` Stage A/C |
| izm not monotonic across boundary in ONE pipeline | Combined-array assembly differs (real finding, not artifact) | Document precisely which pipeline and which rows; do NOT change src/ — record as a genuine divergence for the next session |
| both-finite % low mainly because octave-only ≫ python-only (or vice versa) at specific rows | The two editors mask different cells (e.g. Python sidelobe/w-outlier vs Octave `edit_data.m`) — a masking-policy difference, not velocity divergence | In Task 5, restrict Stage A rms to both-finite cells (already the case) and report mask-disagreement % as its own metric |
| max-diff cell sits at a row whose izm values differ by ≥ one bin length between pipelines | Depth-registration (bin-mapping) difference | Document; compare `d.izm` vs `ens.izm` rms as an extra Stage A row in Task 5 |
| Neighborhoods show plausible velocities both sides, just different | Genuine data difference upstream of editing | Stage A genuinely DIVERGES — carry that verdict into Task 5 |

- [ ] **Step 4: Write up in REPORT.md**

Extend the "P1 — Stage A/C alignment puzzle" section with: the izm boundary numbers, the row-shift table, the mask breakdown, the max-diff neighborhood, and a **named root cause** (or, if truly unresolved after all diagnostics, an explicit statement of what was ruled out and the single most informative next measurement). Honesty rule: "unresolved, but ruled out X/Y/Z" is an acceptable outcome; an unsupported root-cause claim is not.

- [ ] **Step 5: Commit**

```bash
git add octave_harness/diag_stage_a.py octave_harness/REPORT.md
git commit -m "diag: stage A row/mask alignment diagnostics (P1b)"
```

---

### Task 5: Rerun the stage diff with corrected methodology (Priority 2)

Apply whatever comparison fix Task 4 identified, make Stage C depth-bin-aware, and produce the definitive first-DIVERGES verdict.

**Files:**
- Modify: `octave_harness/diff_stages.py`
- Modify: `octave_harness/REPORT.md`

**Interfaces:**
- Consumes: Task 4's diagnosis (REPORT.md P1 section); existing `_row_diff`/`_nearest_match` in `diff_stages.py`.
- Produces: updated `diff_stages.py` output incl. a `FIRST DIVERGES:` line; REPORT.md "P2" section.

- [ ] **Step 1: Apply the Stage A comparison fix from Task 4**

This is conditional on Task 4's diagnosis; the allowed edits are **only** to how `diff_stages.py` pairs/filters cells (never to `src/ladcp/` or the dumps). The concrete edit for each branch of the decision table:
- Row/orientation fix: apply the identified row mapping to the Octave arrays before `_row_diff` (e.g. `oct_arr = oct_arr[row_map]`), with a comment citing the diagnostic.
- Masking-policy difference: keep `_row_diff` as is (it already restricts to both-finite) but add a printed `mask_disagree=` percentage per field: `100 * (np.isfinite(o) != np.isfinite(p)).mean()` computed on the matched arrays, so masking differences are visible separately from velocity differences.
- Genuine divergence / unresolved: no edit; the verdict stands and the caveat text in the module docstring is updated to reflect what Task 4 ruled out.

- [ ] **Step 2: Make Stage C depth-bin-aware**

In `diff_stages.py`, Stage C currently nearest-matches super-ensembles by z with no tolerance (`mean_key_err ≈ 1.55 m` noted as non-trivial). Add a tolerance filter inside `_row_diff` via a new optional parameter:

Replace the current opening of `_row_diff` (diff_stages.py:56-61):

```python
def _row_diff(name, oct_field, py_field, oct_key, py_key, key_label,
              rows=None, angular=False, max_key_err=None):
    """Match columns of oct_field/py_field by nearest oct_key<->py_key, diff
    matched rows. max_key_err, if set, drops pairs whose key distance exceeds
    it (depth-bin-aware pairing for Stage C) and reports the drop count.
    Returns a verdict row (name, max|diff|, rms, %match, verdict)."""
    idx, key_err = _nearest_match(oct_key, py_key)
    keep = np.ones(len(idx), dtype=bool)
    if max_key_err is not None:
        keep = key_err <= max_key_err
        print(f"    [{name}] depth-tolerance {max_key_err:.2f}: "
              f"kept {keep.sum()}/{len(keep)} pairs")
    o = np.asarray(oct_field, dtype=float)[..., keep]
    p = np.asarray(py_field, dtype=float)[..., idx[keep]]
    key_err = key_err[keep]
```

All existing call sites pass 2-D fields, so `[..., keep]` selects columns in every case; the rest of the function body (the `rows` slicing, diff, verdict) is unchanged and now operates on `o`/`p` as before.

For the four Stage C calls, compute the tolerance once and pass it:

```python
    dz_se = float(np.nanmedian(np.abs(np.diff(np.sort(oct_z)))))
    print(f"    median super-ensemble depth spacing dz_se = {dz_se:.2f} m")
    # each Stage C call becomes:
    verdicts.append(_row_diff(field, oct_arr, py_arr, oct_z, py_z, "depth",
                              max_key_err=dz_se / 2))
```

Stage A (time-matched) and Stage D calls stay `max_key_err=None` — unfiltered, identical to before.

- [ ] **Step 3: Rerun the stage diff**

Run:
```bash
uv run python octave_harness/diff_stages.py
```
Expected mechanics: runs to completion, prints Stage A (with whatever fix/metric was added), Stage C (with tolerance filter and kept-count), Stage D (should still show ≈ u rms 0.093 — it had no alignment problem; if Stage D *changed*, the `_row_diff` refactor broke something: fix before proceeding), and a final `FIRST DIVERGES:` (or none) line. Record the full output.

- [ ] **Step 4: Verify the refactor didn't disturb the suite, commit**

Run:
```bash
TEST_DATA_DIR=test_data uv run python -m pytest -q
```
Expected: `197 passed, 8 skipped, 2 xfailed`.

- [ ] **Step 5: Write the P2 section in REPORT.md**

Append:

```markdown
## P2 — stage-diff rerun with corrected methodology (session 2026-07-06)

Methodology changes vs 2026-07-05: <Stage A fix applied, or "none needed">;
Stage C now pairs super-ensembles only within dz_se/2 = <value> m
(<n_kept>/<n_total> pairs kept).

<full summary table from the rerun>

**FIRST DIVERGES: <stage/field, or "no stage diverges beyond thresholds">.**
Read together with the P3 solver-only result (<solver exonerated |
implicated>), the Python pipeline's gap against LDEO is attributed to:
<one-paragraph synthesis naming the stage(s)>.
```

- [ ] **Step 6: Commit**

```bash
git add octave_harness/diff_stages.py octave_harness/REPORT.md
git commit -m "feat: depth-bin-aware stage C pairing + stage A fix; first-DIVERGES verdict (P2)"
```

---

### Task 6: Handover — update the pointers for the next session

**Files:**
- Modify: `octave_harness/REPORT.md` (replace the now-executed "What the next session should investigate" section)
- Modify: `HANDOVER.md` (repo root — add a dated pointer)

**Interfaces:**
- Consumes: P3/P1/P2 sections written in Tasks 2-5.
- Produces: the next session's starting point.

- [ ] **Step 1: Rewrite REPORT.md's next-session section**

Replace the entire `## What the next session should investigate` section (currently REPORT.md:226-250, describing P1/P2/P3 — all now executed) with a new section of the same name containing, based on the actual findings:
- If the solver was implicated: the specific getinv.m-vs-compute_inverse difference to chase (constraint weighting, presolve `ut`/`vt` usage from Task 2 Step 4, lsq conditioning), with the measured numbers.
- If the solver was exonerated: the first-DIVERGES stage from Task 5 as the target, plus which upstream port (prepare_superensembles / rotup2down / offsetup2down / editing) to differential-test next, again feeding Octave dumps in as inputs (the same technique as `solver_only.py`, one stage earlier).
- Any P1 residual that stayed unresolved, with the "single most informative next measurement" recorded in Task 4 Step 4.

- [ ] **Step 2: Update HANDOVER.md**

Add at the top of the relevant status section (read the file's existing structure first and match it):

```markdown
**2026-07-06 session:** solver-only harness + stage-diff methodology fixes
executed — see `octave_harness/REPORT.md` sections P3/P1/P2 for the verdict
on where the Python pipeline diverges from LDEO_IX. Headline: <one sentence
with the measured key number>.
```

- [ ] **Step 3: Final full-suite check and commit**

```bash
TEST_DATA_DIR=test_data uv run python -m pytest -q
```
Expected: `197 passed, 8 skipped, 2 xfailed`.

```bash
git add octave_harness/REPORT.md HANDOVER.md
git commit -m "docs: record P3/P1/P2 outcomes, point next session at <the actual finding>"
```
(Fill the commit message's `<the actual finding>` with the real conclusion.)
