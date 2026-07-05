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
