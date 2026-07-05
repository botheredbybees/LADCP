"""M3 -- stage-by-stage diff: LDEO_IX process_cast.m (Octave, per-step dumps
in octave_harness/work/dumps/) vs our Python pipeline (scripts/diag_rmse_strata.py's
run_pipeline(), factored to expose intermediate stages).

Config note: the Octave run's own log ("rot up2down use mean up/down compass")
confirms prepinv.m's rotup2down/offsetup2down defaulted ON (matching the
recorded p-struct's rotup2down=1/offsetup2down=1) even though
set_cast_params.m never set them explicitly -- so the fair Python
comparison config is rot=True, offset=True (HANDOVER.md's "NEW + rotup2down
+ offsetup2down" run), not the production default.

Convention notes (see PROGRAMMERS_NOTES.md "Combined DL+UL array layout"):
  - Octave's d.izd/izu are 1-based MATLAB row indices; confirmed this cast:
    izd = rows 26..50, izu = rows 25..1 (reversed) -- i.e. UL-reversed-on-top,
    DL-on-bottom, exactly matching our own EnsembleData/SuperEnsemble layout
    (0-based izd = rows 25..49, izu = rows 24..0 reversed). No remapping
    needed beyond the 1-based/0-based offset.
  - Stage A (post-transform) is matched by nearest ensemble time: Octave's
    step01 dump has already been windowed to the profile time range
    (8679 ensembles, "extracting 8679 ensembles as profile" in the log) vs
    our un-windowed per-file arrays -- same nearest-time-match approach as
    M1's diff_ingestion.py, needed because ping intervals are non-constant.
  - Stages C/D (superensembles, final result) are matched by nearest DEPTH
    (z), the natural key once ensembles have been depth-averaged.
"""
from pathlib import Path

import numpy as np
import scipy.io as sio

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from diag_rmse_strata import run_pipeline  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
DUMPS = REPO / "octave_harness" / "work" / "dumps"
DATA_DIR = REPO / "test_data" / "2015_P16N"


def _load(step: int):
    return sio.loadmat(DUMPS / f"step{step:02d}.mat", struct_as_record=False, squeeze_me=True)


def _nearest_match(key_a, key_b):
    """Return indices into b nearest each element of a (1-D sort-based join)."""
    order = np.argsort(key_b)
    kb_sorted = key_b[order]
    idx = np.searchsorted(kb_sorted, key_a)
    idx = np.clip(idx, 1, len(kb_sorted) - 1)
    left = idx - 1
    use_right = np.abs(kb_sorted[idx] - key_a) < np.abs(kb_sorted[left] - key_a)
    nearest = np.where(use_right, idx, left)
    return order[nearest], np.abs(kb_sorted[nearest] - key_a)


def _row_diff(name, oct_field, py_field, oct_key, py_key, key_label, rows=None, angular=False):
    """Match columns of oct_field/py_field by nearest oct_key<->py_key, diff
    matched rows. Returns a verdict row (name, max|diff|, rms, %match, verdict)."""
    idx, key_err = _nearest_match(oct_key, py_key)
    o = np.asarray(oct_field, dtype=float)
    p = np.asarray(py_field, dtype=float)[..., idx]
    if rows is not None:
        o = o[rows]
        p = p[rows]
    diff = o - p
    if angular:
        diff = (diff + 180.0) % 360.0 - 180.0
    finite = np.isfinite(o) & np.isfinite(p)
    pct_finite_match = 100.0 * finite.sum() / o.size if o.size else float("nan")
    if finite.sum() == 0:
        return (name, float("nan"), float("nan"), pct_finite_match, "NO-OVERLAP")
    max_diff = float(np.nanmax(np.abs(diff[finite])))
    rms_diff = float(np.sqrt(np.nanmean(diff[finite] ** 2)))
    mean_key_err = float(np.mean(key_err))
    if max_diff < 1e-6:
        verdict = "MATCH"
    elif rms_diff < 0.02:
        verdict = "NEAR"
    else:
        verdict = "DIVERGES"
    print(
        f"{name:<28}{key_label:<10}max|diff|={max_diff:>10.4g}  rms={rms_diff:>10.4g}  "
        f"%finite-match={pct_finite_match:5.1f}  mean_key_err={mean_key_err:.3g}  [{verdict}]"
    )
    return (name, max_diff, rms_diff, pct_finite_match, verdict)


def main() -> None:
    print("=" * 100)
    print("M3 stage diff: LDEO_IX process_cast.m (Octave) vs Python pipeline")
    print("=" * 100)

    stages = {}
    res, ref = run_pipeline(DATA_DIR, legacy=False, rot=True, offset=True, stages=stages)

    verdicts = []

    # --- Stage A: post-EDIT (Octave step09 vs our post_edit ens) ---
    # NOT step01/post_transform: Octave's step01 dump already has loadrdi's
    # internal outlier() masking applied ("Outlier discarded 6953 bins down
    # looking" etc in the log), while Python's post_transform is captured
    # before any editing -- comparing those two is a different-pipeline-point
    # mismatch, not a real divergence (caught via a >>ocean-velocity 14 m/s
    # "diff" and 37% finite-overlap on the first attempt). step09/post_edit
    # is the first point both pipelines have completed all masking.
    print("\n--- Stage A: post-edit (masking + sidelobe + large-vel + w-outlier QA) ---")
    step9 = _load(9)
    d9 = step9["d"]
    ens_pe = stages["post_edit"]
    oct_time = np.asarray(d9.time_jul, dtype=float)
    py_time = ens_pe.time_jul
    for field, oct_arr, py_arr in [
        ("ru (east vel, all bins)", d9.ru, ens_pe.u),
        ("rv (north vel, all bins)", d9.rv, ens_pe.v),
        ("rw (vert vel, all bins)", d9.rw, ens_pe.w),
        ("weight (post-edit)", d9.weight, ens_pe.weight),
    ]:
        verdicts.append(_row_diff(field, oct_arr, py_arr, oct_time, py_time, "time"))

    # --- Stage C: super-ensembles (Octave step12 di vs our se) ---
    # CAVEAT: matched by nearest depth (z) -- superensemble bin centers don't
    # necessarily coincide between the two pipelines (different averaging
    # windows/edges), so a nonzero mean_key_err here reflects bin-grid
    # misalignment as well as genuine data differences. Not a depth-bin-aware
    # comparison; a DIVERGES verdict here needs that caveat, not a bare claim.
    print("\n--- Stage C: super-ensembles (the solver's actual input) ---")
    print("    CAVEAT: nearest-depth matching only -- see module docstring/comment above")
    step12 = _load(12)
    di = step12["di"]
    se = stages["superensembles"]
    oct_z = np.asarray(di.z, dtype=float)
    py_z = se.z
    for field, oct_arr, py_arr in [
        ("ru (super-ens east vel)", di.ru, se.ru),
        ("rv (super-ens north vel)", di.rv, se.rv),
        ("rw (super-ens vert vel)", di.rw, se.rw),
        ("weight (super-ens)", di.weight, se.weight),
    ]:
        verdicts.append(_row_diff(field, oct_arr, py_arr, oct_z, py_z, "depth"))

    # --- Stage D: final result (Octave step17 dr vs our InverseResult) ---
    print("\n--- Stage D: final inverse solution (u/v profile) ---")
    step17 = _load(17)
    dr = step17["dr"]
    oct_z_final = np.asarray(dr.z, dtype=float)
    py_z_final = res.z
    for field, oct_arr, py_arr in [
        ("u (final east vel)", dr.u, res.u),
        ("v (final north vel)", dr.v, res.v),
    ]:
        verdicts.append(_row_diff(field, oct_arr, py_arr, oct_z_final, py_z_final, "depth"))

    print("\n--- Summary ---")
    print(f"{'stage/field':<28}{'max|diff|':>12}{'rms':>12}{'%match':>10}{'verdict':>12}")
    first_diverge = None
    for name, max_diff, rms_diff, pct, verdict in verdicts:
        print(f"{name:<28}{max_diff:>12.4g}{rms_diff:>12.4g}{pct:>9.1f}%{verdict:>12}")
        if verdict == "DIVERGES" and first_diverge is None:
            first_diverge = name
    if first_diverge:
        print(f"\nFIRST DIVERGES: {first_diverge}")
    else:
        print("\nNo stage DIVERGES (within the thresholds used here).")


if __name__ == "__main__":
    main()
