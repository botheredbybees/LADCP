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
  - Stage A (post-edit, Octave step09 vs Python's post_edit stage -- see the
    inline comment at the Stage A block below for why step09 and not step01)
    is matched by nearest ensemble time: Octave's dumps have already been
    windowed to the profile time range vs our un-windowed per-file arrays --
    same nearest-time-match approach as M1's diff_ingestion.py, needed
    because ping intervals are non-constant.
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


def _oct_label_shift_days():
    """Measure the ADCP time-label shift loadctd.m applied in the Octave run.

    loadctd.m line 443 shifts d.time_jul by its besttlag lagdt (-23.24 s this
    cast), moving Octave's labels off the ADCP clock that Python's labels use.
    Worse, this cast's staggered ping pattern (1.33/1.58 s) repeats every
    2.91 s and 23.24 s = 16 pings almost exactly, so nearest-TIME matching
    silently pairs ensembles 16 pings apart with only ~13 ms apparent error
    (P1's "column/time mismatch ruled out" was fooled by this). Content-match
    the heading series (untouched by loadctd) between the pre-loadctd step03
    dump and step09 to recover the true pairing, and return the label shift
    to subtract from Octave times before time-matching against Python.
    """
    d3 = _load(3)["d"]
    d9 = _load(9)["d"]
    h3 = np.asarray(d3.hdg, dtype=float)
    h9 = np.asarray(d9.hdg, dtype=float)
    if h3.ndim > 1:
        h3 = h3[0]
    if h9.ndim > 1:
        h9 = h9[0]
    c3 = np.exp(1j * np.radians(h3))
    c9 = np.exp(1j * np.radians(h9))
    best_k, best_co = 0, -1.0
    for k in range(h3.size - h9.size + 1):
        co = float(np.abs(np.vdot(c3[k:k + h9.size], c9))) / h9.size
        if co > best_co:
            best_k, best_co = k, co
    t3 = np.asarray(d3.time_jul, dtype=float)
    t9 = np.asarray(d9.time_jul, dtype=float)
    shift_days = float(np.mean(t9 - t3[best_k:best_k + t9.size]))
    print(f"    Octave label shift (loadctd lagdt): {shift_days*86400:+.2f} s "
          f"(heading content-match k={best_k}, coherence={best_co:.4f})")
    return shift_days


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


def _row_diff(name, oct_field, py_field, oct_key, py_key, key_label,
              rows=None, angular=False, max_key_err=None):
    """Match columns of oct_field/py_field by nearest oct_key<->py_key, diff
    matched rows. max_key_err, if set, drops pairs whose key distance exceeds
    it (depth-bin-aware pairing for Stage C) and reports the drop count.
    Also prints mask_disagree -- the % of matched cells where the two
    pipelines disagree on finite-vs-masked (isfinite(o) != isfinite(p)) --
    separately from the rms, which is restricted to both-finite cells only,
    so a masking-policy difference doesn't get conflated with a velocity
    difference (Task 4's P1 finding: the Stage A 14 m/s max|diff| was a
    Python-side unmasked outlier, not a real divergence).
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
    if rows is not None:
        o = o[rows]
        p = p[rows]
    diff = o - p
    if angular:
        diff = (diff + 180.0) % 360.0 - 180.0
    finite = np.isfinite(o) & np.isfinite(p)
    pct_finite_match = 100.0 * finite.sum() / o.size if o.size else float("nan")
    mask_disagree = (
        100.0 * (np.isfinite(o) != np.isfinite(p)).mean() if o.size else float("nan")
    )
    if finite.sum() == 0:
        return (name, float("nan"), float("nan"), pct_finite_match, "NO-OVERLAP")
    max_diff = float(np.nanmax(np.abs(diff[finite])))
    rms_diff = float(np.sqrt(np.nanmean(diff[finite] ** 2)))
    mean_key_err = float(np.mean(key_err)) if key_err.size else float("nan")
    if max_diff < 1e-6:
        verdict = "MATCH"
    elif rms_diff < 0.02:
        verdict = "NEAR"
    else:
        verdict = "DIVERGES"
    print(
        f"{name:<28}{key_label:<10}max|diff|={max_diff:>10.4g}  rms={rms_diff:>10.4g}  "
        f"%finite-match={pct_finite_match:5.1f}  mean_key_err={mean_key_err:.3g}  "
        f"mask_disagree={mask_disagree:5.1f}%  [{verdict}]"
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
    # Undo loadctd.m's ADCP time-label shift so nearest-time matching pairs
    # the SAME physical pings (see _oct_label_shift_days docstring).
    oct_time = np.asarray(d9.time_jul, dtype=float) - _oct_label_shift_days()
    py_time = ens_pe.time_jul
    for field, oct_arr, py_arr in [
        ("ru (east vel, all bins)", d9.ru, ens_pe.u),
        ("rv (north vel, all bins)", d9.rv, ens_pe.v),
        ("rw (vert vel, all bins)", d9.rw, ens_pe.w),
        ("weight (post-edit)", d9.weight, ens_pe.weight),
    ]:
        verdicts.append(_row_diff(field, oct_arr, py_arr, oct_time, py_time, "time"))

    # Task 4's (P1) named next measurement: full-array d.izm - ens.izm
    # statistics, to quantify the depth-registration offset it found at two
    # sample columns (an early ~15 m estimate). P2's full-array measurement
    # superseded that: the offset is depth-varying (rms 47.9 m, max|diff|
    # 108.8 m), growing in magnitude with cast depth rather than being a
    # fixed ~15 m shift. izm is METRES OF DEPTH, not m/s -- a large rms
    # here is the EXPECTED finding (confirming/quantifying the offset), so
    # a DIVERGES verdict on this row is correct and informative, not a
    # velocity divergence. Kept out of the velocity-only FIRST-DIVERGES
    # tally in the Summary below.
    print("    NOTE: izm units are metres of depth, not m/s -- DIVERGES here "
          "means the depth-varying registration offset (P2: rms 47.9 m, "
          "growing in magnitude with cast depth), not a velocity gap.")
    izm_verdict = _row_diff(
        "izm (depth registration, all bins)", d9.izm, ens_pe.izm,
        oct_time, py_time, "time",
    )
    verdicts.append(izm_verdict)

    # Reproducible summary stats for the izm depth-registration offset (P2
    # measurement) -- computed here so REPORT.md's quoted mean/median/corr
    # numbers can be regenerated from this tracked script instead of the
    # ad-hoc scratch code used to first derive them.
    izm_idx, _ = _nearest_match(oct_time, py_time)
    izm_o = np.asarray(d9.izm, dtype=float)
    izm_p = np.asarray(ens_pe.izm, dtype=float)[..., izm_idx]
    izm_finite = np.isfinite(izm_o) & np.isfinite(izm_p)
    izm_diff = izm_o[izm_finite] - izm_p[izm_finite]
    izm_depth = izm_o[izm_finite]
    izm_corr = float(np.corrcoef(izm_diff, izm_depth)[0, 1])
    print(
        f"    izm diff (oct - py): mean={izm_diff.mean():.2f} "
        f"median={np.median(izm_diff):.2f} min={izm_diff.min():.2f} "
        f"max={izm_diff.max():.2f}  corr(diff, depth)={izm_corr:.3f}  "
        "(both negative-down: negative diff = Octave deeper)"
    )

    # --- Stage C: super-ensembles (Octave step12 di vs our se) ---
    # Depth-bin-aware pairing (Task 5, P2): superensemble bin centers don't
    # necessarily coincide between the two pipelines (different averaging
    # windows/edges), so unfiltered nearest-depth matching (mean_key_err was
    # ~1.55 m -- non-trivial) mixed bin-grid misalignment into the diff.
    # Pairs whose depth distance exceeds half the median super-ensemble
    # depth spacing are now dropped instead of force-matched.
    print("\n--- Stage C: super-ensembles (the solver's actual input) ---")
    step12 = _load(12)
    di = step12["di"]
    se = stages["superensembles"]
    oct_z = np.asarray(di.z, dtype=float)
    py_z = se.z
    dz_se = float(np.nanmedian(np.abs(np.diff(np.sort(oct_z)))))
    print(f"    median super-ensemble depth spacing dz_se = {dz_se:.2f} m")
    for field, oct_arr, py_arr in [
        ("ru (super-ens east vel)", di.ru, se.ru),
        ("rv (super-ens north vel)", di.rv, se.rv),
        ("rw (super-ens vert vel)", di.rw, se.rw),
        ("weight (super-ens)", di.weight, se.weight),
    ]:
        verdicts.append(_row_diff(field, oct_arr, py_arr, oct_z, py_z, "depth",
                                  max_key_err=dz_se / 2))

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
    first_diverge_velocity = None
    for name, max_diff, rms_diff, pct, verdict in verdicts:
        print(f"{name:<28}{max_diff:>12.4g}{rms_diff:>12.4g}{pct:>9.1f}%{verdict:>12}")
        if verdict == "DIVERGES" and first_diverge is None:
            first_diverge = name
        # izm is metres of depth, not m/s -- its expected depth-varying
        # DIVERGES (P2's depth-registration finding: rms 47.9 m, growing
        # in magnitude with cast depth) must not stand in for a velocity
        # verdict, so it's tracked separately here.
        if verdict == "DIVERGES" and first_diverge_velocity is None and "izm" not in name:
            first_diverge_velocity = name
    if first_diverge:
        print(f"\nFIRST DIVERGES: {first_diverge}")
    else:
        print("\nNo stage DIVERGES (within the thresholds used here).")
    print(
        "FIRST DIVERGES (velocity fields only, excl. izm depth-registration): "
        f"{first_diverge_velocity if first_diverge_velocity else 'none'}"
    )


if __name__ == "__main__":
    main()
