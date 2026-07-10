"""Localize the residual Stage A divergence (post 3-beam, 2026-07-11).

After the depth-registration fix, editing port, sound-speed correction, and
3-beam solutions, Stage A still shows ru/rv rms ~0.085 on both-finite cells
and weight mask_disagree ~48%. No single named mechanism remains, so this
script maps WHERE the residual lives before any further porting:

  1. per-row (bin) breakdown: finite counts (oct/py/both), ru rms, weight
     finite counts, weight rms on both-finite cells;
  2. per-cast-phase (time decile) ru rms;
  3. weight-value scatter summary (are the VALUES on both-finite cells even
     the same quantity?).

Columns are paired by nearest time after undoing Octave's loadctd label
shift (same approach as diff_stages.py Stage A).
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "octave_harness"))

from diag_rmse_strata import run_pipeline  # noqa: E402
from diff_stages import _load, _nearest_match, _oct_label_shift_days  # noqa: E402

DATA_DIR = REPO / "test_data" / "2015_P16N"


def main() -> None:
    stages = {}
    run_pipeline(DATA_DIR, legacy=False, rot=True, offset=True, stages=stages)
    ens = stages["post_edit"]

    d9 = _load(9)["d"]
    oct_time = np.asarray(d9.time_jul, dtype=float) - _oct_label_shift_days()
    idx, terr = _nearest_match(oct_time, ens.time_jul)
    print(f"matched {len(idx)} columns, mean time err {terr.mean()*86400*1000:.1f} ms")

    o_ru = np.asarray(d9.ru, dtype=float)
    p_ru = ens.u[:, idx]
    o_wt = np.asarray(d9.weight, dtype=float)
    p_wt = ens.weight[:, idx]
    o_izm = np.asarray(d9.izm, dtype=float)

    nrow, ncol = o_ru.shape
    print("\n--- per-row (combined-array bin) breakdown ---")
    print(f"{'row':>4} {'oct_fin':>8} {'py_fin':>8} {'both':>8} {'ru_rms':>8} "
          f"{'wt_oct_fin':>10} {'wt_py_fin':>10} {'wt_both':>8} {'wt_rms':>8}")
    for r in range(nrow):
        of = np.isfinite(o_ru[r])
        pf = np.isfinite(p_ru[r])
        both = of & pf
        rms = float(np.sqrt(np.mean((o_ru[r, both] - p_ru[r, both]) ** 2))) \
            if both.sum() else float("nan")
        owf = np.isfinite(o_wt[r])
        pwf = np.isfinite(p_wt[r])
        wboth = owf & pwf
        wrms = float(np.sqrt(np.mean((o_wt[r, wboth] - p_wt[r, wboth]) ** 2))) \
            if wboth.sum() else float("nan")
        print(f"{r:>4} {of.sum():>8} {pf.sum():>8} {both.sum():>8} {rms:>8.4f} "
              f"{owf.sum():>10} {pwf.sum():>10} {wboth.sum():>8} {wrms:>8.4f}")

    print("\n--- per-cast-phase (time decile) ru rms on both-finite cells ---")
    both = np.isfinite(o_ru) & np.isfinite(p_ru)
    col_dec = np.floor(np.arange(ncol) / ncol * 10).astype(int)
    depth_col = np.nanmean(o_izm, axis=0)
    for d in range(10):
        cols = col_dec == d
        b = both[:, cols]
        diff = (o_ru[:, cols] - p_ru[:, cols])[b]
        rms = float(np.sqrt(np.mean(diff ** 2))) if diff.size else float("nan")
        print(f"decile {d}: n={b.sum():>6}  ru_rms={rms:.4f}  "
              f"mean depth {np.nanmean(depth_col[cols]):8.1f} m")

    print("\n--- weight values on both-finite cells ---")
    wboth = np.isfinite(o_wt) & np.isfinite(p_wt)
    ow = o_wt[wboth]
    pw = p_wt[wboth]
    print(f"n={ow.size}")
    print(f"octave weight: mean {ow.mean():.3f}  median {np.median(ow):.3f}  "
          f"p5 {np.percentile(ow,5):.3f}  p95 {np.percentile(ow,95):.3f}")
    print(f"python weight: mean {pw.mean():.3f}  median {np.median(pw):.3f}  "
          f"p5 {np.percentile(pw,5):.3f}  p95 {np.percentile(pw,95):.3f}")
    print(f"corr(oct, py) = {np.corrcoef(ow, pw)[0,1]:.3f}")
    # where does the mask disagreement live?
    oct_only = np.isfinite(o_wt) & ~np.isfinite(p_wt)
    py_only = ~np.isfinite(o_wt) & np.isfinite(p_wt)
    print(f"weight oct-only-finite: {oct_only.mean()*100:.1f}%  "
          f"py-only-finite: {py_only.mean()*100:.1f}%")
    # velocity-finite cells where weight masking disagrees
    vel_both = np.isfinite(o_ru) & np.isfinite(p_ru)
    wm_dis = vel_both & (np.isfinite(o_wt) != np.isfinite(p_wt))
    print(f"cells with both-finite VELOCITY but disagreeing weight mask: "
          f"{wm_dis.sum()} ({100*wm_dis.sum()/max(vel_both.sum(),1):.1f}% of "
          "both-finite velocity cells)")


def ul_shift_scan() -> None:
    """Hypothesis test: is the UL-block residual a UL->DL ensemble pairing
    difference? loadrdi.m merges UL onto DL via bestlag on w (this cast:
    'shift ADCP timeseries by 1 ensembles'); Python uses nearest-time
    (ul_idx). Scan ul_idx + s for s in -2..2 against Octave step09 UL rows:
    if one shift collapses the rms, the pairing is the root cause.
    (Sound-speed scaling ~0.5% is ignored here -- negligible vs the 0.15
    rms effect under test.)
    """
    from ladcp.ingestion.rdi import load_rdi
    from ladcp.transforms.beam2earth import beam2earth, uvrot

    rdi = load_rdi(DATA_DIR / "003DL000.000")
    rdi_ul = load_rdi(DATA_DIR / "003UL000.000")
    u_ul, v_ul, w_ul = beam2earth(
        rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
        rdi_ul.heading, rdi_ul.pitch, rdi_ul.roll,
        20.0, gimbaled=False, beams_up=True, allow_3beam=True,
    )
    u_ul, _ = uvrot(u_ul, v_ul, -12.318441)
    ul_idx = np.argmin(
        np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0
    )

    d9 = _load(9)["d"]
    oct_time = np.asarray(d9.time_jul, dtype=float) - _oct_label_shift_days()
    idx, _ = _nearest_match(oct_time, rdi.time_julian)
    o_ru = np.asarray(d9.ru, dtype=float)

    n_ul = rdi_ul.nbin
    print("\n--- UL column-shift scan (octave UL rows vs python UL bins) ---")
    for s in (-2, -1, 0, 1, 2):
        cols = np.clip(ul_idx + s, 0, rdi_ul.nens - 1)
        # combined row r (0..n_ul-1) = UL bin n_ul-1-r, matched DL columns
        py_ul = u_ul[::-1, :][:, cols][:, idx]  # (n_ul, n_matched)
        o_ul = o_ru[:n_ul, :]
        both = np.isfinite(o_ul) & np.isfinite(py_ul)
        rms = float(np.sqrt(np.mean((o_ul[both] - py_ul[both]) ** 2)))
        print(f"  shift {s:+d}: both-finite n={both.sum():>7}  ru_rms={rms:.4f}")


if __name__ == "__main__":
    main()
    ul_shift_scan()
