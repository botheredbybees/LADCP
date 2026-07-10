"""Formation-only differential harness (REPORT.md P5 handoff item 2).

Feeds the exact input LDEO_IX's prepinv.m received -- the `d` struct dumped
after step09 (EDIT DATA) -- into Python's rotup2down() +
prepare_superensembles(), and diffs the result against Octave's own step10
`di` (FORM SUPER ENSEMBLES). Both implementations see IDENTICAL input, so
any difference isolates the super-ensemble formation logic -- the last
stage that diverges (Stage A is NEAR since P5, and P3 exonerated the
solver).

Run from the LADCP repo root:
    uv run python octave_harness/formation_only.py
"""
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from ladcp.solution.inverse import (  # noqa: E402
    EnsembleData,
    prepare_superensembles,
    rotup2down,
)

DUMPS = REPO / "octave_harness" / "work" / "dumps"


def _load(step: int):
    return sio.loadmat(
        DUMPS / f"step{step:02d}.mat", struct_as_record=False, squeeze_me=True
    )


def _f(x):
    return np.asarray(x, dtype=float)


def build_ens(d) -> tuple[EnsembleData, np.ndarray, np.ndarray]:
    """Map the Octave post-edit d struct onto EnsembleData."""
    bvel = np.atleast_2d(_f(d.bvel))
    if bvel.shape[0] in (3, 4):
        bvel = bvel.T                      # -> (nens, 3|4)
    ens = EnsembleData(
        u=_f(d.ru), v=_f(d.rv), w=_f(d.rw),
        weight=_f(d.weight), izm=_f(d.izm),
        z=_f(d.z), time_jul=_f(d.time_jul),
        bvel=bvel[:, :3], bvels=np.full_like(bvel[:, :3], 0.02),
        hbot=_f(d.hbot),
        izd=np.asarray(d.izd, dtype=int).ravel() - 1,
        izu=np.asarray(d.izu, dtype=int).ravel() - 1,
        slat=_f(d.slat), slon=_f(d.slon),
    )
    hdg = _f(d.hdg)
    return ens, hdg[0], hdg[1]


def diff_field(name, o, p, both_shape=True):
    o = _f(o)
    p = _f(p)
    if o.shape != p.shape:
        print(f"{name:<10} SHAPE MISMATCH oct{o.shape} vs py{p.shape}")
        return
    both = np.isfinite(o) & np.isfinite(p)
    mask_dis = 100.0 * (np.isfinite(o) != np.isfinite(p)).mean()
    if both.sum() == 0:
        print(f"{name:<10} no overlap")
        return
    diff = o[both] - p[both]
    print(f"{name:<10} rms={np.sqrt(np.mean(diff**2)):>10.4g}  "
          f"max|diff|={np.abs(diff).max():>10.4g}  "
          f"n_both={both.sum():>7}  mask_disagree={mask_dis:5.1f}%")


def main() -> None:
    d9 = _load(9)["d"]
    di = _load(10)["di"]

    ens, hdg_dl, hdg_ul = build_ens(d9)
    ens_rot, off = rotup2down(ens, hdg_dl, hdg_ul)
    print(f"rotup2down offset: {off if np.isscalar(off) else 'series'}")
    # superens_std_min = Single_Ping_Err/sqrt(Pings_per_Ensemble) for the
    # P16N WH300s -- quoted by the Octave run's own log ("set 1360 values to
    # minimum super ensemble std 0.083833"). outlier_nblock = LDEO's
    # p.outlier_n, set at loadrdi from the RAW ping rate and reused verbatim
    # by prepinv's outlier(di, p) call.
    import math
    dt_min = float(np.nanmean(np.diff(_f(d9.time_jul)))) * 24.0 * 60.0
    nblock = int(math.ceil(5.0 / dt_min))
    print(f"outlier_nblock (raw ping rate): {nblock}")
    se = prepare_superensembles(
        ens_rot, superens_std_min=0.083833,
        outlier_nblock=nblock, tilt_deg=_f(d9.tilt),
    )

    oct_z = _f(di.z)
    print(f"\nn_se: octave {oct_z.size}  python {se.z.size}")
    if oct_z.size == se.z.size:
        print("index-paired comparison:")
        diff_field("z", oct_z, se.z)
        diff_field("dt", di.dt, se.dt)
        diff_field("izm", di.izm, se.izm)
        diff_field("ru", di.ru, se.ru)
        diff_field("rv", di.rv, se.rv)
        diff_field("rw", di.rw, se.rw)
        diff_field("ruvs", di.ruvs, se.ruvs)
        diff_field("weight", di.weight, se.weight)
    else:
        # nearest-depth pairing fallback
        order = np.argsort(se.z)
        zs = se.z[order]
        pos = np.clip(np.searchsorted(zs, oct_z), 1, len(zs) - 1)
        nearer = np.where(
            np.abs(zs[pos] - oct_z) < np.abs(zs[pos - 1] - oct_z), pos, pos - 1
        )
        idx = order[nearer]
        zerr = np.abs(se.z[idx] - oct_z)
        print(f"nearest-z pairing: mean|dz|={zerr.mean():.2f} m  "
              f"median={np.median(zerr):.2f}  max={zerr.max():.2f}")
        diff_field("z", oct_z, se.z[idx])
        diff_field("ru", di.ru, se.ru[:, idx])
        diff_field("rv", di.rv, se.rv[:, idx])
        diff_field("rw", di.rw, se.rw[:, idx])
        diff_field("ruvs", di.ruvs, se.ruvs[:, idx])
        diff_field("weight", di.weight, se.weight[:, idx])


if __name__ == "__main__":
    main()
