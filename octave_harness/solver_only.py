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
