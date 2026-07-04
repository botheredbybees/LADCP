"""E2 diagnostic: LDEO's own UL vs DL agreement (VALIDATION_PLAN Phase 1 E2).

test_data/2018_S4P/*.nc embed LDEO_IX per-instrument profiles: u_do/v_do
(downlooker-only) and u_up/v_up (uplooker-only).  If LDEO's transform is
consistent, the two instrument profiles must agree in overlapping depths
WITHOUT any compass-offset correction (the S4P and P16N processing logs show
none was applied).  This establishes the target behaviour for our transform.

Usage:  TEST_DATA_DIR=test_data uv run python scripts/diag_e2_s4p_uldl.py
"""

from __future__ import annotations

import os
from pathlib import Path

import netCDF4
import numpy as np


def main():
    base = Path(os.environ.get("TEST_DATA_DIR", "test_data")) / "2018_S4P"
    for fname in ("001.nc", "002.nc", "003.nc"):
        path = base / fname
        if not path.exists():
            print(f"{fname}: not found, skipping")
            continue
        ds = netCDF4.Dataset(path)
        z = np.asarray(ds.variables["z"][:], dtype=float)
        ud = np.asarray(ds.variables["u_do"][:], dtype=float)
        vd = np.asarray(ds.variables["v_do"][:], dtype=float)
        uu = np.asarray(ds.variables["u_up"][:], dtype=float)
        vu = np.asarray(ds.variables["v_up"][:], dtype=float)
        ds.close()

        ok = np.isfinite(ud) & np.isfinite(uu) & np.isfinite(vd) & np.isfinite(vu)
        n = int(ok.sum())
        if n < 5:
            print(f"{fname}: only {n} overlapping bins, skipping")
            continue
        du = ud[ok] - uu[ok]
        dv = vd[ok] - vu[ok]
        rmse_u = float(np.sqrt(np.mean(du**2)))
        rmse_v = float(np.sqrt(np.mean(dv**2)))
        # rotation angle mapping UL profile onto DL profile
        wd = ud[ok] + 1j * vd[ok]
        wu = uu[ok] + 1j * vu[ok]
        s = np.sum(wd * np.conj(wu))
        rho = np.abs(s) / np.sqrt(np.sum(np.abs(wd) ** 2) * np.sum(np.abs(wu) ** 2))
        theta = np.degrees(np.angle(s))
        print(
            f"{fname}: {n:3d} overlap bins ({z[ok].min():.0f}-{z[ok].max():.0f} m); "
            f"u RMSE(do-up)={rmse_u:.4f}  v RMSE={rmse_v:.4f} m/s; "
            f"rotation fit UL->DL: theta={theta:+.2f} deg (rho={rho:.3f})"
        )


if __name__ == "__main__":
    main()
