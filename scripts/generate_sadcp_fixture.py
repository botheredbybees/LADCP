"""Generate synthetic SADCP fixture from LDEO_IX reference NetCDF.

Usage:
    python scripts/generate_sadcp_fixture.py

Requires TEST_DATA_DIR env var pointing to the directory containing
2015_P16N/003.nc. Writes 2015_P16N/sadcp_003.npz.
"""
import os
import sys
from pathlib import Path

import netCDF4
import numpy as np


def main() -> None:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        sys.exit("ERROR: TEST_DATA_DIR env var not set")
    base = Path(env) / "2015_P16N"
    src = base / "003.nc"
    if not src.exists():
        sys.exit(f"ERROR: reference file not found: {src}")

    ds = netCDF4.Dataset(src)
    z = np.array(ds.variables["z"][:], dtype=np.float64)      # positive m
    u = np.array(ds.variables["u"][:], dtype=np.float64)      # m/s
    v = np.array(ds.variables["v"][:], dtype=np.float64)      # m/s
    nvel = np.array(ds.variables["nvel"][:], dtype=np.int32)
    ds.close()

    sel = np.isfinite(u) & np.isfinite(v) & (nvel >= 3)
    z_sel = z[sel]
    u_sel = u[sel]
    v_sel = v[sel]
    err_sel = np.full_like(u_sel, 0.05)   # matches InverseParams.velerr default

    out = base / "sadcp_003.npz"
    np.savez(out, z=z_sel, u=u_sel, v=v_sel, err=err_sel)
    print(f"Written: {out}")
    print(f"  {sel.sum()} depth bins, depth range {z_sel.min():.0f}–{z_sel.max():.0f} m")


if __name__ == "__main__":
    main()
