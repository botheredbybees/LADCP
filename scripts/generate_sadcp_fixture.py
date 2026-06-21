"""Generate SADCP fixture for the P16N cast 003 integration test.

Reads the real OS75 shipboard ADCP NetCDF, averages water velocity over the
cast 003 time window, and writes sadcp_003.npz with independent SADCP data
covering ~29–733 m depth.  The result provides an independent (non-circular)
constraint for the inverse solver integration test.

Usage:
    python scripts/generate_sadcp_fixture.py

Requires TEST_DATA_DIR env var pointing to the directory containing
2015_P16N/. Writes 2015_P16N/sadcp_003.npz.

Note: the previous version of this script used the 003.nc LDEO_IX reference
output as synthetic SADCP.  That was circular (the answer used as constraint)
and has been replaced with this approach using the real OS75 ADCP file.
"""
import os
import sys
from pathlib import Path

import netCDF4
import numpy as np

# P16N cast 003 time window in Julian days (2015-01-01 epoch)
# Derived from 003DL000.000 ensemble time range; ±0.01 JD buffer applied.
CAST_T_START_JD = 2457124.7336
CAST_T_END_JD = 2457124.8817
JAN1_2015_JD = 2457023.5  # Julian day of 2015-01-01 00:00 UTC

# SADCP measurement uncertainty (m/s) — conservative estimate for OS75
SADCP_ERR = 0.05


def _jan1_jd() -> float:
    """Julian day of 2015-01-01 00:00 UTC, computed from the ingestion helper."""
    from ladcp.ingestion._pd0 import _to_julian
    return _to_julian(2015, 1, 1, 0.0)


def main() -> None:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        sys.exit("ERROR: TEST_DATA_DIR env var not set")
    base = Path(env) / "2015_P16N"
    src = base / "2015_P16N_leg1_os75nb_short.nc"
    if not src.exists():
        sys.exit(f"ERROR: OS75 SADCP file not found: {src}")

    jan1_jd = _jan1_jd()

    ds = netCDF4.Dataset(src)
    # time is in "days since 2015-01-01 00:00:00" → add epoch to get JD
    t_jd = np.asarray(ds.variables["time"][:], dtype=float) + jan1_jd
    # u, v are water-column velocity (Earth frame, m/s); masked arrays
    u_raw = np.ma.filled(ds.variables["u"][:], np.nan)   # (ntime, ndepth)
    v_raw = np.ma.filled(ds.variables["v"][:], np.nan)
    depth_raw = np.ma.filled(ds.variables["depth"][:], np.nan)  # (ntime, ndepth)
    ds.close()

    # Select time records that overlap the cast window
    in_window = (t_jd >= CAST_T_START_JD - 0.01) & (t_jd <= CAST_T_END_JD + 0.01)
    n_records = int(in_window.sum())
    if n_records == 0:
        sys.exit("ERROR: No SADCP records found in cast 003 time window")
    print(f"Records in window: {n_records}")

    # Average water velocity over the window per depth bin
    u_avg = np.nanmean(u_raw[in_window], axis=0)   # (ndepth,)
    v_avg = np.nanmean(v_raw[in_window], axis=0)
    depth_avg = np.nanmean(depth_raw[in_window], axis=0)

    # Keep only bins where all three quantities are finite
    valid = np.isfinite(u_avg) & np.isfinite(v_avg) & np.isfinite(depth_avg)
    z_sel = depth_avg[valid]
    u_sel = u_avg[valid]
    v_sel = v_avg[valid]
    err_sel = np.full_like(u_sel, SADCP_ERR)

    out = base / "sadcp_003.npz"
    np.savez(out, z=z_sel, u=u_sel, v=v_sel, err=err_sel)
    print(f"Written: {out}")
    print(f"  {valid.sum()} depth bins, depth range {z_sel.min():.0f}–{z_sel.max():.0f} m")
    print(f"  u mean={u_sel.mean():.4f}, v mean={v_sel.mean():.4f} m/s")
    print(f"  u std={u_sel.std():.4f}, v std={v_sel.std():.4f} m/s")


if __name__ == "__main__":
    main()
