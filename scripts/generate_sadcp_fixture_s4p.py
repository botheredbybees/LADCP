"""Generate SADCP fixture for S4P casts 001, 002, 003 integration tests.

Reads os75nb_short.nc, averages water velocity over each cast time window
(extracted from the reference NC outputs), and writes sadcp_NNN.npz files.

Usage:
    python scripts/generate_sadcp_fixture_s4p.py

Requires TEST_DATA_DIR env var pointing to a directory containing 2018_S4P/.
"""
import os
import sys
from pathlib import Path

import netCDF4
import numpy as np

# Cast metadata: (cast_label, approx_lat, approx_lon)
CASTS = [
    ("001", -70.45, 168.47),
    ("002", -70.36, 168.63),
    ("003", -70.10, 169.13),
]


def main() -> None:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        sys.exit("ERROR: TEST_DATA_DIR env var not set")

    base = Path(env) / "2018_S4P"
    sadcp_nc = base / "SADCP/os75nb/contour/os75nb_short.nc"

    if not sadcp_nc.exists():
        sys.exit(f"ERROR: SADCP file not found: {sadcp_nc}")

    from ladcp.ingestion.sadcp import load_sadcp_nc

    for cast, lat, lon in CASTS:
        ref_nc = base / f"{cast}.nc"
        if not ref_nc.exists():
            print(f"SKIP {cast}: reference {ref_nc} not found")
            continue

        # Extract cast time window from reference output (ensemble Julian times)
        ds = netCDF4.Dataset(str(ref_nc))
        tim = np.asarray(ds.variables["tim"][:], dtype=float)
        ds.close()
        t_start = float(tim[0])
        t_end = float(tim[-1])

        profile = load_sadcp_nc(
            sadcp_nc, t_start, t_end,
            lat=lat, lon=lon,
        )
        if profile is None:
            print(f"WARN {cast}: no SADCP records in time window")
            continue

        out = base / f"sadcp_{cast}.npz"
        np.savez(str(out), z=profile.z, u=profile.u, v=profile.v, err=profile.err)
        print(f"Written {out}")
        print(f"  {len(profile.z)} bins, {profile.z.min():.0f}–{profile.z.max():.0f} m")
        print(f"  u={profile.u.mean():.4f} v={profile.v.mean():.4f} m/s")


if __name__ == "__main__":
    main()
