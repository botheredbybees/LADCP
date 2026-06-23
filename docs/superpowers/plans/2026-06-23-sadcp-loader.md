# SADCP Loader Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `load_sadcp_nc()` in `src/ladcp/ingestion/sadcp.py` that reads OS75 shipboard ADCP NetCDF files (JASADCP format), time-windows records to a cast, averages to a single profile, and returns arrays ready for `compute_inverse()`.

**Architecture:** The existing `scripts/generate_sadcp_fixture.py` already does the core operation for P16N — this plan generalises it into a proper library function. A `SADCPProfile` dataclass carries `z/u/v/err`. The loader handles the JASADCP longitude offset convention (`lon % 360`), optional position sanity-check, and per-bin error scaling that mirrors `loadsadcp.m` (scale by `max(nvel)/nvel` so bins with fewer valid pings get larger error). A fixture-generator script for S4P and an integration test against `001.nc` close the loop.

**Tech Stack:** Python 3.11, netCDF4, numpy, pytest

## Global Constraints

- Package name: `ladcp`
- Test env var: `TEST_DATA_DIR` pointing to a directory that contains `2018_S4P/` and `2015_P16N/`
- Integration tests carry `@pytest.mark.integration` and `scope="module"` fixtures
- `pytest.skip()` (not `pytest.xfail`) when test data files are absent
- `err_default = 0.05` m/s matches `SADCP_ERR` in the existing P16N fixture script
- Do not modify `scripts/generate_sadcp_fixture.py` (P16N fixture still works)

---

### Task 1: `SADCPProfile` dataclass and `load_sadcp_nc()` function

**Files:**
- Create: `src/ladcp/ingestion/sadcp.py`
- Create: `tests/unit/test_sadcp_loader.py`

**Interfaces:**
- Produces: `SADCPProfile(z, u, v, err)` and `load_sadcp_nc(path, t_start_jd, t_end_jd, *, lat, lon, pos_tol, err_default, dt_slack_days, max_depth) -> SADCPProfile | None`
- Consumed by Task 2 (unit tests) and Task 3 (S4P fixture generator)

- [ ] **Step 1: Write failing unit tests**

Create `tests/unit/test_sadcp_loader.py`:

```python
"""Unit tests for load_sadcp_nc()."""
from __future__ import annotations
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.sadcp import SADCPProfile, load_sadcp_nc

# Julian day of 2020-01-01 00:00 UTC (well-known reference date)
JAN1_2020_JD = 2458849.5


@pytest.fixture
def sadcp_nc(tmp_path: Path) -> Path:
    """Minimal SADCP NetCDF with 10 time steps, 5 depth bins."""
    nc_path = tmp_path / "test_sadcp.nc"
    ds = netCDF4.Dataset(str(nc_path), "w")
    ds.createDimension("time", 10)
    ds.createDimension("depth_cell", 5)

    t = ds.createVariable("time", "f8", ("time",))
    t.units = "days since 2020-01-01 00:00:00"
    # Indices 2-7 fall in window [1.0, 1.5]; rest outside
    t[:] = [0.0, 0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 2.5]

    lon_v = ds.createVariable("lon", "f4", ("time",))
    lon_v[:] = [170.0] * 10

    lat_v = ds.createVariable("lat", "f4", ("time",))
    lat_v[:] = [-45.0] * 10

    depth_v = ds.createVariable("depth", "f4", ("depth_cell",))
    depth_v[:] = [50.0, 100.0, 150.0, 200.0, 250.0]

    u_v = ds.createVariable("u", "f4", ("time", "depth_cell"), fill_value=1e35)
    u_v[:] = 0.2  # constant for easy assertion

    v_v = ds.createVariable("v", "f4", ("time", "depth_cell"), fill_value=1e35)
    v_v[:] = 0.1

    ds.close()
    return nc_path


def test_load_returns_profile(sadcp_nc: Path) -> None:
    profile = load_sadcp_nc(sadcp_nc, JAN1_2020_JD + 1.0, JAN1_2020_JD + 1.5)
    assert isinstance(profile, SADCPProfile)
    assert len(profile.z) == 5
    np.testing.assert_allclose(profile.u, 0.2, atol=1e-4)
    np.testing.assert_allclose(profile.v, 0.1, atol=1e-4)
    assert np.all(profile.err > 0)


def test_load_returns_none_outside_window(sadcp_nc: Path) -> None:
    profile = load_sadcp_nc(sadcp_nc, JAN1_2020_JD + 5.0, JAN1_2020_JD + 6.0)
    assert profile is None


def test_position_check_raises_on_mismatch(sadcp_nc: Path) -> None:
    with pytest.raises(ValueError, match="position"):
        load_sadcp_nc(
            sadcp_nc, JAN1_2020_JD + 1.0, JAN1_2020_JD + 1.5,
            lat=-30.0, lon=170.0,  # wrong latitude
        )


def test_position_check_passes_when_close(sadcp_nc: Path) -> None:
    profile = load_sadcp_nc(
        sadcp_nc, JAN1_2020_JD + 1.0, JAN1_2020_JD + 1.5,
        lat=-45.0, lon=170.0,
    )
    assert profile is not None


def test_lon_normalization(tmp_path: Path) -> None:
    """lon offset by -360 (e.g. 168°E stored as -192) normalises correctly."""
    nc_path = tmp_path / "test_lon.nc"
    ds = netCDF4.Dataset(str(nc_path), "w")
    ds.createDimension("time", 3)
    ds.createDimension("depth_cell", 2)
    t = ds.createVariable("time", "f8", ("time",))
    t.units = "days since 2020-01-01 00:00:00"
    t[:] = [1.0, 1.2, 1.4]
    lon_v = ds.createVariable("lon", "f4", ("time",))
    lon_v[:] = [-192.0, -191.9, -191.8]  # 168.x°E stored offset by -360
    lat_v = ds.createVariable("lat", "f4", ("time",))
    lat_v[:] = [-70.4] * 3
    d = ds.createVariable("depth", "f4", ("depth_cell",))
    d[:] = [50.0, 100.0]
    u_v = ds.createVariable("u", "f4", ("time", "depth_cell"), fill_value=1e35)
    u_v[:] = 0.1
    v_v = ds.createVariable("v", "f4", ("time", "depth_cell"), fill_value=1e35)
    v_v[:] = 0.05
    ds.close()

    profile = load_sadcp_nc(
        nc_path, JAN1_2020_JD + 0.9, JAN1_2020_JD + 1.5,
        lat=-70.4, lon=168.1,
    )
    assert profile is not None


def test_max_depth_filter(sadcp_nc: Path) -> None:
    profile = load_sadcp_nc(
        sadcp_nc, JAN1_2020_JD + 1.0, JAN1_2020_JD + 1.5,
        max_depth=150.0,  # should exclude 200 and 250 m bins
    )
    assert profile is not None
    assert float(profile.z.max()) <= 150.0


def test_error_scales_with_nvel(tmp_path: Path) -> None:
    """A bin with fewer valid records gets proportionally larger error."""
    nc_path = tmp_path / "test_nvel.nc"
    ds = netCDF4.Dataset(str(nc_path), "w")
    ds.createDimension("time", 4)
    ds.createDimension("depth_cell", 2)
    t = ds.createVariable("time", "f8", ("time",))
    t.units = "days since 2020-01-01 00:00:00"
    t[:] = [1.0, 1.1, 1.2, 1.3]
    lon_v = ds.createVariable("lon", "f4", ("time",))
    lon_v[:] = [170.0] * 4
    lat_v = ds.createVariable("lat", "f4", ("time",))
    lat_v[:] = [-45.0] * 4
    d = ds.createVariable("depth", "f4", ("depth_cell",))
    d[:] = [50.0, 100.0]
    # Bin 0: all 4 records valid; Bin 1: only 2 records valid (rest NaN)
    u_v = ds.createVariable("u", "f4", ("time", "depth_cell"), fill_value=1e35)
    u_data = np.array([[0.2, 0.3], [0.2, np.nan], [0.2, 0.3], [0.2, np.nan]])
    u_v[:] = u_data
    v_v = ds.createVariable("v", "f4", ("time", "depth_cell"), fill_value=1e35)
    v_v[:] = 0.1
    ds.close()

    profile = load_sadcp_nc(nc_path, JAN1_2020_JD + 0.9, JAN1_2020_JD + 1.4)
    assert profile is not None
    # Bin 1 (depth 100m) has half the observations → error should be larger
    assert profile.err[1] > profile.err[0]
```

- [ ] **Step 2: Run to verify failure**

```
pytest tests/unit/test_sadcp_loader.py -v
```

Expected: `ImportError: cannot import name 'SADCPProfile' from 'ladcp.ingestion.sadcp'`

- [ ] **Step 3: Implement `src/ladcp/ingestion/sadcp.py`**

```python
"""SADCP (shipboard ADCP) NetCDF loader.

Reads OS75 JASADCP-format NetCDF files, windows to a cast time range,
averages water-column velocity to a single profile, and returns arrays
ready for compute_inverse().

Mirrors the logic of loadsadcp.m (LDEO_IX) and generalises
scripts/generate_sadcp_fixture.py.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import netCDF4
import numpy as np
from numpy.typing import NDArray


@dataclass
class SADCPProfile:
    """Cast-averaged SADCP profile, ready for compute_inverse()."""
    z: NDArray[np.float64]    # depth, m, positive downward
    u: NDArray[np.float64]    # eastward velocity, m/s
    v: NDArray[np.float64]    # northward velocity, m/s
    err: NDArray[np.float64]  # per-bin velocity error estimate, m/s


def _parse_epoch_jd(units: str) -> float:
    """Parse 'days since YYYY-MM-DD [...]' → Julian day of that epoch."""
    from ladcp.ingestion._pd0 import _to_julian  # reuse existing helper

    date_str = units.split("since")[-1].strip().split()[0]  # "YYYY-MM-DD"
    yyyy, mm, dd = (int(x) for x in date_str.split("-"))
    return _to_julian(yyyy, mm, dd, 0.0)


def load_sadcp_nc(
    path: str | Path,
    t_start_jd: float,
    t_end_jd: float,
    *,
    lat: float | None = None,
    lon: float | None = None,
    pos_tol: float = 0.1,
    err_default: float = 0.05,
    dt_slack_days: float = 0.0,
    max_depth: float | None = None,
) -> SADCPProfile | None:
    """Load an OS75 SADCP NetCDF and return a cast-averaged profile.

    Parameters
    ----------
    path:
        Path to the NetCDF file (JASADCP format with time/lon/lat/depth/u/v).
    t_start_jd, t_end_jd:
        Cast time window in Julian days (same epoch as ladcp.ingestion._pd0).
    lat, lon:
        Expected cast position in decimal degrees (N positive, E positive).
        If provided, raises ValueError if the mean SADCP position in the
        window differs by more than *pos_tol* degrees.
    pos_tol:
        Maximum allowed position mismatch in degrees (default 0.1°, matching
        loadsadcp.m).
    err_default:
        Fallback per-bin error when std dev cannot be computed (single record
        or std dev is zero).  Default 0.05 m/s (OS75 conservative estimate).
    dt_slack_days:
        Expand the time window symmetrically by this many days (default 0,
        matching loadsadcp.m ``p.sadcp_dtok=0``).
    max_depth:
        Discard bins below this depth (m).  Default: keep all.

    Returns
    -------
    SADCPProfile or None
        None if no records fall within the time window.

    Raises
    ------
    ValueError
        If *lat* / *lon* are given and the mean SADCP position is too far.
    """
    path = Path(path)
    ds = netCDF4.Dataset(str(path))
    try:
        time_var = ds.variables["time"]
        epoch_jd = _parse_epoch_jd(str(time_var.units))
        t_jd = np.asarray(time_var[:], dtype=float) + epoch_jd

        lon_raw = np.asarray(ds.variables["lon"][:], dtype=float)
        lon_norm = lon_raw % 360.0  # handle JASADCP negative-offset convention

        lat_raw = np.asarray(ds.variables["lat"][:], dtype=float)

        depth_raw = np.asarray(ds.variables["depth"][:], dtype=float)
        u_raw = np.ma.filled(ds.variables["u"][:], np.nan).astype(float)
        v_raw = np.ma.filled(ds.variables["v"][:], np.nan).astype(float)
    finally:
        ds.close()

    # --- Time window ---
    in_window = (t_jd >= t_start_jd - dt_slack_days) & (
        t_jd <= t_end_jd + dt_slack_days
    )
    n_rec = int(in_window.sum())
    if n_rec == 0:
        return None

    # --- Position sanity check (mirrors loadsadcp.m) ---
    if lat is not None and lon is not None:
        mean_lat = float(np.nanmean(lat_raw[in_window]))
        mean_lon = float(np.nanmean(lon_norm[in_window]))
        lat_err = abs(mean_lat - lat)
        lon_scale = max(np.cos(np.radians(lat)), 0.01)
        lon_err = abs(mean_lon - lon) / lon_scale
        if lat_err > pos_tol or lon_err > pos_tol:
            raise ValueError(
                f"SADCP position ({mean_lat:.3f}°N, {mean_lon:.3f}°E) differs "
                f"from station ({lat:.3f}°N, {lon:.3f}°E) by "
                f"{lat_err:.3f}°lat / {lon_err:.3f}°lon (tol {pos_tol}°)"
            )

    u_sel = u_raw[in_window]  # (n_rec, n_depth)
    v_sel = v_raw[in_window]

    u_avg = np.nanmean(u_sel, axis=0)  # (n_depth,)
    v_avg = np.nanmean(v_sel, axis=0)

    # --- Per-bin error (mirrors loadsadcp.m nstd + nvel scaling) ---
    if n_rec > 1:
        nvel = np.sum(np.isfinite(u_sel + v_sel), axis=0).astype(float)
        u_std = np.nanstd(u_sel, axis=0, ddof=1)
        v_std = np.nanstd(v_sel, axis=0, ddof=1)
        v_err = (u_std + v_std) / 2.0
        max_nvel = float(np.nanmax(nvel)) if np.any(nvel > 0) else 1.0
        # Scale: bins with fewer records get proportionally larger error
        with np.errstate(divide="ignore", invalid="ignore"):
            v_err = np.where(nvel > 0, v_err * max_nvel / nvel, err_default)
        v_err = np.where(v_err == 0.0, err_default, v_err)
    else:
        v_err = np.full_like(u_avg, err_default)

    # --- Depth array: handle (n_time, n_depth) or (n_depth,) layout ---
    if depth_raw.ndim == 2:
        z_avg = np.nanmean(depth_raw[in_window], axis=0)
    else:
        z_avg = depth_raw.copy()

    # --- Filter valid bins ---
    valid = (
        np.isfinite(u_avg)
        & np.isfinite(v_avg)
        & np.isfinite(z_avg)
        & np.isfinite(v_err)
    )
    if max_depth is not None:
        valid &= z_avg < max_depth
    if not np.any(valid):
        return None

    return SADCPProfile(
        z=z_avg[valid],
        u=u_avg[valid],
        v=v_avg[valid],
        err=v_err[valid],
    )
```

- [ ] **Step 4: Run tests**

```
pytest tests/unit/test_sadcp_loader.py -v
```

Expected: all 7 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ladcp/ingestion/sadcp.py tests/unit/test_sadcp_loader.py
git commit -m "feat: add SADCP NetCDF loader (load_sadcp_nc)"
```

---

### Task 2: S4P fixture generator script

**Files:**
- Create: `scripts/generate_sadcp_fixture_s4p.py`

**Interfaces:**
- Consumes: `load_sadcp_nc()` from Task 1; `001.nc` tim variable for cast time window
- Produces: `test_data/2018_S4P/sadcp_001.npz` (same schema as `sadcp_003.npz`)

- [ ] **Step 1: Write the script**

Create `scripts/generate_sadcp_fixture_s4p.py`:

```python
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
```

- [ ] **Step 2: Verify script runs (requires TEST_DATA_DIR)**

```bash
TEST_DATA_DIR="C:/Users/peter_sha/Documents/sourcecode/LADCP/test_data" \
    python scripts/generate_sadcp_fixture_s4p.py
```

Expected output for each cast: `Written .../sadcp_001.npz` with depth range 31–975 m.

- [ ] **Step 3: Commit**

```bash
git add scripts/generate_sadcp_fixture_s4p.py
git commit -m "feat: S4P SADCP fixture generator script"
```

---

### Task 3: S4P integration test

**Files:**
- Create: `tests/integration/test_sadcp_loader_s4p.py`

**Interfaces:**
- Consumes: `load_sadcp_nc()`, `001.nc` (tim, z_sadcp, u_sadcp, v_sadcp), `os75nb_short.nc`
- Produces: validated SADCP profile comparable to reference embedded profile

- [ ] **Step 1: Write the integration test**

Create `tests/integration/test_sadcp_loader_s4p.py`:

```python
"""Integration test: SADCP loader vs S4P cast 001 reference.

Requires TEST_DATA_DIR env var pointing to a directory containing 2018_S4P/.

Note: the reference 001.nc embeds u_sadcp/v_sadcp extracted from SADCP.mat
(a different product from os75nb_short.nc).  We compare mean absolute difference
rather than exact equality; tolerance 0.10 m/s.
"""
from __future__ import annotations

import os
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.sadcp import SADCPProfile, load_sadcp_nc


@pytest.fixture(scope="module")
def s4p_dir() -> Path:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        pytest.skip("TEST_DATA_DIR not set — see test_data/sources.md")
    p = Path(env) / "2018_S4P"
    if not p.exists():
        pytest.skip(f"2018_S4P directory not found at {p}")
    return p


@pytest.fixture(scope="module")
def sadcp_nc_path(s4p_dir: Path) -> Path:
    p = s4p_dir / "SADCP/os75nb/contour/os75nb_short.nc"
    if not p.exists():
        pytest.skip(f"SADCP NC not found: {p}")
    return p


@pytest.fixture(scope="module")
def ref_001(s4p_dir: Path) -> Path:
    p = s4p_dir / "001.nc"
    if not p.exists():
        pytest.skip(f"Reference 001.nc not found: {p}")
    return p


@pytest.fixture(scope="module")
def cast_001_window(ref_001: Path) -> tuple[float, float]:
    ds = netCDF4.Dataset(str(ref_001))
    tim = np.asarray(ds.variables["tim"][:], dtype=float)
    ds.close()
    return float(tim[0]), float(tim[-1])


@pytest.fixture(scope="module")
def sadcp_profile(sadcp_nc_path: Path, cast_001_window: tuple[float, float]) -> SADCPProfile:
    t_start, t_end = cast_001_window
    profile = load_sadcp_nc(
        sadcp_nc_path, t_start, t_end,
        lat=-70.45, lon=168.47,
    )
    if profile is None:
        pytest.skip("No SADCP records in cast 001 time window")
    return profile


@pytest.mark.integration
def test_sadcp_load_returns_profile(sadcp_profile: SADCPProfile) -> None:
    assert isinstance(sadcp_profile, SADCPProfile)


@pytest.mark.integration
def test_sadcp_profile_depth_range(sadcp_profile: SADCPProfile) -> None:
    """OS75 covers ~31–975 m; expect at least 100–500 m."""
    assert float(sadcp_profile.z.min()) < 100.0
    assert float(sadcp_profile.z.max()) > 500.0


@pytest.mark.integration
def test_sadcp_profile_all_finite(sadcp_profile: SADCPProfile) -> None:
    assert np.all(np.isfinite(sadcp_profile.u))
    assert np.all(np.isfinite(sadcp_profile.v))
    assert np.all(np.isfinite(sadcp_profile.err))
    assert np.all(sadcp_profile.err > 0)


@pytest.mark.integration
def test_sadcp_vs_reference_mae(sadcp_profile: SADCPProfile, ref_001: Path) -> None:
    """Compare against u_sadcp/v_sadcp embedded in reference NC.

    These came from SADCP.mat (different processing path), so we use a
    generous 0.10 m/s MAE tolerance rather than requiring exact agreement.
    """
    ds = netCDF4.Dataset(str(ref_001))
    ref_z = np.asarray(ds.variables["z_sadcp"][:], dtype=float)
    ref_u = np.asarray(ds.variables["u_sadcp"][:], dtype=float)
    ref_v = np.asarray(ds.variables["v_sadcp"][:], dtype=float)
    ds.close()

    our_u = np.interp(ref_z, sadcp_profile.z, sadcp_profile.u,
                      left=np.nan, right=np.nan)
    our_v = np.interp(ref_z, sadcp_profile.z, sadcp_profile.v,
                      left=np.nan, right=np.nan)

    valid = np.isfinite(ref_u) & np.isfinite(our_u)
    assert valid.sum() >= 10, f"Too few overlapping bins: {valid.sum()}"

    mae_u = float(np.mean(np.abs(our_u[valid] - ref_u[valid])))
    mae_v = float(np.mean(np.abs(our_v[valid] - ref_v[valid])))

    assert mae_u < 0.10, f"u MAE {mae_u:.4f} m/s exceeds 0.10 m/s (different SADCP products)"
    assert mae_v < 0.10, f"v MAE {mae_v:.4f} m/s exceeds 0.10 m/s"
```

- [ ] **Step 2: Run (skips cleanly without data)**

```
pytest tests/integration/test_sadcp_loader_s4p.py -v
```

Expected without `TEST_DATA_DIR`: all tests SKIP.

- [ ] **Step 3: Run with data**

```
TEST_DATA_DIR="C:/Users/peter_sha/Documents/sourcecode/LADCP/test_data" \
    pytest tests/integration/test_sadcp_loader_s4p.py -v -m integration
```

Expected: all 4 tests PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/integration/test_sadcp_loader_s4p.py
git commit -m "test: S4P SADCP loader integration test vs 001.nc reference"
```

---

## Self-Review

**Spec coverage:**
- ✅ `SADCPProfile` dataclass with z/u/v/err
- ✅ `load_sadcp_nc()` with time window, `lon % 360`, position check, error scaling
- ✅ `max_depth` filter for depth limiting
- ✅ S4P fixture generator script (generalises P16N approach)
- ✅ Integration test with 0.10 m/s MAE tolerance and explicit caveat on SADCP.mat vs os75nb difference
- ✅ Unit tests cover: nominal, out-of-window→None, position mismatch→ValueError, lon normalisation, max_depth, nvel-scaled error

**No placeholders present.**

**Type consistency:** `SADCPProfile` defined in Task 1, consumed by Tasks 2–3 using exact field names `z`, `u`, `v`, `err`.
