# GPS Barotropic + SADCP Velocity Constraint Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire GPS-derived ship velocity (barotropic constraint) and synthetic SADCP profile into the LADCP inverse solver integration test, replacing the current zero/NaN placeholders.

**Architecture:** Three changes in concert — (1) extend `CTDTimeSeries` to parse lat/lon columns embedded in SBE CNV files and add a `compute_ship_velocity()` function for linear-regression drift estimation, (2) create a one-off fixture-generator script that derives a synthetic SADCP profile from the `003.nc` reference output, (3) update the integration test to load both constraints and pass them to `compute_inverse()`.

**Tech Stack:** NumPy (array ops, `polyfit`), Python dataclasses, `np.savez`/`np.load` for the NPZ fixture, `netCDF4` for the fixture generator.

## Global Constraints

- Python ≥ 3.11; `NDArray[np.float64]` for all array type hints.
- SBE bad-flag sentinel `-9.990e-29` masked via `np.isclose(arr, bad_flag, rtol=1e-3, atol=0)` (existing whole-array mask applied before lat/lon extraction).
- Lat/lon header flag string: `Store Lat/Lon Data = Append to Every Scan` (appears on a `*`-prefixed line).
- Equirectangular displacement constant: `111320.0` m/degree.
- Unit tests: `tests/test_ctd_loader.py`.
- Integration test: `tests/integration/test_inverse_p16n_cast003.py`.
- SADCP fixture path: `$TEST_DATA_DIR/2015_P16N/sadcp_003.npz`.
- Unit test run command: `python -m pytest tests/test_ctd_loader.py -v`
- Full unit test run command: `python -m pytest tests/ -v --ignore=tests/integration`

---

## File Map

| Action | File | What changes |
|--------|------|-------------|
| Modify | `src/ladcp/ingestion/ctd.py` | Add `lat`/`lon` fields to `CTDTimeSeries`; detect header flag; slice columns in binary + ASCII readers; add `compute_ship_velocity()` |
| Modify | `src/ladcp/ingestion/__init__.py` | Export `compute_ship_velocity` |
| Create | `scripts/generate_sadcp_fixture.py` | One-off data prep: reads `003.nc`, writes `sadcp_003.npz` |
| Modify | `tests/test_ctd_loader.py` | 5 new unit tests (lat/lon parsing × 3, ship velocity × 2) |
| Modify | `tests/integration/test_inverse_p16n_cast003.py` | Wire GPS + SADCP into `inverse_result` fixture; update xfail reasons |

---

## Task 1: CTD lat/lon parsing

**Files:**
- Modify: `src/ladcp/ingestion/ctd.py`
- Modify: `tests/test_ctd_loader.py`

**Interfaces:**
- Produces: `CTDTimeSeries.lat: NDArray[np.float64] | None`, `CTDTimeSeries.lon: NDArray[np.float64] | None`
- Produces: `_parse_sbe_header` sets `header_info['lat_lon_appended']: bool`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ctd_loader.py`:

```python
def _make_binary_cnv_latlon(
    tmp_path: Path,
    lat_vals: list[float],
    lon_vals: list[float],
    bad_lon_idx: int | None = None,
) -> Path:
    """Binary CNV with 'Store Lat/Lon Data' flag and nquan=6."""
    BAD = -9.990e-29
    rows = []
    for i, (la, lo) in enumerate(zip(lat_vals, lon_vals)):
        lo_val = BAD if (bad_lon_idx is not None and i == bad_lon_idx) else lo
        rows.append([100.0 + i * 100, 15.0 - i, 35.0 + i * 0.1,
                     2451545.0 + i * 0.1, la, lo_val])
    data = np.array(rows, dtype=np.float64)
    nquan = 6
    header = (
        f"# nquan = {nquan}\n"
        "# name 0 = prDM: Pressure [db]\n"
        "# name 1 = t090C: Temperature\n"
        "# name 2 = sal00: Salinity\n"
        "# name 3 = timeJ: Julian Days\n"
        "# name 4 = latitude: Latitude\n"
        "# name 5 = longitude: Longitude\n"
        f"# bad_flag = {BAD}\n"
        "# file_type = binary\n"
        "* Store Lat/Lon Data = Append to Every Scan\n"
        "*END*\n"
    ).encode()
    p = tmp_path / "test_latlon.cnv"
    p.write_bytes(header + data.astype("<f4").tobytes())
    return p


def test_ctd_loads_lat_lon(tmp_path: Path):
    """CNV with 'Store Lat/Lon Data' flag populates lat/lon of correct length."""
    p = _make_binary_cnv_latlon(tmp_path, [-30.0, -30.01], [-140.0, -139.99])
    result = load_ctd(p)
    assert result.lat is not None
    assert result.lon is not None
    assert len(result.lat) == 2
    assert len(result.lon) == 2
    assert abs(result.lat[0] - (-30.0)) < 0.01
    assert abs(result.lon[0] - (-140.0)) < 0.01
    # Standard columns must still be populated correctly
    assert abs(result.pressure_dbar[0] - 100.0) < 1.0
    assert abs(result.pressure_dbar[1] - 200.0) < 1.0


def test_ctd_no_lat_lon_returns_none(tmp_path: Path):
    """CNV without 'Store Lat/Lon Data' flag returns lat=None, lon=None."""
    data = np.array([[100.0, 15.0, 35.0, 2451545.0]], dtype=np.float64)
    p = _make_binary_cnv(tmp_path, data)
    result = load_ctd(p)
    assert result.lat is None
    assert result.lon is None


def test_ctd_bad_flag_lat_lon_becomes_nan(tmp_path: Path):
    """SBE bad-flag sentinel in a lat/lon column becomes NaN."""
    p = _make_binary_cnv_latlon(
        tmp_path,
        [-30.0, -30.01],
        [-140.0, -139.99],
        bad_lon_idx=0,  # row 0 lon is bad
    )
    result = load_ctd(p)
    assert result.lon is not None
    assert np.isnan(result.lon[0])
    assert not np.isnan(result.lon[1])
    assert not np.isnan(result.lat[0])
```

- [ ] **Step 2: Run tests to verify they fail**

```
python -m pytest tests/test_ctd_loader.py::test_ctd_loads_lat_lon tests/test_ctd_loader.py::test_ctd_no_lat_lon_returns_none tests/test_ctd_loader.py::test_ctd_bad_flag_lat_lon_becomes_nan -v
```

Expected: `AttributeError: 'CTDTimeSeries' object has no attribute 'lat'` (or similar).

- [ ] **Step 3: Implement the changes in `src/ladcp/ingestion/ctd.py`**

**3a. Add fields to `CTDTimeSeries`** (after the `salinity` field):

```python
@dataclass
class CTDTimeSeries:
    time_julian: NDArray[np.float64]    # (nctd,) Julian days
    pressure_dbar: NDArray[np.float64]  # (nctd,) positive down
    temp_c: NDArray[np.float64]         # (nctd,) NaN if absent
    salinity: NDArray[np.float64]       # (nctd,) NaN if absent
    lat: NDArray[np.float64] | None = None  # (nctd,) degrees N; NaN = bad fix
    lon: NDArray[np.float64] | None = None  # (nctd,) degrees E; NaN = bad fix
```

**3b. In `_parse_sbe_header`, add `lat_lon_appended` to result dict and detect the flag.**

In the `result` dict initialization, add:
```python
'lat_lon_appended': False,
```

In the `for line in lines:` loop, before `if not line.startswith('#'): continue`, add:
```python
if 'Store Lat/Lon Data = Append to Every Scan' in line:
    result['lat_lon_appended'] = True
```

Full updated loop block:
```python
for line in lines:
    if line.strip() == '*END*':
        break
    if 'Store Lat/Lon Data = Append to Every Scan' in line:
        result['lat_lon_appended'] = True
    if not line.startswith('#'):
        continue
    content = line[1:].strip()
    if content.startswith('nquan ='):
        result['nquan'] = int(content.split('=', 1)[1].strip())
    elif m := re.match(r'name (\d+) = (\w+)', content):
        result['columns'][int(m.group(1))] = m.group(2)
    elif content.startswith('bad_flag ='):
        result['bad_flag'] = float(content.split('=', 1)[1].strip())
    elif 'file_type = binary' in content:
        result['file_type'] = 'binary'
    elif content.startswith('start_time ='):
        raw_date = content.split('=', 1)[1].strip()
        result['start_time_julian'] = _parse_start_time(raw_date)
        try:
            dt = datetime.strptime(raw_date.strip()[:20], '%b %d %Y %H:%M:%S')
            result['start_year'] = dt.year
        except ValueError:
            pass
```

**3c. Add a `_extract_latlon` helper** (place just before `_build_ctd_time_series`):

```python
def _extract_latlon(
    arr: np.ndarray,
    col_roles: dict[int, str | None],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[int, str | None]]:
    """Slice lat/lon off the last two columns; return trimmed arr and col_roles."""
    lat = arr[:, -2].copy()
    lon = arr[:, -1].copy()
    arr = arr[:, :-2]
    n_std = arr.shape[1]
    col_roles = {k: v for k, v in col_roles.items() if k < n_std}
    return lat, lon, arr, col_roles
```

**3d. Update `_build_ctd_time_series` signature and return** — add `lat` and `lon` keyword args:

```python
def _build_ctd_time_series(
    arr: np.ndarray,
    col_roles: dict[int, str | None],
    time_start_julian: float | None,
    start_year: int | None = None,
    lat: np.ndarray | None = None,
    lon: np.ndarray | None = None,
) -> CTDTimeSeries:
```

Change the final return to:
```python
    return CTDTimeSeries(
        time_julian=time_julian,
        pressure_dbar=pressure,
        temp_c=temp if temp is not None else np.full(n, np.nan),
        salinity=salinity if salinity is not None else np.full(n, np.nan),
        lat=lat,
        lon=lon,
    )
```

**3e. Update `_read_sbe_binary`** — after the bad_flag masking block and before the return, add:

```python
    lat = lon = None
    if header_info.get('lat_lon_appended') and arr.shape[1] >= 2:
        lat, lon, arr, col_roles = _extract_latlon(arr, col_roles)

    return _build_ctd_time_series(
        arr, col_roles, time_start_julian,
        start_year=header_info.get('start_year'),
        lat=lat, lon=lon,
    )
```

**3f. Update `_read_sbe_ascii`** — same pattern after the bad_flag masking block and before the return:

```python
    lat = lon = None
    if header_info.get('lat_lon_appended') and arr.shape[1] >= 2:
        lat, lon, arr, col_roles = _extract_latlon(arr, col_roles)

    return _build_ctd_time_series(
        arr, col_roles, time_start_julian,
        start_year=header_info.get('start_year'),
        lat=lat, lon=lon,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

```
python -m pytest tests/test_ctd_loader.py -v
```

Expected: all existing tests pass + 3 new tests pass. Zero failures.

- [ ] **Step 5: Commit**

```bash
git add src/ladcp/ingestion/ctd.py tests/test_ctd_loader.py
git commit -m "feat: parse lat/lon columns from SBE CNV files; extend CTDTimeSeries"
```

---

## Task 2: `compute_ship_velocity` function + export

**Files:**
- Modify: `src/ladcp/ingestion/ctd.py`
- Modify: `src/ladcp/ingestion/__init__.py`
- Modify: `tests/test_ctd_loader.py`

**Interfaces:**
- Consumes: `CTDTimeSeries.lat`, `CTDTimeSeries.lon`, `CTDTimeSeries.time_julian`
- Produces: `compute_ship_velocity(lat, lon, time_jul) -> tuple[float, float]` — `(u_ship_mps, v_ship_mps)`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ctd_loader.py`:

```python
from ladcp.ingestion.ctd import compute_ship_velocity


def test_compute_ship_velocity_linear_track():
    """Straight eastward track at 1.0 m/s returns (u≈1.0, v≈0.0)."""
    import math
    # At lat=0 deg, 1 deg longitude = 111320 m. Moving east at 1 m/s:
    # Δlon_deg/s = 1.0 / 111320
    lat0 = 0.0
    speed_mps = 1.0
    dt_s = np.array([0.0, 100.0, 200.0, 300.0])
    east_m = speed_mps * dt_s
    lon = lon0 = -140.0
    lon_arr = lon0 + east_m / (math.cos(math.radians(lat0)) * 111320.0)
    lat_arr = np.full(4, lat0)
    t0 = 2457100.0  # arbitrary Julian day
    time_jul = t0 + dt_s / 86400.0
    u_ship, v_ship = compute_ship_velocity(lat_arr, lon_arr, time_jul)
    assert abs(u_ship - 1.0) < 0.01, f"u_ship={u_ship:.4f} expected ≈1.0"
    assert abs(v_ship) < 0.01, f"v_ship={v_ship:.4f} expected ≈0.0"


def test_compute_ship_velocity_insufficient_data():
    """All-NaN lat/lon returns (0.0, 0.0)."""
    lat = np.full(5, np.nan)
    lon = np.full(5, np.nan)
    t0 = 2457100.0
    time_jul = t0 + np.arange(5) / 86400.0
    u_ship, v_ship = compute_ship_velocity(lat, lon, time_jul)
    assert u_ship == 0.0
    assert v_ship == 0.0
```

- [ ] **Step 2: Run tests to verify they fail**

```
python -m pytest tests/test_ctd_loader.py::test_compute_ship_velocity_linear_track tests/test_ctd_loader.py::test_compute_ship_velocity_insufficient_data -v
```

Expected: `ImportError: cannot import name 'compute_ship_velocity'`.

- [ ] **Step 3: Implement `compute_ship_velocity` in `src/ladcp/ingestion/ctd.py`**

Add after `assign_bin_depths` (at the end of the file):

```python
def compute_ship_velocity(
    lat: NDArray[np.float64],
    lon: NDArray[np.float64],
    time_jul: NDArray[np.float64],
) -> tuple[float, float]:
    """Estimate mean ship velocity from GPS fixes via linear regression.

    Returns (u_ship, v_ship) in m/s (eastward, northward).
    Returns (0.0, 0.0) when fewer than 2 valid fixes exist.
    """
    valid = np.isfinite(lat) & np.isfinite(lon)
    if int(valid.sum()) < 2:
        return 0.0, 0.0
    lat_v = lat[valid]
    lon_v = lon[valid]
    t_v = time_jul[valid]
    lat0 = float(lat_v[0])
    east_m = (lon_v - lon_v[0]) * math.cos(math.radians(lat0)) * 111320.0
    north_m = (lat_v - lat_v[0]) * 111320.0
    t_s = (t_v - t_v[0]) * 86400.0
    u_ship = float(np.polyfit(t_s, east_m, 1)[0])
    v_ship = float(np.polyfit(t_s, north_m, 1)[0])
    return u_ship, v_ship
```

- [ ] **Step 4: Export from `src/ladcp/ingestion/__init__.py`**

Replace the file contents with:

```python
from ladcp.ingestion.ctd import compute_ship_velocity

__all__ = ["compute_ship_velocity"]
```

- [ ] **Step 5: Run tests to verify they pass**

```
python -m pytest tests/test_ctd_loader.py -v
```

Expected: all tests pass (original + 5 new). Zero failures.

- [ ] **Step 6: Commit**

```bash
git add src/ladcp/ingestion/ctd.py src/ladcp/ingestion/__init__.py tests/test_ctd_loader.py
git commit -m "feat: add compute_ship_velocity(); export from ingestion package"
```

---

## Task 3: SADCP fixture generator script

**Files:**
- Create: `scripts/generate_sadcp_fixture.py`

**Interfaces:**
- Consumes: `$TEST_DATA_DIR/2015_P16N/003.nc` (LDEO_IX reference output with variables `z`, `u`, `v`, `nvel`)
- Produces: `$TEST_DATA_DIR/2015_P16N/sadcp_003.npz` (keys: `z`, `u`, `v`, `err`)

- [ ] **Step 1: Create `scripts/` directory and write the script**

Create `scripts/generate_sadcp_fixture.py`:

```python
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
```

- [ ] **Step 2: Run the script (requires TEST_DATA_DIR)**

```bash
python scripts/generate_sadcp_fixture.py
```

Expected output (exact numbers will vary):
```
Written: /path/to/2015_P16N/sadcp_003.npz
  N depth bins, depth range A–B m
```

If TEST_DATA_DIR is not set, skip this step and note the fixture is optional in the integration test (it's guarded by `if sadcp_path.exists()`).

- [ ] **Step 3: Commit**

```bash
git add scripts/generate_sadcp_fixture.py
git commit -m "feat: add SADCP fixture generator script for P16N cast 003"
```

---

## Task 4: Wire integration test

**Files:**
- Modify: `tests/integration/test_inverse_p16n_cast003.py`

**Interfaces:**
- Consumes: `CTDTimeSeries.lat`, `CTDTimeSeries.lon` (from Task 1)
- Consumes: `compute_ship_velocity()` (from Task 2)
- Consumes: `sadcp_003.npz` (from Task 3, optional)

- [ ] **Step 1: Update the `inverse_result` fixture**

Replace the entire `inverse_result` fixture in `tests/integration/test_inverse_p16n_cast003.py` with:

```python
@pytest.fixture(scope="module")
def inverse_result(dl_path: Path, cnv_path: Path, test_data_dir: Path) -> InverseResult:
    """Run full pipeline on P16N cast 003 raw data."""
    from ladcp.ingestion.ctd import compute_ship_velocity

    rdi = load_rdi(dl_path)
    ctd = load_ctd(cnv_path)

    # beam2earth: file is in beam coordinates (EX byte 0x04), need explicit transform.
    # The rdi.u/v/w/e fields hold raw beam data for beam-coord files.
    u_earth, v_earth, w_earth = beam2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG,
        gimbaled=True,
    )

    # assign_bin_depths returns positive-down z_m (nens,) and izm (nbin, nens).
    z_m, izm_pos = assign_bin_depths(rdi, ctd, looker="down")

    # EnsembleData depth convention: negative = below surface.
    z_neg = -z_m           # (nens,) negative down
    izm_neg = -izm_pos     # (nbin, nens) negative down

    # Correlation-based weight: mean over 4 beams, normalised to 0–1.
    weight = np.nanmean(rdi.corr.astype(np.float64), axis=2) / 128.0

    # Bottom-track velocity in Earth frame.  btrack_vel_ms is (4, nens) beam-frame;
    # apply the same beam→earth rotation used for water-track data.
    bt_u_e, bt_v_e, bt_w_e = beam2earth(
        rdi.btrack_vel_ms[0],
        rdi.btrack_vel_ms[1],
        rdi.btrack_vel_ms[2],
        rdi.btrack_vel_ms[3],
        rdi.heading,
        rdi.pitch,
        rdi.roll,
        THETA_DEG,
        gimbaled=True,
    )
    bvel = np.stack([bt_u_e, bt_v_e, bt_w_e], axis=1)  # (nens, 3) Earth frame
    bvels = np.full_like(bvel, 0.02)                      # 2 cm/s nominal std
    hbot = np.nanmean(rdi.btrack_range_m, axis=0)  # (nens,) mean of 4-beam ranges

    # GPS nav interpolated onto ensemble times
    if ctd.lat is not None:
        slat = np.interp(
            rdi.time_julian, ctd.time_julian, ctd.lat, left=np.nan, right=np.nan
        )
        slon = np.interp(
            rdi.time_julian, ctd.time_julian, ctd.lon, left=np.nan, right=np.nan
        )
        u_ship, v_ship = compute_ship_velocity(ctd.lat, ctd.lon, ctd.time_julian)
    else:
        slat = np.full(rdi.nens, np.nan)
        slon = np.full(rdi.nens, np.nan)
        u_ship, v_ship = 0.0, 0.0

    # SADCP fixture (optional — skipped gracefully if not generated yet)
    sadcp_path = test_data_dir / "sadcp_003.npz"
    if sadcp_path.exists():
        npz = np.load(sadcp_path)
        sadcp_z, sadcp_u, sadcp_v, sadcp_err = (
            npz["z"], npz["u"], npz["v"], npz["err"]
        )
    else:
        sadcp_z = sadcp_u = sadcp_v = sadcp_err = None

    ens = EnsembleData(
        u=u_earth,
        v=v_earth,
        w=w_earth,
        weight=weight,
        izm=izm_neg,
        z=z_neg,
        time_jul=rdi.time_julian,
        bvel=bvel,
        bvels=bvels,
        hbot=hbot,
        izd=np.arange(rdi.nbin),
        izu=np.array([], dtype=int),
        slat=slat,
        slon=slon,
    )

    ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)

    se = prepare_superensembles(ens, dz=16.0)
    return compute_inverse(
        se,
        u_ship=u_ship,
        v_ship=v_ship,
        sadcp_z=sadcp_z,
        sadcp_u=sadcp_u,
        sadcp_v=sadcp_v,
        sadcp_err=sadcp_err,
    )
```

- [ ] **Step 2: Update xfail reason strings**

Replace both xfail decorators (on `test_inverse_u_rmse` and `test_inverse_v_rmse`):

Change:
```python
@pytest.mark.xfail(strict=False, reason="Pipeline gap: no GPS/SADCP constraint; remove once pipeline complete")
```

To:
```python
@pytest.mark.xfail(strict=False, reason="RMSE target not yet met; remove once pipeline complete")
```

- [ ] **Step 3: Run the unit test suite to verify no regressions**

```
python -m pytest tests/ --ignore=tests/integration -v
```

Expected: all tests pass (≥ 95 tests, zero failures).

- [ ] **Step 4: Run integration tests if TEST_DATA_DIR is available**

```
python -m pytest tests/integration/ -v
```

Expected: non-RMSE integration tests pass; RMSE tests are `xfail` (may XPASS or XFAIL depending on GPS/SADCP quality). No `FAILED` results.

- [ ] **Step 5: Commit**

```bash
git add tests/integration/test_inverse_p16n_cast003.py
git commit -m "feat: wire GPS barotropic + SADCP constraints into integration test"
```

---

## Self-Review

**Spec coverage check:**

| Spec requirement | Task covering it |
|-----------------|-----------------|
| `CTDTimeSeries.lat`/`.lon` optional fields | Task 1 Step 3a |
| Detect `Store Lat/Lon Data = Append to Every Scan` | Task 1 Step 3b |
| Replace SBE sentinel with NaN | Task 1 Step 3e/3f (via existing bad_flag mask applied before slicing) |
| `compute_ship_velocity()` with polyfit | Task 2 Step 3 |
| Export from `ladcp.ingestion` | Task 2 Step 4 |
| `scripts/generate_sadcp_fixture.py` | Task 3 |
| Interpolate CTD lat/lon onto ensemble times | Task 4 Step 1 |
| Pass `u_ship`/`v_ship` to `compute_inverse` | Task 4 Step 1 |
| Load SADCP fixture (optional) | Task 4 Step 1 |
| Update xfail reason strings | Task 4 Step 2 |
| Unit test: CNV with flag → non-NaN lat/lon | Task 1 Step 1 |
| Unit test: CNV without flag → None | Task 1 Step 1 |
| Unit test: linear track at known speed | Task 2 Step 1 |
| Unit test: all-NaN → (0.0, 0.0) | Task 2 Step 1 |

**Placeholder scan:** No TBDs, no "add appropriate error handling", all code blocks are complete.

**Type consistency:**
- `CTDTimeSeries.lat: NDArray[np.float64] | None` — used in Task 4 as `ctd.lat`  ✓
- `compute_ship_velocity(lat, lon, time_jul) -> tuple[float, float]` — used in Task 4 ✓
- `_extract_latlon(arr, col_roles) -> (lat, lon, arr, col_roles)` — called in Task 1 Steps 3e/3f ✓
- `compute_inverse(se, u_ship=..., v_ship=..., sadcp_z=..., ...)` — matches existing signature ✓
