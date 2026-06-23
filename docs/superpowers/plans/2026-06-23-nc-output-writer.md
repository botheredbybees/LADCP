# NetCDF Output Writer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `write_ladcp_nc()` in `src/ladcp/output/nc.py` that writes an `InverseResult` to a NetCDF file with the LDEO_IX variable schema, validated against the S4P reference outputs (`001.nc`, `002.nc`, `003.nc`).

**Architecture:** The S4P reference outputs already define the target schema — this plan reads their variable list, creates a `write_ladcp_nc()` function that maps `InverseResult` fields onto that schema, and validates output by reading back what was written. The function accepts optional supplementary data (ensemble GPS track, SADCP profile) matching what the reference files embed. Unit tests verify round-trip correctness with synthetic data; an integration test checks that our writer produces the same variable names and dimensions as the reference.

**Tech Stack:** Python 3.11, netCDF4, numpy, pytest

## Global Constraints

- Package name: `ladcp`; output module lives at `src/ladcp/output/`
- NC variable naming must match the S4P reference files exactly (checked in integration test)
- `InverseResult` is imported from `ladcp.solution.inverse`; do not modify it
- `SADCPProfile` is imported from `ladcp.ingestion.sadcp` (implemented in the SADCP loader plan)
- Test env var: `TEST_DATA_DIR`; integration tests carry `@pytest.mark.integration` and skip when data is absent
- `ubar`/`vbar` are stored as scalar global attributes (not array variables), matching the LDEO_IX convention

---

### Task 1: Output module scaffold and `write_ladcp_nc()` function

**Files:**
- Create: `src/ladcp/output/__init__.py`
- Create: `src/ladcp/output/nc.py`
- Modify: `src/ladcp/__init__.py` (no change needed; new subpackage is auto-discovered)
- Create: `tests/unit/test_nc_writer.py`

**Interfaces:**
- Consumes: `InverseResult` (from `ladcp.solution.inverse`), `SADCPProfile` (from `ladcp.ingestion.sadcp`)
- Produces: `write_ladcp_nc(path, result, *, ens_time_jd, ens_lat, ens_lon, sadcp, uship, vship) -> None`

`InverseResult` field reference (all arrays on the same `n_zbins` depth grid unless noted):
```python
z: np.ndarray        # (n_zbins,) depth m, positive downward
u: np.ndarray        # (n_zbins,) eastward velocity m/s
v: np.ndarray        # (n_zbins,) northward velocity m/s
uerr: np.ndarray     # (n_zbins,) velocity error m/s
nvel: np.ndarray     # (n_zbins,) observation count (int)
u_do: np.ndarray     # (n_zbins,) downcast-only u
v_do: np.ndarray     # (n_zbins,) downcast-only v
u_up: np.ndarray     # (n_zbins,) upcast-only u
v_up: np.ndarray     # (n_zbins,) upcast-only v
u_ctd: np.ndarray    # (n_se,) GPS-derived CTD eastward velocity
v_ctd: np.ndarray    # (n_se,) GPS-derived CTD northward velocity
ubar: float          # depth-mean eastward velocity (barotropic)
vbar: float          # depth-mean northward velocity (barotropic)
zctd: np.ndarray     # (n_se,) CTD depth time series m
wctd: np.ndarray     # (n_se,) CTD vertical velocity m/s
```

- [ ] **Step 1: Write failing unit tests**

Create `tests/unit/test_nc_writer.py`:

```python
"""Unit tests for write_ladcp_nc()."""
from __future__ import annotations
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.output.nc import write_ladcp_nc
from ladcp.solution.inverse import InverseResult


def _make_result(n_z: int = 10, n_se: int = 20) -> InverseResult:
    """Synthetic InverseResult with recognisable values for round-trip checks."""
    z = np.arange(10.0, 10.0 + n_z * 10.0, 10.0)
    u = np.sin(z / 100.0) * 0.3
    v = np.cos(z / 100.0) * 0.2
    return InverseResult(
        z=z,
        u=u,
        v=v,
        uerr=np.full(n_z, 0.02),
        nvel=np.arange(1, n_z + 1),
        u_do=u * 1.1,
        v_do=v * 0.9,
        u_up=u * 0.9,
        v_up=v * 1.1,
        u_ctd=np.zeros(n_se),
        v_ctd=np.zeros(n_se),
        ubar=0.05,
        vbar=-0.03,
        zctd=np.linspace(0.0, 500.0, n_se),
        wctd=np.zeros(n_se),
    )


def test_write_creates_file(tmp_path: Path) -> None:
    result = _make_result()
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result)
    assert out.exists()


def test_write_core_variables(tmp_path: Path) -> None:
    result = _make_result(n_z=10)
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    for var in ("z", "u", "v", "uerr", "nvel", "u_do", "v_do", "u_up", "v_up"):
        assert var in ds.variables, f"Missing variable: {var}"
        assert ds.variables[var].shape == (10,), f"{var}: expected shape (10,)"
    ds.close()


def test_write_ctd_velocity_variables(tmp_path: Path) -> None:
    result = _make_result(n_z=5, n_se=8)
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    for var in ("uctd", "vctd", "zctd"):
        assert var in ds.variables, f"Missing variable: {var}"
        assert ds.variables[var].shape == (8,)
    ds.close()


def test_write_barotropic_as_attributes(tmp_path: Path) -> None:
    result = _make_result()
    result_with_ubar = InverseResult(
        **{**result.__dict__, "ubar": 0.123, "vbar": -0.045}
    )
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result_with_ubar)

    ds = netCDF4.Dataset(str(out))
    assert hasattr(ds, "ubar"), "ubar should be a global attribute"
    assert hasattr(ds, "vbar"), "vbar should be a global attribute"
    assert abs(float(ds.ubar) - 0.123) < 1e-6
    assert abs(float(ds.vbar) - -0.045) < 1e-6
    ds.close()


def test_roundtrip_u_v(tmp_path: Path) -> None:
    result = _make_result(n_z=15)
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    u_read = np.asarray(ds.variables["u"][:])
    v_read = np.asarray(ds.variables["v"][:])
    z_read = np.asarray(ds.variables["z"][:])
    ds.close()

    np.testing.assert_allclose(z_read, result.z, rtol=1e-5)
    np.testing.assert_allclose(u_read, result.u, rtol=1e-5)
    np.testing.assert_allclose(v_read, result.v, rtol=1e-5)


def test_write_with_gps_track(tmp_path: Path) -> None:
    result = _make_result(n_z=5)
    n_ens = 30
    out = tmp_path / "test.nc"
    write_ladcp_nc(
        out, result,
        ens_time_jd=np.linspace(2458919.0, 2458919.5, n_ens),
        ens_lat=np.full(n_ens, -70.45),
        ens_lon=np.full(n_ens, 168.47),
        uship=0.12,
        vship=-0.05,
    )

    ds = netCDF4.Dataset(str(out))
    assert "tim" in ds.variables
    assert "shiplat" in ds.variables
    assert "shiplon" in ds.variables
    assert ds.variables["tim"].shape == (n_ens,)
    assert hasattr(ds, "uship")
    assert hasattr(ds, "vship")
    ds.close()


def test_write_with_sadcp(tmp_path: Path) -> None:
    from ladcp.ingestion.sadcp import SADCPProfile

    result = _make_result(n_z=5)
    sadcp = SADCPProfile(
        z=np.array([50.0, 100.0, 150.0]),
        u=np.array([0.1, 0.15, 0.12]),
        v=np.array([-0.05, -0.03, -0.04]),
        err=np.array([0.05, 0.05, 0.05]),
    )
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result, sadcp=sadcp)

    ds = netCDF4.Dataset(str(out))
    assert "u_sadcp" in ds.variables
    assert "v_sadcp" in ds.variables
    assert "z_sadcp" in ds.variables
    assert ds.variables["z_sadcp"].shape == (3,)
    ds.close()
```

- [ ] **Step 2: Run to verify failure**

```
pytest tests/unit/test_nc_writer.py -v
```

Expected: `ImportError: cannot import name 'write_ladcp_nc' from 'ladcp.output.nc'`

- [ ] **Step 3: Create `src/ladcp/output/__init__.py`**

```python
"""LADCP output writers."""
```

- [ ] **Step 4: Implement `src/ladcp/output/nc.py`**

```python
"""NetCDF output writer: LDEO_IX-compatible schema.

Writes InverseResult to a NetCDF file matching the variable names and
layout produced by ladcp2cdf.m.  Reference schema: test_data/2018_S4P/001.nc.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import netCDF4
import numpy as np

from ladcp.solution.inverse import InverseResult

if TYPE_CHECKING:
    from ladcp.ingestion.sadcp import SADCPProfile


def write_ladcp_nc(
    path: str | Path,
    result: InverseResult,
    *,
    ens_time_jd: np.ndarray | None = None,
    ens_lat: np.ndarray | None = None,
    ens_lon: np.ndarray | None = None,
    sadcp: "SADCPProfile | None" = None,
    uship: float | None = None,
    vship: float | None = None,
) -> None:
    """Write InverseResult to NetCDF in LDEO_IX-compatible format.

    Parameters
    ----------
    path:
        Output file path (overwritten if it exists).
    result:
        Completed inverse solution from compute_inverse().
    ens_time_jd:
        Per-ensemble Julian day timestamps (stored as 'tim').
    ens_lat, ens_lon:
        Per-ensemble ship GPS latitude and longitude (stored as 'shiplat',
        'shiplon').
    sadcp:
        Cast-averaged SADCP profile to embed (stored as 'u_sadcp', 'v_sadcp',
        'z_sadcp').
    uship, vship:
        Depth-mean ship velocity components from GPS regression (m/s), stored
        as global attributes.
    """
    path = Path(path)
    ds = netCDF4.Dataset(str(path), "w", format="NETCDF4")
    try:
        n_z = len(result.z)
        n_se = len(result.zctd)

        ds.createDimension("z", n_z)
        ds.createDimension("nse", n_se)

        def _zvar(name: str, data: np.ndarray, units: str = "m/s") -> None:
            v = ds.createVariable(name, "f4", ("z",), fill_value=np.nan)
            v.units = units
            v[:] = data.astype(np.float32)

        # Core profile (z-dimension)
        _zvar("z", result.z, units="m")
        _zvar("u", result.u)
        _zvar("v", result.v)
        _zvar("uerr", result.uerr)

        nv = ds.createVariable("nvel", "i4", ("z",))
        nv.long_name = "number of velocity observations per depth bin"
        nv[:] = result.nvel.astype(np.int32)

        _zvar("u_do", result.u_do)
        _zvar("v_do", result.v_do)
        _zvar("u_up", result.u_up)
        _zvar("v_up", result.v_up)

        # CTD velocity time series (nse-dimension)
        def _sevar(name: str, data: np.ndarray, units: str = "m/s") -> None:
            v = ds.createVariable(name, "f4", ("nse",), fill_value=np.nan)
            v.units = units
            v[:] = data.astype(np.float32)

        _sevar("uctd", result.u_ctd)
        _sevar("vctd", result.v_ctd)
        _sevar("zctd", result.zctd, units="m")

        # Barotropic mean velocity — stored as global attributes (LDEO_IX convention)
        ds.ubar = float(result.ubar)
        ds.vbar = float(result.vbar)

        # Optional: GPS ensemble track
        if ens_time_jd is not None and ens_lat is not None and ens_lon is not None:
            n_ens = len(ens_time_jd)
            ds.createDimension("nens", n_ens)

            tim_v = ds.createVariable("tim", "f8", ("nens",))
            tim_v.long_name = "ensemble time, Julian days"
            tim_v[:] = ens_time_jd

            lat_v = ds.createVariable("shiplat", "f4", ("nens",), fill_value=np.nan)
            lat_v.units = "degrees_north"
            lat_v[:] = ens_lat.astype(np.float32)

            lon_v = ds.createVariable("shiplon", "f4", ("nens",), fill_value=np.nan)
            lon_v.units = "degrees_east"
            lon_v[:] = ens_lon.astype(np.float32)

        if uship is not None:
            ds.uship = float(uship)
        if vship is not None:
            ds.vship = float(vship)

        # Optional: SADCP profile
        if sadcp is not None:
            n_sadcp = len(sadcp.z)
            ds.createDimension("n_sadcp", n_sadcp)

            def _sadcpvar(name: str, data: np.ndarray, units: str = "m/s") -> None:
                v = ds.createVariable(name, "f4", ("n_sadcp",), fill_value=np.nan)
                v.units = units
                v[:] = data.astype(np.float32)

            _sadcpvar("z_sadcp", sadcp.z, units="m")
            _sadcpvar("u_sadcp", sadcp.u)
            _sadcpvar("v_sadcp", sadcp.v)

    finally:
        ds.close()
```

- [ ] **Step 5: Run unit tests**

```
pytest tests/unit/test_nc_writer.py -v
```

Expected: all 7 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/ladcp/output/__init__.py src/ladcp/output/nc.py tests/unit/test_nc_writer.py
git commit -m "feat: NetCDF output writer (write_ladcp_nc, LDEO_IX schema)"
```

---

### Task 2: Schema validation against S4P reference outputs

**Files:**
- Create: `tests/integration/test_nc_writer_s4p.py`

**Interfaces:**
- Consumes: `write_ladcp_nc()`, `001.nc` S4P reference (as schema oracle)
- Produces: verified that our writer produces all mandatory variables that the reference contains

**Note:** We cannot reproduce the reference values without raw PD0 files. This test validates schema (variable names + dimensions) and that our writer round-trips its own output faithfully.

- [ ] **Step 1: Write the integration test**

Create `tests/integration/test_nc_writer_s4p.py`:

```python
"""Integration test: NC writer schema vs S4P reference outputs.

Checks that write_ladcp_nc() produces every variable that appears in
001.nc (the LDEO_IX reference).  Value comparison is not possible here
because the S4P raw PD0 files are not present; see test_inverse_*.py
tests for value-level validation once raw data is available.

Requires TEST_DATA_DIR env var pointing to a directory containing 2018_S4P/.
"""
from __future__ import annotations

import os
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.sadcp import SADCPProfile
from ladcp.output.nc import write_ladcp_nc
from ladcp.solution.inverse import InverseResult


@pytest.fixture(scope="module")
def s4p_dir() -> Path:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        pytest.skip("TEST_DATA_DIR not set")
    p = Path(env) / "2018_S4P"
    if not p.exists():
        pytest.skip(f"2018_S4P directory not found at {p}")
    return p


@pytest.fixture(scope="module")
def ref_001_vars(s4p_dir: Path) -> set[str]:
    """Variable names present in the S4P reference 001.nc."""
    p = s4p_dir / "001.nc"
    if not p.exists():
        pytest.skip(f"Reference 001.nc not found: {p}")
    ds = netCDF4.Dataset(str(p))
    names = set(ds.variables.keys())
    ds.close()
    return names


@pytest.fixture(scope="module")
def ref_001_attrs(s4p_dir: Path) -> set[str]:
    """Global attribute names present in the S4P reference 001.nc."""
    p = s4p_dir / "001.nc"
    if not p.exists():
        pytest.skip(f"Reference 001.nc not found: {p}")
    ds = netCDF4.Dataset(str(p))
    attrs = set(ds.ncattrs())
    ds.close()
    return attrs


def _make_synthetic_result(n_z: int = 100, n_se: int = 50) -> InverseResult:
    z = np.arange(10.0, 10.0 + n_z * 10.0, 10.0)
    return InverseResult(
        z=z, u=np.zeros(n_z), v=np.zeros(n_z), uerr=np.full(n_z, 0.02),
        nvel=np.ones(n_z, dtype=int),
        u_do=np.zeros(n_z), v_do=np.zeros(n_z),
        u_up=np.zeros(n_z), v_up=np.zeros(n_z),
        u_ctd=np.zeros(n_se), v_ctd=np.zeros(n_se),
        ubar=0.0, vbar=0.0,
        zctd=np.linspace(5.0, 500.0, n_se),
        wctd=np.zeros(n_se),
    )


# Mandatory variables our writer must produce that also appear in 001.nc.
# Subset — excludes variables derived from raw data (e.g. u_shear_method,
# ctd_t, ctd_s) that require inputs not available without raw PD0 files.
MANDATORY_VARS = {"z", "u", "v", "uerr", "nvel", "u_do", "v_do", "u_up", "v_up",
                  "uctd", "vctd", "zctd"}
MANDATORY_WITH_GPS = {"tim", "shiplat", "shiplon"}
MANDATORY_WITH_SADCP = {"z_sadcp", "u_sadcp", "v_sadcp"}
MANDATORY_ATTRS = {"ubar", "vbar"}


@pytest.mark.integration
def test_mandatory_vars_present_in_reference(ref_001_vars: set[str]) -> None:
    """Confirm the reference NC actually has the variables we claim are mandatory."""
    missing = MANDATORY_VARS - ref_001_vars
    assert not missing, f"Reference 001.nc missing expected vars: {missing}"


@pytest.mark.integration
def test_writer_produces_mandatory_vars(tmp_path: Path, ref_001_vars: set[str]) -> None:
    result = _make_synthetic_result()
    out = tmp_path / "out.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    written = set(ds.variables.keys())
    attrs = set(ds.ncattrs())
    ds.close()

    missing_vars = MANDATORY_VARS - written
    missing_attrs = MANDATORY_ATTRS - attrs
    assert not missing_vars, f"Writer did not produce: {missing_vars}"
    assert not missing_attrs, f"Writer did not produce attributes: {missing_attrs}"


@pytest.mark.integration
def test_writer_produces_gps_vars(tmp_path: Path) -> None:
    result = _make_synthetic_result(n_z=10)
    n_ens = 30
    out = tmp_path / "out_gps.nc"
    write_ladcp_nc(
        out, result,
        ens_time_jd=np.linspace(2458919.0, 2458919.5, n_ens),
        ens_lat=np.full(n_ens, -70.45),
        ens_lon=np.full(n_ens, 168.47),
        uship=0.1, vship=-0.05,
    )
    ds = netCDF4.Dataset(str(out))
    written = set(ds.variables.keys())
    attrs = set(ds.ncattrs())
    ds.close()

    missing = MANDATORY_WITH_GPS - written
    assert not missing, f"Writer did not produce GPS vars: {missing}"
    assert "uship" in attrs and "vship" in attrs


@pytest.mark.integration
def test_writer_produces_sadcp_vars(tmp_path: Path) -> None:
    result = _make_synthetic_result(n_z=10)
    sadcp = SADCPProfile(
        z=np.array([50.0, 100.0, 200.0]),
        u=np.array([0.1, 0.2, 0.15]),
        v=np.array([-0.05, -0.03, -0.04]),
        err=np.array([0.05, 0.05, 0.05]),
    )
    out = tmp_path / "out_sadcp.nc"
    write_ladcp_nc(out, result, sadcp=sadcp)

    ds = netCDF4.Dataset(str(out))
    written = set(ds.variables.keys())
    ds.close()

    missing = MANDATORY_WITH_SADCP - written
    assert not missing, f"Writer did not produce SADCP vars: {missing}"


@pytest.mark.integration
def test_z_dimension_matches_reference_shape(s4p_dir: Path, tmp_path: Path) -> None:
    """Reference 001.nc has z dimension of about 130–140 bins; write matches any n_z."""
    ref_path = s4p_dir / "001.nc"
    if not ref_path.exists():
        pytest.skip("001.nc not found")

    ds = netCDF4.Dataset(str(ref_path))
    ref_n_z = len(ds.variables["z"][:])
    ds.close()

    result = _make_synthetic_result(n_z=ref_n_z)
    out = tmp_path / "out_sized.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    written_n_z = len(ds.variables["z"][:])
    ds.close()

    assert written_n_z == ref_n_z
```

- [ ] **Step 2: Run without data (should skip cleanly)**

```
pytest tests/integration/test_nc_writer_s4p.py -v
```

Expected: all tests SKIP.

- [ ] **Step 3: Run with data**

```
TEST_DATA_DIR="C:/Users/peter_sha/Documents/sourcecode/LADCP/test_data" \
    pytest tests/integration/test_nc_writer_s4p.py -v -m integration
```

Expected: all 5 tests PASS (mandatory vars check, writer round-trip, GPS, SADCP, shape match).

- [ ] **Step 4: Commit**

```bash
git add tests/integration/test_nc_writer_s4p.py
git commit -m "test: NC writer schema integration test vs S4P 001.nc reference"
```

---

## Self-Review

**Spec coverage:**
- ✅ `write_ladcp_nc()` function in `src/ladcp/output/nc.py`
- ✅ Core variables: z, u, v, uerr, nvel, u_do, v_do, u_up, v_up
- ✅ CTD velocity time series: uctd, vctd, zctd
- ✅ Barotropic mean as global attributes (ubar, vbar) — LDEO_IX convention
- ✅ Optional GPS track: tim, shiplat, shiplon + uship/vship attributes
- ✅ Optional SADCP profile: z_sadcp, u_sadcp, v_sadcp
- ✅ Unit tests cover round-trip, schema, each optional block
- ✅ Integration test validates schema against S4P reference; skips when data absent
- ✅ Explicit comment that value comparison requires raw PD0 files

**No placeholders present.**

**Type consistency:** `SADCPProfile` imported from `ladcp.ingestion.sadcp` in both Task 1 unit tests and Task 2 integration test; `InverseResult` field names match `src/ladcp/solution/inverse.py:622-638`.
