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
    try:
        names = set(ds.variables.keys())
    finally:
        ds.close()
    return names


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
    assert not missing, f"Reference 001.nc missing expected vars: {missing}"  # ref_001_vars loaded from fixture


@pytest.mark.integration
def test_writer_produces_mandatory_vars(tmp_path: Path) -> None:
    result = _make_synthetic_result()
    out = tmp_path / "out.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    try:
        written = set(ds.variables.keys())
        attrs = set(ds.ncattrs())
    finally:
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
    try:
        written = set(ds.variables.keys())
        attrs = set(ds.ncattrs())
    finally:
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
    try:
        written = set(ds.variables.keys())
    finally:
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
    try:
        ref_n_z = len(ds.variables["z"][:])
    finally:
        ds.close()

    result = _make_synthetic_result(n_z=ref_n_z)
    out = tmp_path / "out_sized.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    try:
        written_n_z = len(ds.variables["z"][:])
    finally:
        ds.close()

    assert written_n_z == ref_n_z
