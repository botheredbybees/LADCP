"""Unit tests for write_ladcp_nc()."""
from __future__ import annotations
from pathlib import Path

import netCDF4
import numpy as np

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
    try:
        for var in ("z", "u", "v", "uerr", "nvel", "u_do", "v_do", "u_up", "v_up"):
            assert var in ds.variables, f"Missing variable: {var}"
            assert ds.variables[var].shape == (10,), f"{var}: expected shape (10,)"
    finally:
        ds.close()


def test_write_ctd_velocity_variables(tmp_path: Path) -> None:
    result = _make_result(n_z=5, n_se=8)
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result)

    ds = netCDF4.Dataset(str(out))
    try:
        for var in ("uctd", "vctd", "zctd"):
            assert var in ds.variables, f"Missing variable: {var}"
            assert ds.variables[var].shape == (8,)
    finally:
        ds.close()


def test_write_barotropic_as_attributes(tmp_path: Path) -> None:
    result = _make_result()
    result_with_ubar = InverseResult(
        **{**result.__dict__, "ubar": 0.123, "vbar": -0.045}
    )
    out = tmp_path / "test.nc"
    write_ladcp_nc(out, result_with_ubar)

    ds = netCDF4.Dataset(str(out))
    try:
        assert hasattr(ds, "ubar"), "ubar should be a global attribute"
        assert hasattr(ds, "vbar"), "vbar should be a global attribute"
        assert abs(float(ds.ubar) - 0.123) < 1e-6
        assert abs(float(ds.vbar) - -0.045) < 1e-6
    finally:
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
    try:
        assert "tim" in ds.variables
        assert "shiplat" in ds.variables
        assert "shiplon" in ds.variables
        assert ds.variables["tim"].shape == (n_ens,)
        assert hasattr(ds, "uship")
        assert hasattr(ds, "vship")
    finally:
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
    try:
        assert "u_sadcp" in ds.variables
        assert "v_sadcp" in ds.variables
        assert "z_sadcp" in ds.variables
        assert ds.variables["z_sadcp"].shape == (3,)
    finally:
        ds.close()
