"""Unit tests for load_sadcp_nc()."""
from __future__ import annotations
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion._pd0 import _to_julian
from ladcp.ingestion.sadcp import SADCPProfile, load_sadcp_nc

# Julian day of 2020-01-01 00:00 UTC in the midnight-based _to_julian convention
# (matches the dialect used by the PD0 parser and _parse_epoch_jd)
JAN1_2020_JD = _to_julian(2020, 1, 1, 0.0)  # 2458850.0


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
