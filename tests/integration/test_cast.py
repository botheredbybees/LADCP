"""Integration tests against I7N GO-SHIP cast 002 reference data."""

import pytest
import xarray as xr


@pytest.mark.integration
def test_reference_netcdf_loads(test_data_dir):
    """Reference LDEO_IX output for cast 002 is readable and has data."""
    nc_path = test_data_dir / "data" / "002.nc"
    assert nc_path.exists(), f"Reference NetCDF not found at {nc_path}"

    ds = xr.open_dataset(nc_path)
    assert len(ds.data_vars) > 0, "NetCDF has no data variables"
    ds.close()


@pytest.mark.integration
def test_reference_netcdf_has_velocity(test_data_dir):
    """Reference NetCDF contains at least one velocity variable."""
    ds = xr.open_dataset(test_data_dir / "data" / "002.nc")
    vel_vars = [
        v
        for v in ds.data_vars
        if any(tok in v.lower() for tok in ("u", "v", "vel", "east", "north"))
    ]
    assert vel_vars, (
        f"No velocity variable found. Variables present: {list(ds.data_vars)}"
    )
    ds.close()


@pytest.mark.integration
def test_reference_cast_depth(test_data_dir):
    """Cast 002 reaches approximately 4892 m (within 50 m tolerance)."""
    ds = xr.open_dataset(test_data_dir / "data" / "002.nc")
    # depth coordinate may be named 'depth', 'z', 'pressure', or similar
    depth_coords = [
        c
        for c in ds.coords
        if any(tok in c.lower() for tok in ("depth", "z", "pressure", "p"))
    ]
    assert depth_coords, f"No depth coordinate found. Coords: {list(ds.coords)}"
    import numpy as np

    max_depth = float(np.nanmax(ds[depth_coords[0]].values))
    assert abs(max_depth - 4892) < 50, (
        f"Expected max depth ~4892 m, got {max_depth:.1f} m"
    )
    ds.close()
