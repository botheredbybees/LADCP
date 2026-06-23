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

    _gps = (ens_time_jd, ens_lat, ens_lon)
    if any(x is not None for x in _gps) and not all(x is not None for x in _gps):
        raise ValueError(
            "ens_time_jd, ens_lat, and ens_lon must all be provided together or all omitted"
        )

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
            tim_v.units = "Julian days"
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
