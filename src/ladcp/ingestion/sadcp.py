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
        # Floor before nvel scaling so zero-std bins still get inflated for sparse records
        v_err = np.where(v_err == 0.0, err_default, v_err)
        max_nvel = float(np.nanmax(nvel)) if np.any(nvel > 0) else 1.0
        # Scale: bins with fewer records get proportionally larger error
        with np.errstate(divide="ignore", invalid="ignore"):
            v_err = np.where(nvel > 0, v_err * max_nvel / nvel, err_default)
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
