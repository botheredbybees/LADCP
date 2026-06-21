"""CTD time-series loading and ADCP bin depth assignment."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from ladcp.ingestion._pd0 import _to_julian
from ladcp.ingestion._types import RDIData


@dataclass
class CTDTimeSeries:
    time_julian: NDArray[np.float64]    # (nctd,) Julian days
    pressure_dbar: NDArray[np.float64]  # (nctd,) positive down
    temp_c: NDArray[np.float64]         # (nctd,) NaN if absent
    salinity: NDArray[np.float64]       # (nctd,) NaN if absent
    lat: NDArray[np.float64] | None = None  # (nctd,) degrees N; NaN = bad fix
    lon: NDArray[np.float64] | None = None  # (nctd,) degrees E; NaN = bad fix


def _read_header_lines(path: Path) -> tuple[list[str], int]:
    """Open file in binary mode, decode header until *END*."""
    lines: list[str] = []
    with open(path, 'rb') as f:
        while True:
            raw = f.readline()
            if not raw:
                raise ValueError(f"No *END* marker found in {path}")
            line = raw.decode('latin-1', errors='replace')
            lines.append(line)
            if line.strip() == '*END*':
                return lines, f.tell()


def _parse_sbe_header(lines: list[str]) -> dict:
    """Extract structured metadata from SBE CNV header lines."""
    result: dict = {
        'nquan': None,
        'columns': {},
        'bad_flag': None,
        'file_type': 'ascii',
        'start_time_julian': None,
        'start_year': None,
        'lat_lon_appended': False,
    }
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
    return result


def _detect_format(header_info: dict) -> str:
    if header_info['file_type'] == 'binary':
        return 'binary'
    if header_info['columns']:
        return 'sbe_ascii'
    return 'generic'


_COL_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r'^pr'), 'pressure'),
    (re.compile(r'^t[0-9]'), 'temp'),
    (re.compile(r'^sal'), 'salinity'),
    (re.compile(r'^timeJ$'), 'time_julian'),
    (re.compile(r'^timeS$'), 'time_elapsed_s'),
]


def _map_column(name: str) -> str | None:
    for pattern, role in _COL_PATTERNS:
        if pattern.match(name):
            return role
    return None


def _parse_start_time(date_str: str) -> float:
    dt = datetime.strptime(date_str.strip()[:20], '%b %d %Y %H:%M:%S')
    frac_hour = dt.hour + dt.minute / 60.0 + dt.second / 3600.0
    return _to_julian(dt.year, dt.month, dt.day, frac_hour)


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


def _build_ctd_time_series(
    arr: np.ndarray,
    col_roles: dict[int, str | None],
    time_start_julian: float | None,
    start_year: int | None = None,
    lat: np.ndarray | None = None,
    lon: np.ndarray | None = None,
) -> CTDTimeSeries:
    """Assemble CTDTimeSeries from a (n_scans, nquan) float64 array."""
    time_raw: np.ndarray | None = None
    time_role: str | None = None
    pressure: np.ndarray | None = None
    temp: np.ndarray | None = None
    salinity: np.ndarray | None = None

    for col_idx, role in col_roles.items():
        col = arr[:, col_idx]
        if role == 'time_julian' and time_raw is None:
            time_raw, time_role = col, 'julian'
        elif role == 'time_elapsed_s' and time_raw is None:
            time_raw, time_role = col, 'elapsed_s'
        elif role == 'pressure' and pressure is None:
            pressure = col
        elif role == 'temp' and temp is None:
            temp = col
        elif role == 'salinity' and salinity is None:
            salinity = col

    if pressure is None:
        raise ValueError("No pressure column found (expected prefix: prDM, prSM, prE)")
    if time_raw is None:
        raise ValueError("No time column found (expected: timeJ or timeS)")

    if time_role == 'julian':
        if start_year is not None and len(time_raw) > 0 and time_raw[0] <= 366:
            # timeJ is day-of-year (1 = Jan 1 00:00); convert to absolute Julian day
            jan1_julian = _to_julian(start_year, 1, 1, 0.0)
            time_julian = jan1_julian + (time_raw - 1.0)
        else:
            time_julian = time_raw
    else:
        if time_start_julian is None:
            raise ValueError(
                "timeS column requires time_start_julian (from header # start_time or kwarg)"
            )
        time_julian = time_raw / 86400.0 + time_start_julian

    n = len(pressure)
    return CTDTimeSeries(
        time_julian=time_julian,
        pressure_dbar=pressure,
        temp_c=temp if temp is not None else np.full(n, np.nan),
        salinity=salinity if salinity is not None else np.full(n, np.nan),
        lat=lat,
        lon=lon,
    )


def _read_sbe_ascii(
    path: Path, data_offset: int, header_info: dict, **kwargs
) -> CTDTimeSeries:
    nquan = header_info['nquan']
    if nquan is None:
        raise ValueError(f"# nquan not found in header of {path}")
    col_roles = {i: _map_column(name) for i, name in header_info['columns'].items()}
    bad_flag = header_info['bad_flag']
    _hdr_t = header_info.get('start_time_julian')
    time_start_julian: float | None = _hdr_t if _hdr_t is not None else kwargs.get('time_start_julian')

    # data_offset is unused here; np.loadtxt skips the header via comments=['*','#']
    arr = np.loadtxt(path, skiprows=0, comments=['*', '#']).reshape(-1, nquan)
    arr = arr.astype(np.float64)
    if bad_flag is not None:
        mask = np.isclose(arr, bad_flag, rtol=1e-3, atol=0)
        arr[mask] = np.nan

    lat = lon = None
    if header_info.get('lat_lon_appended') and arr.shape[1] >= 2:
        lat, lon, arr, col_roles = _extract_latlon(arr, col_roles)

    return _build_ctd_time_series(
        arr, col_roles, time_start_julian,
        start_year=header_info.get('start_year'),
        lat=lat, lon=lon,
    )


def _read_sbe_binary(
    path: Path, data_offset: int, header_info: dict, **kwargs
) -> CTDTimeSeries:
    nquan = header_info['nquan']
    if nquan is None:
        raise ValueError(f"# nquan not found in header of {path}")
    col_roles = {i: _map_column(name) for i, name in header_info['columns'].items()}
    bad_flag = header_info['bad_flag']
    _hdr_t = header_info.get('start_time_julian')
    time_start_julian: float | None = _hdr_t if _hdr_t is not None else kwargs.get('time_start_julian')

    with open(path, 'rb') as f:
        f.seek(data_offset)
        raw = f.read()

    n_scans = len(raw) // (nquan * 4)
    arr = (
        np.frombuffer(raw[: n_scans * nquan * 4], dtype='<f4')
        .reshape(n_scans, nquan)
        .astype(np.float64)
    )
    if bad_flag is not None:
        mask = np.isclose(arr, bad_flag, rtol=1e-3, atol=0)
        arr[mask] = np.nan

    lat = lon = None
    if header_info.get('lat_lon_appended') and arr.shape[1] >= 2:
        lat, lon, arr, col_roles = _extract_latlon(arr, col_roles)

    return _build_ctd_time_series(
        arr, col_roles, time_start_julian,
        start_year=header_info.get('start_year'),
        lat=lat, lon=lon,
    )


def _read_generic_ascii(path: Path, **kwargs) -> CTDTimeSeries:
    """Read generic fixed-column ASCII CTD file with configurable column layout.

    kwargs:
        skip_rows: int = 0
        col_time: int = 0
        col_pressure: int = 1
        col_temp: int = 2
        col_salinity: int | None = None
        time_base: str = 'elapsed_s'  ('elapsed_s' or 'julian')
        time_start_julian: float | None = None (required if time_base='elapsed_s')
    """
    skip_rows: int = kwargs.get('skip_rows', 0)
    col_time: int = kwargs.get('col_time', 0)
    col_pressure: int = kwargs.get('col_pressure', 1)
    col_temp: int = kwargs.get('col_temp', 2)
    col_salinity: int | None = kwargs.get('col_salinity', None)
    time_base: str = kwargs.get('time_base', 'elapsed_s')
    time_start_julian: float | None = kwargs.get('time_start_julian')

    arr = np.loadtxt(path, skiprows=skip_rows)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    time_raw = arr[:, col_time]
    pressure = arr[:, col_pressure]
    n = len(pressure)

    if time_base == 'julian':
        time_julian = time_raw
    elif time_base == 'elapsed_s':
        if time_start_julian is None:
            raise ValueError("time_base='elapsed_s' requires time_start_julian kwarg")
        time_julian = time_raw / 86400.0 + time_start_julian
    else:
        raise ValueError(f"Unknown time_base: {time_base!r}; expected 'elapsed_s' or 'julian'")

    temp = arr[:, col_temp] if col_temp < arr.shape[1] else np.full(n, np.nan)
    if col_salinity is not None and col_salinity < arr.shape[1]:
        salinity = arr[:, col_salinity]
    else:
        salinity = np.full(n, np.nan)

    return CTDTimeSeries(
        time_julian=time_julian,
        pressure_dbar=pressure,
        temp_c=temp,
        salinity=salinity,
    )


def load_ctd(path: Path | str, **kwargs) -> CTDTimeSeries:
    """Load a CTD time-series file; auto-detect SBE binary, SBE ASCII, or generic ASCII.

    For generic ASCII files, kwargs configure the reader:
        skip_rows (int): header lines to skip (default 0)
        col_time, col_pressure, col_temp, col_salinity (int): 0-based column indices
        time_base ('elapsed_s' | 'julian'): how to interpret the time column
        time_start_julian (float): required when time_base='elapsed_s'

    For SBE files, time_start_julian can override the header-derived value.
    Note: ADCP ensembles outside the CTD time range receive the first/last CTD
    pressure via flat extrapolation in assign_bin_depths(); mask those ensembles
    in QA before trusting the interpolated depth.
    """
    path = Path(path)
    try:
        lines, data_offset = _read_header_lines(path)
        header_info = _parse_sbe_header(lines)
        fmt = _detect_format(header_info)

        if fmt == 'binary':
            return _read_sbe_binary(path, data_offset, header_info, **kwargs)
        elif fmt == 'sbe_ascii':
            return _read_sbe_ascii(path, data_offset, header_info, **kwargs)
        else:
            return _read_generic_ascii(path, **kwargs)
    except ValueError as e:
        if "No *END* marker found" in str(e):
            return _read_generic_ascii(path, **kwargs)
        raise


def assign_bin_depths(
    rdi: RDIData,
    ctd: CTDTimeSeries,
    *,
    looker: str = "down",
    lat_deg: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute instrument depth and bin absolute depths from CTD and ADCP geometry.

    Returns:
        z_m: (nens,) instrument depth in metres, positive down
        izm: (nbin, nens) absolute depth of each bin in metres, positive down
    """
    p_interp = np.interp(rdi.time_julian, ctd.time_julian, ctd.pressure_dbar)

    if lat_deg is not None:
        sin2 = math.sin(math.radians(lat_deg)) ** 2
        g = 9.780318 * (1 + 5.2788e-3 * sin2) + 1.092e-6 * p_interp
        z_m = (
            9.72659 * p_interp
            - 2.2512e-5 * p_interp**2
            + 2.279e-10 * p_interp**3
            - 1.82e-15 * p_interp**4
        ) / g
    else:
        z_m = p_interp * 1.00445

    bin_offsets = rdi.dist_m + np.arange(rdi.nbin) * rdi.blen_m  # (nbin,)
    sign = 1.0 if looker == "down" else -1.0
    izm = z_m[np.newaxis, :] + sign * bin_offsets[:, np.newaxis]  # (nbin, nens)

    return z_m, izm
