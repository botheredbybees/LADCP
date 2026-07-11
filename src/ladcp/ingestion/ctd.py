"""CTD time-series loading and ADCP bin depth assignment."""
from __future__ import annotations

import math
import re
import warnings
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
        'nvalues': None,
        'columns': {},
        'bad_flag': None,
        'file_type': 'ascii',
        'start_time_julian': None,
        'start_year': None,
        'lat_lon_appended': False,
        'interval_s': None,
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
        elif content.startswith('nvalues ='):
            result['nvalues'] = int(content.split('=', 1)[1].strip())
        elif m := re.match(r'name (\d+) = (\w+)', content):
            result['columns'][int(m.group(1))] = m.group(2)
        elif content.startswith('bad_flag ='):
            result['bad_flag'] = float(content.split('=', 1)[1].strip())
        elif content.startswith('interval ='):
            # "# interval = seconds: 0.0416667" -- scan interval for files
            # with no explicit time column (e.g. I7N 2018 raw 24 Hz cnv).
            m_int = re.match(r'interval\s*=\s*seconds:\s*([0-9.eE+-]+)', content)
            if m_int:
                result['interval_s'] = float(m_int.group(1))
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
    interval_s: float | None = None,
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
        if interval_s is not None and time_start_julian is not None:
            # No time column, but the header gives a constant scan interval
            # and a start time (e.g. I7N 2018 raw 24 Hz cnv: "# interval =
            # seconds: 0.0416667" + "# start_time = ... [NMEA time,
            # header]"). Synthesize the axis; a constant offset in the NMEA
            # header time relative to the ADCP clock is later measured and
            # corrected by estimate_ctd_adcp_lag() (+/-75 s search window
            # with a min-correlation tripwire), so start_time only needs to
            # be accurate to better than that window. Assumes a gap-free
            # constant-rate stream, the same assumption LDEO_IX's own
            # 24Hz->1Hz decimation makes for these files.
            n_scans = arr.shape[0]
            time_julian = time_start_julian + np.arange(n_scans) * interval_s / 86400.0
            return CTDTimeSeries(
                time_julian=time_julian,
                pressure_dbar=pressure,
                temp_c=temp if temp is not None else np.full(n_scans, np.nan),
                salinity=salinity if salinity is not None else np.full(n_scans, np.nan),
                lat=lat,
                lon=lon,
            )
        raise ValueError(
            "No time column found (expected: timeJ or timeS), and no "
            "'# interval = seconds' + '# start_time' header pair to "
            "synthesize a time axis from"
        )

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
        interval_s=header_info.get('interval_s'),
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
        interval_s=header_info.get('interval_s'),
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


def _uh_timeseries_columns(path: Path) -> dict[str, int] | None:
    """Detect the UH/SOEST CLIVAR-archive CTD time-series format.

    These files (NCEI raw-level0 'ctd_timeseries_SSSSS_CCC_gps.txt', e.g.
    A16N 2013 accession 0205839) have two '#' comment lines -- the first
    naming the columns (seconds pressure temperature salinity ...
    timestamp latitude longitude), the second giving units -- followed by
    plain numeric rows. Returns {name: column index} when the first line
    matches, else None.
    """
    with open(path, 'rb') as fh:
        first = fh.readline(4096).decode('ascii', errors='replace')
    if not first.startswith('#'):
        return None
    names = [t.lower() for t in first.lstrip('#').split()]
    col = {n: i for i, n in enumerate(names)}
    if 'pressure' in col and ('timestamp' in col or 'seconds' in col):
        return col
    return None


def _read_uh_timeseries(
    path: Path,
    col: dict[str, int],
    *,
    dday_yearbase: int = 2000,
    time_start_julian: float | None = None,
    **_ignored,
) -> CTDTimeSeries:
    """Read a UH/SOEST CLIVAR-archive CTD time series (see _uh_timeseries_columns).

    Absolute time comes from the 'timestamp' column -- GPS decimal days
    (dday) since dday_yearbase Jan 1 00:00 (A16N 2013 uses yearbase
    2000). GPS timestamps update in ~1 s steps while scans are 0.5 s, so
    the smooth elapsed 'seconds' column carries the per-scan spacing and
    the GPS dday only anchors it (median offset). Without a timestamp
    column, time_start_julian is required to anchor the elapsed axis.
    """
    arr = np.loadtxt(path, comments='#')
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    n = arr.shape[0]

    def _get(name: str) -> np.ndarray | None:
        i = col.get(name)
        return arr[:, i] if i is not None and i < arr.shape[1] else None

    pressure = _get('pressure')
    if pressure is None:
        raise ValueError(f"No pressure column in {path}")
    elapsed = _get('seconds')
    dday = _get('timestamp')
    if dday is not None:
        base = _to_julian(dday_yearbase, 1, 1, 0.0)
        if elapsed is not None:
            anchor = base + float(np.nanmedian(dday - elapsed / 86400.0))
            time_julian = anchor + elapsed / 86400.0
        else:
            time_julian = base + dday
    elif elapsed is not None and time_start_julian is not None:
        time_julian = time_start_julian + elapsed / 86400.0
    else:
        raise ValueError(
            f"{path}: no 'timestamp' (dday) column and no time_start_julian "
            "anchor for the elapsed 'seconds' column"
        )

    temp = _get('temperature')
    sal = _get('salinity')
    return CTDTimeSeries(
        time_julian=time_julian,
        pressure_dbar=pressure,
        temp_c=temp if temp is not None else np.full(n, np.nan),
        salinity=sal if sal is not None else np.full(n, np.nan),
        lat=_get('latitude'),
        lon=_get('longitude'),
    )


def _detect_binary_by_filesize(path: Path, data_offset: int, header_info: dict) -> None:
    """Override file_type to 'binary' when the data section size matches nvalues×nquan×4.

    Also clears lat_lon_appended when the size matches nquan (not nquan+2) floats per scan,
    which happens when Seasave writes binary without embedding the appended lat/lon columns.
    """
    nquan = header_info.get('nquan')
    nvalues = header_info.get('nvalues')
    if not nquan or not nvalues:
        return
    data_bytes = path.stat().st_size - data_offset
    if abs(data_bytes - nvalues * nquan * 4) <= 4:
        header_info['file_type'] = 'binary'
        header_info['lat_lon_appended'] = False
    elif abs(data_bytes - nvalues * (nquan + 2) * 4) <= 4:
        header_info['file_type'] = 'binary'


def load_ctd(path: Path | str, **kwargs) -> CTDTimeSeries:
    """Load a CTD time-series file; auto-detects SBE binary, SBE ASCII,
    UH/SOEST CLIVAR-archive time series (named-column '#' header, GPS dday
    timestamps -- see _read_uh_timeseries), or generic ASCII.

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
        _detect_binary_by_filesize(path, data_offset, header_info)
        fmt = _detect_format(header_info)

        if fmt == 'binary':
            return _read_sbe_binary(path, data_offset, header_info, **kwargs)
        elif fmt == 'sbe_ascii':
            return _read_sbe_ascii(path, data_offset, header_info, **kwargs)
        else:
            return _read_generic_ascii(path, **kwargs)
    except ValueError as e:
        if "No *END* marker found" in str(e):
            uh_cols = _uh_timeseries_columns(path)
            if uh_cols is not None:
                return _read_uh_timeseries(path, uh_cols, **kwargs)
            return _read_generic_ascii(path, **kwargs)
        raise


def pressure_to_depth(
    p_dbar: NDArray[np.float64], lat_deg: float
) -> NDArray[np.float64]:
    """Pressure (dbar) to depth (m), Saunders & Fofonoff (1976), EOS-80 refit.

    Port of loadctd.m::p2z (LDEO_IX). Check value: 9712.654 m at p = 10000
    dbar, lat = 30. Uses z=0 at p=0 (no 1-atm surface offset).
    """
    p = np.asarray(p_dbar, dtype=np.float64) / 10.0  # dbar -> bars, as in p2z
    x = math.sin(math.radians(lat_deg)) ** 2
    g = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-5 * p
    depth = (((-1.82e-11 * p + 2.279e-7) * p - 2.2512e-3) * p + 97.2659) * p
    return depth / g


def assign_bin_depths(
    rdi: RDIData,
    ctd: CTDTimeSeries,
    *,
    looker: str = "down",
    lat_deg: float | None = None,
    time_offset_days: float = 0.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute instrument depth and bin absolute depths from CTD and ADCP geometry.

    lat_deg selects the Saunders & Fofonoff pressure->depth conversion
    (pressure_to_depth(), matching LDEO_IX loadctd.m); without it a crude
    shallow-water fallback z = p * 1.00445 is used, which reads ~2% too deep
    below ~2000 dbar -- always pass lat_deg when it is known.

    time_offset_days shifts the CTD pressure sampling to ADCP time + offset,
    correcting a CTD-ADCP clock offset (see estimate_ctd_adcp_lag()).

    Returns:
        z_m: (nens,) instrument depth in metres, positive down
        izm: (nbin, nens) absolute depth of each bin in metres, positive down
    """
    p_interp = np.interp(
        rdi.time_julian + time_offset_days, ctd.time_julian, ctd.pressure_dbar
    )

    if lat_deg is not None:
        z_m = pressure_to_depth(p_interp, lat_deg)
    else:
        z_m = p_interp * 1.00445

    bin_offsets = rdi.dist_m + np.arange(rdi.nbin) * rdi.blen_m  # (nbin,)
    sign = 1.0 if looker == "down" else -1.0
    izm = z_m[np.newaxis, :] + sign * bin_offsets[:, np.newaxis]  # (nbin, nens)

    return z_m, izm


def ctd_in_water_window(
    ctd: CTDTimeSeries,
    *,
    cut_dbar: float = 10.0,
) -> tuple[float, float]:
    """CTD in-water time window (loadctd.m 'Cut CTD profile', p.cut).

    Port of loadctd.m lines 236-276: the cast start is the LAST scan
    before the pressure maximum with 0 < p < cut_dbar (package just below
    the surface starting its descent), the end is the FIRST scan after
    the maximum with 0 < p < cut_dbar (package back at the surface).
    On-deck scans (p <= 0 or NaN) never qualify. LDEO then discards all
    ADCP ensembles outside this window (loadctd.m:517-520) -- without
    the cut, pre-deployment and post-recovery ensembles enter the
    solution with clamped surface pressure (I7N 003: ~7 min of on-deck
    pinging before the CTD record even starts).

    Returns:
        (t_start_julian, t_end_julian) in the CTD time base. Falls back
        to the full record boundary on either side where no qualifying
        scan exists (same as loadctd.m's istart=1 / iend=end branches).
    """
    p = ctd.pressure_dbar
    t = ctd.time_julian
    finite = np.isfinite(p) & np.isfinite(t)
    if not finite.any():
        raise ValueError("CTD record has no finite pressure/time scans")
    ipmax = int(np.nanargmax(np.where(finite, p, -np.inf)))
    near_surface = finite & (np.abs(p - cut_dbar / 2.0) < cut_dbar / 2.0)

    before = np.flatnonzero(near_surface[: ipmax + 1])
    istart = int(before[-1]) if before.size else 0
    after = np.flatnonzero(near_surface[ipmax:])
    iend = int(after[0]) + ipmax if after.size else len(t) - 1
    return float(t[istart]), float(t[iend])


def estimate_ctd_adcp_lag(
    time_adcp_julian: NDArray[np.float64],
    w_adcp: NDArray[np.float64],
    ctd: CTDTimeSeries,
    *,
    lat_deg: float | None = None,
    max_lag_scans: int = 150,
    coarse_window_s: float = 600.0,
) -> tuple[int, float, float]:
    """Estimate the CTD-ADCP clock offset by w cross-correlation.

    Equivalent of LDEO_IX loadctd.m's besttlag() step (whole-series variant):
    the package sinking rate is observed independently by the CTD (as the
    pressure time-derivative) and by the ADCP (earth-frame vertical velocity),
    so cross-correlating the two time series exposes any clock offset.

    Args:
        time_adcp_julian: (nens,) ADCP ensemble times (Julian days).
        w_adcp: (nens,) earth-frame vertical velocity per ensemble, e.g.
            nanmedian over bins (positive during descent).
        ctd: CTD time series (its own clock).
        lat_deg: latitude for pressure_to_depth; fallback 1.00445 when None.
        max_lag_scans: search window in resampled CTD scans (LDEO ctdmaxlag
            default 150).
        coarse_window_s: when > max_lag_scans * scan interval, a coarse
            pre-alignment pass first scans +/- this many seconds to find
            the neighbourhood of the correlation peak, and the fine
            +/-max_lag_scans search runs centred there. Needed when the
            CTD time base is a synthesized NMEA-header start time
            (I7N-style cnv with no time column), where the deck-box vs
            ADCP clock offset (~3 min measured on I7N 2018) far exceeds
            the LDEO-default fine window. 0 disables the coarse pass.

    Returns:
        (lag_scans, lagdt_days, corr): CTD pressure evaluated at
        time_adcp_julian + lagdt_days aligns with the ADCP samples --
        i.e. pass lagdt_days as assign_bin_depths(time_offset_days=...).
        lag_scans counts resampled-grid intervals of lagdt_days.
    """
    # Resample pressure onto a regular >=0.5 s grid before differentiating:
    # high-rate SBE files quantize the time column coarsely (24 Hz data with
    # ~0.7 s timestamp resolution -> most consecutive dt are exactly 0), so a
    # per-scan gradient is undefined. LDEO's workflow has the same property
    # via its 2 Hz files + loadctd.m's time-jitter fix.
    finite_p = np.isfinite(ctd.pressure_dbar) & np.isfinite(ctd.time_julian)
    t_p = ctd.time_julian[finite_p]
    span_days = float(t_p[-1] - t_p[0])
    nominal_dt_s = span_days * 86400.0 / max(len(t_p) - 1, 1)
    dtctd_s = max(0.5, nominal_dt_s)
    dtctd_days = dtctd_s / 86400.0
    t_grid = t_p[0] + np.arange(int(span_days / dtctd_days) + 1) * dtctd_days
    p_grid = np.interp(t_grid, t_p, ctd.pressure_dbar[finite_p])

    if lat_deg is not None:
        z_ctd = -pressure_to_depth(p_grid, lat_deg)
    else:
        z_ctd = -p_grid * 1.00445
    # w_ctd positive during descent (z negative-down, decreasing), loadctd.m
    w_ctd = -np.gradient(z_ctd, dtctd_s)

    finite_w = np.isfinite(w_adcp) & np.isfinite(time_adcp_julian)
    t_a = time_adcp_julian[finite_w]
    w_a = w_adcp[finite_w]
    w_a = w_a - np.nanmean(w_a)

    t_c = t_grid
    w_c = w_ctd

    def _scan(lags: range) -> tuple[int, float]:
        best = (lags[0], -np.inf)
        for lag in lags:
            wc = np.interp(t_a + lag * dtctd_days, t_c, w_c)
            wc = wc - wc.mean()
            denom = math.sqrt(float(np.dot(wc, wc)) * float(np.dot(w_a, w_a)))
            if denom == 0.0:
                continue
            c = float(np.dot(wc, w_a)) / denom
            if c > best[1]:
                best = (lag, c)
        return best

    # Coarse pre-alignment: the heave-driven correlation peak is sharp
    # (falls off within ~1 s on 2018 I7N data), so the coarse step must
    # stay near 1 s regardless of the resampled scan interval.
    center = 0
    coarse_max = int(coarse_window_s / dtctd_s)
    if coarse_max > max_lag_scans:
        step = max(1, int(round(1.0 / dtctd_s)))
        center, c_corr = _scan(range(-coarse_max, coarse_max + 1, step))
        if abs(center) + step > coarse_max:
            warnings.warn(
                f"CTD-ADCP coarse lag search hit the +/-{coarse_window_s:.0f}-s"
                f" window edge (lag={center * dtctd_s:.1f} s, "
                f"corr={c_corr:.3f}); the CTD time base may be wrong by more "
                "than that",
                stacklevel=2,
            )

    lag_scans, corr = _scan(
        range(center - max_lag_scans, center + max_lag_scans + 1))
    if abs(lag_scans - center) >= max_lag_scans:
        warnings.warn(
            f"CTD-ADCP lag search hit the +/-{max_lag_scans}-scan window edge "
            f"(lag={lag_scans}, corr={corr:.3f}); a coarse pre-alignment may "
            "be needed",
            stacklevel=2,
        )
    return lag_scans, lag_scans * dtctd_days, corr


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
