"""CTD time-series loading and ADCP bin depth assignment."""
from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from ladcp.ingestion._pd0 import _to_julian


@dataclass
class CTDTimeSeries:
    time_julian: NDArray[np.float64]    # (nctd,) Julian days
    pressure_dbar: NDArray[np.float64]  # (nctd,) positive down
    temp_c: NDArray[np.float64]         # (nctd,) NaN if absent
    salinity: NDArray[np.float64]       # (nctd,) NaN if absent


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
    }
    for line in lines:
        if line.strip() == '*END*':
            break
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
            result['start_time_julian'] = _parse_start_time(
                content.split('=', 1)[1].strip()
            )
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


def _build_ctd_time_series(
    arr: np.ndarray,
    col_roles: dict[int, str | None],
    time_start_julian: float | None,
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
    )


def _read_sbe_ascii(
    path: Path, data_offset: int, header_info: dict, **kwargs
) -> CTDTimeSeries:
    nquan = header_info['nquan']
    if nquan is None:
        raise ValueError(f"# nquan not found in header of {path}")
    col_roles = {i: _map_column(name) for i, name in header_info['columns'].items()}
    bad_flag = header_info['bad_flag']
    time_start_julian: float | None = (
        header_info.get('start_time_julian') or kwargs.get('time_start_julian')
    )

    arr = np.loadtxt(path, skiprows=0, comments=['*', '#']).reshape(-1, nquan)
    arr = arr.astype(np.float64)
    if bad_flag is not None:
        mask = np.isclose(arr, bad_flag, rtol=1e-3, atol=0)
        arr[mask] = np.nan

    return _build_ctd_time_series(arr, col_roles, time_start_julian)


def _read_sbe_binary(
    path: Path, data_offset: int, header_info: dict, **kwargs
) -> CTDTimeSeries:
    nquan = header_info['nquan']
    if nquan is None:
        raise ValueError(f"# nquan not found in header of {path}")
    col_roles = {i: _map_column(name) for i, name in header_info['columns'].items()}
    bad_flag = header_info['bad_flag']
    time_start_julian: float | None = (
        header_info.get('start_time_julian') or kwargs.get('time_start_julian')
    )

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

    return _build_ctd_time_series(arr, col_roles, time_start_julian)


def load_ctd(path: Path | str, **kwargs) -> CTDTimeSeries:
    """Load a CTD time-series file; auto-detect SBE binary, SBE ASCII, or generic ASCII."""
    path = Path(path)
    lines, data_offset = _read_header_lines(path)
    header_info = _parse_sbe_header(lines)
    fmt = _detect_format(header_info)

    if fmt == 'binary':
        return _read_sbe_binary(path, data_offset, header_info, **kwargs)
    elif fmt == 'sbe_ascii':
        return _read_sbe_ascii(path, data_offset, header_info, **kwargs)
    else:
        raise NotImplementedError("Generic ASCII reader not yet implemented")
