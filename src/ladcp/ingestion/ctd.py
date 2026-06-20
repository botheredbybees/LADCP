"""CTD time-series loading and ADCP bin depth assignment."""
from __future__ import annotations

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
    dt = datetime.strptime(date_str.strip(), '%b %d %Y %H:%M:%S')
    frac_hour = dt.hour + dt.minute / 60.0 + dt.second / 3600.0
    return _to_julian(dt.year, dt.month, dt.day, frac_hour)
