"""Integration test: SBE hex decoder vs S4P cast 001 CTD data.

Requires TEST_DATA_DIR env var pointing to a directory containing 2018_S4P/.
Defaults to 'test_data' (the repo's bundled data directory).

Run:
    uv run pytest tests/integration/test_sbe_hex_s4p.py -v -m integration
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from ladcp.ingestion.sbe_hex import load_sbe_hex


@pytest.fixture(scope="module")
def s4p_dir() -> Path:
    env = os.environ.get("TEST_DATA_DIR", "test_data")
    p = Path(env) / "2018_S4P"
    if not p.exists():
        pytest.skip(f"2018_S4P not found at {p}")
    return p


@pytest.fixture(scope="module")
def hex_path(s4p_dir: Path) -> Path:
    p = s4p_dir / "CTD/00101.hex"
    if not p.exists():
        pytest.skip(f"00101.hex not found: {p}")
    return p


@pytest.fixture(scope="module")
def xmlcon_path(s4p_dir: Path) -> Path:
    p = s4p_dir / "CTD/00101.XMLCON"
    if not p.exists():
        pytest.skip(f"00101.XMLCON not found: {p}")
    return p


@pytest.fixture(scope="module")
def decoded(hex_path: Path, xmlcon_path: Path):
    return load_sbe_hex(hex_path, xmlcon_path)


@pytest.mark.integration
def test_decoded_has_scans(decoded) -> None:
    assert len(decoded.time_julian) > 1000, "Expected >1000 scans"


@pytest.mark.integration
def test_decoded_pressure_reaches_bottom(decoded) -> None:
    """Cast 001 reaches ~1361 m; decoded pressure should exceed 1000 dbar."""
    assert np.nanmax(decoded.pressure_dbar) > 1000.0


@pytest.mark.integration
def test_decoded_temperature_range(decoded) -> None:
    """Southern Ocean temps: -2°C to +5°C in the water column; warm on-deck allowed."""
    T_valid = decoded.temp_c[np.isfinite(decoded.temp_c)]
    assert len(T_valid) > 0
    assert T_valid.min() > -3.0
    assert T_valid.max() < 35.0


@pytest.mark.integration
def test_gps_lat_lon_present(decoded) -> None:
    """GPS lat/lon must be present with many valid readings."""
    assert decoded.lat is not None, "No GPS lat decoded"
    assert decoded.lon is not None, "No GPS lon decoded"
    lat_valid = decoded.lat[np.isfinite(decoded.lat)]
    assert len(lat_valid) > 100, "Too few valid GPS lat values"


@pytest.mark.integration
def test_gps_position_near_expected(decoded) -> None:
    """Median GPS position should be near cast 001 header position (−70.45°S, 168.47°E)."""
    lat_valid = decoded.lat[np.isfinite(decoded.lat)]
    lon_valid = decoded.lon[np.isfinite(decoded.lon)]
    assert abs(np.nanmedian(lat_valid) - (-70.45)) < 0.5
    assert abs(np.nanmedian(lon_valid) - 168.47) < 0.5


@pytest.mark.integration
def test_gps_position_tight(decoded) -> None:
    """GPS byte offsets confirmed: mean position within 0.1° of header NMEA fix."""
    lat_valid = decoded.lat[np.isfinite(decoded.lat)]
    lon_valid = decoded.lon[np.isfinite(decoded.lon)]
    assert abs(np.nanmedian(lat_valid) - (-70.45)) < 0.1
    assert abs(np.nanmedian(lon_valid) - 168.47) < 0.1
