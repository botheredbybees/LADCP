# tests/integration/test_ctd_p16n_cast003.py
"""Integration tests for CTD loading against P16N cast 003.

Requires TEST_DATA_DIR env var pointing to a directory containing:
  2015_P16N/003_01.cnv   — CTD time-series (binary SBE format)
  2015_P16N/003/003DL000.000  — Downlooker PD0 binary
"""
import os
from pathlib import Path

import numpy as np
import pytest

from ladcp.ingestion.ctd import CTDTimeSeries, assign_bin_depths, load_ctd
from ladcp.ingestion.rdi import load_rdi


@pytest.fixture
def test_data_dir() -> Path:
    path = Path(os.environ.get("TEST_DATA_DIR", "test_data"))
    if not path.exists():
        pytest.skip("TEST_DATA_DIR not populated — see test_data/sources.md")
    return path


@pytest.fixture
def cnv_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003_01.cnv"
    if not p.exists():
        pytest.skip(f"CTD file not found: {p}")
    return p


@pytest.fixture
def dl_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003" / "003DL000.000"
    if not p.exists():
        pytest.skip(f"DL PD0 file not found: {p}")
    return p


@pytest.fixture
def ctd(cnv_path: Path) -> CTDTimeSeries:
    return load_ctd(cnv_path)


@pytest.mark.integration
def test_load_cnv_returns_ctd_time_series(ctd: CTDTimeSeries):
    assert isinstance(ctd, CTDTimeSeries)
    assert ctd.pressure_dbar.ndim == 1
    assert ctd.time_julian.ndim == 1

@pytest.mark.integration
def test_cnv_pressure_range(ctd: CTDTimeSeries):
    # Header states max = 4367.032 dbar
    assert np.nanmax(ctd.pressure_dbar) > 4000.0
    assert np.nanmax(ctd.pressure_dbar) < 5000.0
    # Allow small sensor offset at surface/deck (~-0.56 dbar observed end-of-cast)
    assert np.nanmin(ctd.pressure_dbar) >= -1.0

@pytest.mark.integration
def test_cnv_scan_count(ctd: CTDTimeSeries):
    # Header: # nvalues = 307198
    assert len(ctd.time_julian) == 307198

@pytest.mark.integration
def test_cnv_time_monotone(ctd: CTDTimeSeries):
    diffs = np.diff(ctd.time_julian)
    assert np.all(diffs >= 0), f"time not monotone; {(diffs < 0).sum()} violations"

@pytest.mark.integration
def test_cnv_temperature_plausible(ctd: CTDTimeSeries):
    finite = ctd.temp_c[np.isfinite(ctd.temp_c)]
    assert len(finite) > 0
    # Allow a tiny fraction of pre-deployment/startup scans with bad values
    bad_low = (finite < -2.0).sum()
    assert bad_low / len(finite) < 0.001, f"{bad_low} samples below -2°C"
    assert np.all(finite <= 32.0)

@pytest.mark.integration
def test_assign_depths_shape(ctd: CTDTimeSeries, dl_path: Path):
    rdi = load_rdi(dl_path)
    z_m, izm = assign_bin_depths(rdi, ctd, looker="down")
    assert z_m.shape == (rdi.nens,)
    assert izm.shape == (rdi.nbin, rdi.nens)
    assert rdi.nbin == 25  # P16N cast003 DL has 25 bins

@pytest.mark.integration
def test_assign_depths_instrument_range(ctd: CTDTimeSeries, dl_path: Path):
    rdi = load_rdi(dl_path)
    z_m, _ = assign_bin_depths(rdi, ctd, looker="down")
    valid = z_m[np.isfinite(z_m)]
    assert np.max(valid) > 3000.0   # deep cast
    assert np.max(valid) < 5000.0
    assert np.min(valid) >= 0.0

@pytest.mark.integration
def test_assign_depths_bin0_deeper_than_instrument(ctd: CTDTimeSeries, dl_path: Path):
    rdi = load_rdi(dl_path)
    z_m, izm = assign_bin_depths(rdi, ctd, looker="down")
    valid = np.isfinite(z_m) & np.isfinite(izm[0, :])
    assert np.all(izm[0, valid] > z_m[valid])
