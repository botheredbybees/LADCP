"""Integration tests: load_rdi() against I7N cast 002 raw PD0 files."""

import numpy as np
import pytest

from ladcp.ingestion.rdi import load_rdi


@pytest.fixture
def dl_path(test_data_dir):
    p = test_data_dir / "raw" / "002DL000.000"
    if not p.exists():
        pytest.skip(f"Downlooker PD0 not found at {p}. Run Task 1 first.")
    return p


@pytest.fixture
def ul_path(test_data_dir):
    p = test_data_dir / "raw" / "002UL000.000"
    if not p.exists():
        pytest.skip(f"Uplooker PD0 not found at {p}. Run Task 1 first.")
    return p


@pytest.mark.integration
def test_dl_loads(dl_path):
    """Downlooker file loads without error."""
    d = load_rdi(dl_path)
    assert d.nens > 100, f"Expected >100 ensembles, got {d.nens}"
    assert d.nbin > 0


@pytest.mark.integration
def test_dl_ensemble_count(dl_path):
    """Downlooker has expected number of ensembles (~600 for cast 002)."""
    d = load_rdi(dl_path)
    assert 400 < d.nens < 1200, f"Unexpected ensemble count: {d.nens}"


@pytest.mark.integration
def test_dl_bin_geometry(dl_path):
    """Bin length ~8 m, dist_m ~6 m for 300 kHz Workhorse."""
    d = load_rdi(dl_path)
    assert abs(d.blen_m - 8.0) < 1.0, f"blen_m={d.blen_m}"
    assert 4.0 < d.dist_m < 20.0, f"dist_m={d.dist_m}"


@pytest.mark.integration
def test_dl_heading_in_range(dl_path):
    """All headings are in [0, 360)."""
    d = load_rdi(dl_path)
    valid = d.heading[np.isfinite(d.heading)]
    assert len(valid) > 0
    assert np.all(valid >= 0) and np.all(valid < 360)


@pytest.mark.integration
def test_dl_velocity_finite_fraction(dl_path):
    """At least 50% of velocity data is finite (not NaN)."""
    d = load_rdi(dl_path)
    frac = np.isfinite(d.u).mean()
    assert frac > 0.5, f"Too many NaN velocities: {frac:.1%} finite"


@pytest.mark.integration
def test_dl_time_monotone(dl_path):
    """Ensemble times are monotonically increasing."""
    d = load_rdi(dl_path)
    dt = np.diff(d.time_julian)
    assert np.all(dt >= 0), "Time is not monotonically increasing"


@pytest.mark.integration
def test_ul_loads(ul_path):
    """Uplooker file loads without error."""
    d = load_rdi(ul_path)
    assert d.nens > 100
    assert d.nbin > 0


@pytest.mark.integration
def test_dl_bottom_track_finite(dl_path):
    """Bottom track ranges are mostly finite for deep cast with BT data."""
    d = load_rdi(dl_path)
    frac = np.isfinite(d.btrack_range_m).mean()
    assert frac > 0.1, f"Expected some BT data, got {frac:.1%} finite"
