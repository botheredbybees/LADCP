"""Integration tests: load_rdi() against 2015 P16N cast 003 raw PD0 files.

Data source: NCEI archive 0221195 (2015_P16N GO-SHIP cruise)
Files: test_data/2015_P16N/003DL000.000, 003UL000.000
Reference output: test_data/2015_P16N/003.nc (LDEO_IX processed)
"""

import numpy as np
import pytest

from ladcp.ingestion.rdi import load_rdi


@pytest.fixture
def dl_path(test_data_dir):
    p = test_data_dir / "2015_P16N" / "003DL000.000"
    if not p.exists():
        pytest.skip(f"P16N downlooker PD0 not found at {p}")
    return p


@pytest.fixture
def ul_path(test_data_dir):
    p = test_data_dir / "2015_P16N" / "003UL000.000"
    if not p.exists():
        pytest.skip(f"P16N uplooker PD0 not found at {p}")
    return p


@pytest.mark.integration
def test_dl_loads(dl_path):
    d = load_rdi(dl_path)
    assert d.nens > 100
    assert d.nbin > 0


@pytest.mark.integration
def test_dl_ensemble_count(dl_path):
    """~8970 ensembles for cast 003 (~3.6-hour cast at 1-second pings)."""
    d = load_rdi(dl_path)
    assert 6000 < d.nens < 12000, f"Unexpected ensemble count: {d.nens}"


@pytest.mark.integration
def test_dl_bin_geometry(dl_path):
    """300 kHz Workhorse: 8 m bins, first bin centre ~8 m."""
    d = load_rdi(dl_path)
    assert abs(d.blen_m - 8.0) < 0.5, f"blen_m={d.blen_m}"
    assert d.nbin == 25, f"nbin={d.nbin}"
    assert 6.0 < d.dist_m < 12.0, f"dist_m={d.dist_m}"


@pytest.mark.integration
def test_dl_heading_in_range(dl_path):
    """Headings are in [0, 360]."""
    d = load_rdi(dl_path)
    valid = d.heading[np.isfinite(d.heading)]
    assert len(valid) > 0
    assert np.all(valid >= 0) and np.all(valid <= 360)


@pytest.mark.integration
def test_dl_velocity_finite_fraction(dl_path):
    """At least 50% of DL velocities are finite."""
    d = load_rdi(dl_path)
    frac = np.isfinite(d.u).mean()
    assert frac > 0.5, f"Too many NaN velocities: {frac:.0%} finite"


@pytest.mark.integration
def test_dl_time_monotone(dl_path):
    d = load_rdi(dl_path)
    assert np.all(np.diff(d.time_julian) >= 0), "Time is not monotonically increasing"


@pytest.mark.integration
def test_dl_bottom_track_finite(dl_path):
    """Some BT data present — cast reaches near-bottom."""
    d = load_rdi(dl_path)
    frac = np.isfinite(d.btrack_range_m).mean()
    assert frac > 0.1, f"Expected some BT data, got {frac:.0%} finite"


@pytest.mark.integration
def test_ul_loads(ul_path):
    d = load_rdi(ul_path)
    assert d.nens > 100
    assert d.nbin > 0


@pytest.mark.integration
def test_dl_ul_time_overlap(dl_path, ul_path):
    """DL and UL time spans should overlap (same cast, same wire)."""
    dl = load_rdi(dl_path)
    ul = load_rdi(ul_path)
    dl_start, dl_end = dl.time_julian[0], dl.time_julian[-1]
    ul_start, ul_end = ul.time_julian[0], ul.time_julian[-1]
    overlap_start = max(dl_start, ul_start)
    overlap_end = min(dl_end, ul_end)
    assert overlap_end > overlap_start, "DL and UL time ranges do not overlap"


@pytest.mark.integration
def test_dl_coord_transform_beam_gimbaled(dl_path):
    """P16N 003DL000.000 recorded in beam-frame with gimbaled tilt (EX=0x04)."""
    d = load_rdi(dl_path)
    assert d.coord_transform == 4, (
        f"Expected EX=4 (beam+gimbaled), got {d.coord_transform}"
    )
    assert (d.coord_transform >> 3) & 0x03 == 0, (
        "Expected beam-frame (coord bits = 0b00)"
    )
    assert bool(d.coord_transform & 0x04), "Expected gimbaled bit set"
