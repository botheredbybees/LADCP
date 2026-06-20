"""Integration tests: beam2earth() applied to P16N cast 003 DL data.

Verifies Earth-frame velocities are oceanographically plausible.
This is a sanity-check level test — not numerical validation against 002.nc.

Data source: NCEI archive 0221195 (2015_P16N GO-SHIP cruise)
File: test_data/2015_P16N/003DL000.000
EX byte: 0x04 (Beam coordinates, gimbaled=True, binmap=False)
"""

import numpy as np
import pytest

from ladcp.ingestion.rdi import load_rdi
from ladcp.transforms.beam2earth import beam2earth

THETA_DEG = 20.0  # RDI Workhorse 300 kHz


@pytest.fixture
def dl_data(test_data_dir):
    p = test_data_dir / "2015_P16N" / "003DL000.000"
    if not p.exists():
        pytest.skip(f"P16N downlooker PD0 not found at {p}")
    return load_rdi(p)


@pytest.mark.integration
def test_transform_runs_without_error(dl_data):
    d = dl_data
    u, v, w = beam2earth(
        d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG, gimbaled=True
    )
    assert u.shape == d.u.shape
    assert v.shape == d.v.shape
    assert w.shape == d.w.shape


@pytest.mark.integration
def test_earth_frame_finite_fraction_comparable_to_beam(dl_data):
    """Transform preserves NaN pattern: finite fraction ≈ same as input."""
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    beam_finite = np.isfinite(d.u).mean()
    earth_finite = np.isfinite(u).mean()
    # Vz uses all four beams, so NaN spread increases in Earth frame (typically ~12%).
    # Allow up to 15% change.
    assert abs(earth_finite - beam_finite) < 0.15, (
        f"beam finite={beam_finite:.2f}, earth finite={earth_finite:.2f}"
    )


@pytest.mark.integration
def test_horizontal_velocity_plausible(dl_data):
    """Mean horizontal speed should be oceanographically plausible (0–2 m/s)."""
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    mean_spd = np.sqrt(np.nanmean(u**2) + np.nanmean(v**2))
    assert mean_spd < 2.0, f"Mean horizontal speed too large: {mean_spd:.3f} m/s"
    # Should see *some* velocity (not degenerate/all-zero)
    assert mean_spd > 0.001, f"Mean speed suspiciously low: {mean_spd:.6f} m/s"


@pytest.mark.integration
def test_vertical_velocity_plausible_magnitude(dl_data):
    """Raw beam2earth w includes instrument descent rate (~1 m/s); values stay bounded.

    Note: RMS w > RMS horizontal is expected for raw LADCP data — the CTD rosette
    descends/ascends at ~1 m/s, dominating the w signal. The signed mean is small
    because descent and ascent phases cancel across the cast.
    """
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    rms_w = np.sqrt(np.nanmean(w**2))
    assert rms_w < 3.0, f"RMS vertical speed too large: {rms_w:.3f} m/s"
    # Signed mean is small — descent and ascent roughly cancel over a full cast
    assert abs(np.nanmean(w)) < 0.5, (
        f"Mean w suspiciously large: {np.nanmean(w):.3f} m/s"
    )


@pytest.mark.integration
def test_velocity_magnitudes_not_nan_dominated(dl_data):
    """At least 50% of Earth-frame cells have finite values."""
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    assert np.isfinite(u).mean() > 0.5
