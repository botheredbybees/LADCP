"""Integration test: SADCP loader vs S4P cast 001 reference.

Requires TEST_DATA_DIR env var pointing to a directory containing 2018_S4P/.

Note: the reference 001.nc embeds u_sadcp/v_sadcp extracted from SADCP.mat
(a different product from os75nb_short.nc).  We compare mean absolute difference
rather than exact equality; tolerance 0.10 m/s.
"""
from __future__ import annotations

import os
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.sadcp import SADCPProfile, load_sadcp_nc


@pytest.fixture(scope="module")
def s4p_dir() -> Path:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        pytest.skip("TEST_DATA_DIR not set — see test_data/sources.md")
    p = Path(env) / "2018_S4P"
    if not p.exists():
        pytest.skip(f"2018_S4P directory not found at {p}")
    return p


@pytest.fixture(scope="module")
def sadcp_nc_path(s4p_dir: Path) -> Path:
    p = s4p_dir / "SADCP/os75nb/contour/os75nb_short.nc"
    if not p.exists():
        pytest.skip(f"SADCP NC not found: {p}")
    return p


@pytest.fixture(scope="module")
def ref_001(s4p_dir: Path) -> Path:
    p = s4p_dir / "001.nc"
    if not p.exists():
        pytest.skip(f"Reference 001.nc not found: {p}")
    return p


@pytest.fixture(scope="module")
def cast_001_window(ref_001: Path) -> tuple[float, float]:
    ds = netCDF4.Dataset(str(ref_001))
    tim = np.asarray(ds.variables["tim"][:], dtype=float)
    ds.close()
    return float(tim[0]), float(tim[-1])


@pytest.fixture(scope="module")
def sadcp_profile(sadcp_nc_path: Path, cast_001_window: tuple[float, float]) -> SADCPProfile:
    t_start, t_end = cast_001_window
    profile = load_sadcp_nc(
        sadcp_nc_path, t_start, t_end,
        lat=-70.45, lon=168.47,
    )
    if profile is None:
        pytest.skip("No SADCP records in cast 001 time window")
    return profile


@pytest.mark.integration
def test_sadcp_load_returns_profile(sadcp_profile: SADCPProfile) -> None:
    assert isinstance(sadcp_profile, SADCPProfile)


@pytest.mark.integration
def test_sadcp_profile_depth_range(sadcp_profile: SADCPProfile) -> None:
    """OS75 valid coverage during the cast window; expect at least 100–400 m."""
    assert float(sadcp_profile.z.min()) < 100.0
    assert float(sadcp_profile.z.max()) > 400.0


@pytest.mark.integration
def test_sadcp_profile_all_finite(sadcp_profile: SADCPProfile) -> None:
    assert np.all(np.isfinite(sadcp_profile.u))
    assert np.all(np.isfinite(sadcp_profile.v))
    assert np.all(np.isfinite(sadcp_profile.err))
    assert np.all(sadcp_profile.err > 0)


@pytest.mark.integration
def test_sadcp_vs_reference_mae(sadcp_profile: SADCPProfile, ref_001: Path) -> None:
    """Compare against u_sadcp/v_sadcp embedded in reference NC.

    These came from SADCP.mat (different processing path), so we use a
    generous 0.10 m/s MAE tolerance rather than requiring exact agreement.
    """
    ds = netCDF4.Dataset(str(ref_001))
    ref_vars = ds.variables
    if "z_sadcp" not in ref_vars or "u_sadcp" not in ref_vars:
        ds.close()
        pytest.skip("Reference 001.nc lacks z_sadcp/u_sadcp variables")
    ref_z = np.asarray(ref_vars["z_sadcp"][:], dtype=float)
    ref_u = np.asarray(ref_vars["u_sadcp"][:], dtype=float)
    ref_v = np.asarray(ref_vars["v_sadcp"][:], dtype=float)
    ds.close()

    our_u = np.interp(ref_z, sadcp_profile.z, sadcp_profile.u,
                      left=np.nan, right=np.nan)
    our_v = np.interp(ref_z, sadcp_profile.z, sadcp_profile.v,
                      left=np.nan, right=np.nan)

    valid = np.isfinite(ref_u) & np.isfinite(our_u)
    assert valid.sum() >= 10, f"Too few overlapping bins: {valid.sum()}"

    mae_u = float(np.mean(np.abs(our_u[valid] - ref_u[valid])))
    mae_v = float(np.mean(np.abs(our_v[valid] - ref_v[valid])))

    assert mae_u < 0.10, f"u MAE {mae_u:.4f} m/s exceeds 0.10 m/s (different SADCP products)"
    assert mae_v < 0.10, f"v MAE {mae_v:.4f} m/s exceeds 0.10 m/s"
