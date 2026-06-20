"""Integration tests: shear solution against P16N cast 003 LDEO reference.

Requires TEST_DATA_DIR env var pointing to a directory containing:
  2015_P16N/003DL000.000   — Downlooker PD0 binary
  2015_P16N/003_01.cnv     — CTD time-series (binary SBE)
  2015_P16N/003.nc         — LDEO_IX processed reference output
"""
import os
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.ctd import assign_bin_depths, load_ctd
from ladcp.ingestion.rdi import load_rdi
from ladcp.solution.shear import ShearProfile, compute_shear


@pytest.fixture
def test_data_dir() -> Path:
    path = Path(os.environ.get("TEST_DATA_DIR", "test_data"))
    if not path.exists():
        pytest.skip("TEST_DATA_DIR not populated — see test_data/sources.md")
    return path


@pytest.fixture
def dl_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003DL000.000"
    if not p.exists():
        pytest.skip(f"DL PD0 file not found: {p}")
    return p


@pytest.fixture
def cnv_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003_01.cnv"
    if not p.exists():
        pytest.skip(f"CTD file not found: {p}")
    return p


@pytest.fixture
def ref_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003.nc"
    if not p.exists():
        pytest.skip(f"Reference NetCDF not found: {p}")
    return p


@pytest.fixture
def shear_result(dl_path: Path, cnv_path: Path) -> ShearProfile:
    rdi = load_rdi(dl_path)
    ctd = load_ctd(cnv_path)
    _, izm = assign_bin_depths(rdi, ctd, looker="down")
    weight = np.ones((rdi.nbin, rdi.nens), dtype=np.float64)
    return compute_shear(rdi.u, rdi.v, rdi.w, izm, weight, dz=10.0)


@pytest.mark.integration
def test_shear_profile_is_shear_profile(shear_result: ShearProfile):
    assert isinstance(shear_result, ShearProfile)


@pytest.mark.integration
def test_shear_depth_axis_starts_at_5m(shear_result: ShearProfile):
    assert shear_result.z[0] == pytest.approx(5.0, abs=0.1)


@pytest.mark.integration
def test_shear_depth_axis_covers_cast_depth(shear_result: ShearProfile, ref_path: Path):
    ds = netCDF4.Dataset(ref_path)
    ref_z_max = float(np.max(ds.variables["z"][:]))
    ds.close()
    # Our profile must reach at least 80 % of the reference maximum depth
    assert shear_result.z[-1] >= 0.8 * ref_z_max, (
        f"z_max={shear_result.z[-1]:.0f} m vs ref {ref_z_max:.0f} m"
    )


@pytest.mark.integration
def test_integrated_profile_zero_mean(shear_result: ShearProfile):
    assert np.mean(shear_result.u_rel) == pytest.approx(0.0, abs=1e-12)
    assert np.mean(shear_result.v_rel) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.integration
def test_shear_magnitude_plausible(shear_result: ShearProfile, ref_path: Path):
    """Integrated velocity should be the same order of magnitude as the reference."""
    ds = netCDF4.Dataset(ref_path)
    ref_u = np.array(ds.variables["u_shear_method"][:])
    ds.close()
    ref_rms = float(np.sqrt(np.nanmean(ref_u**2)))
    our_rms = float(np.sqrt(np.nanmean(shear_result.u_rel**2)))
    # Without full weight editing our RMS may differ by up to 5×, but same order
    assert our_rms > ref_rms / 10.0, f"Our u_rel RMS {our_rms:.4f} is too small vs ref {ref_rms:.4f}"
    assert our_rms < ref_rms * 10.0, f"Our u_rel RMS {our_rms:.4f} is too large vs ref {ref_rms:.4f}"


@pytest.mark.integration
def test_shear_profile_n_populated(shear_result: ShearProfile):
    """More than half the depth bins should have valid shear estimates."""
    assert (shear_result.n > 2).sum() > len(shear_result.n) // 2
