"""Integration test: inverse solver vs P16N cast 003 LDEO_IX reference.

Requires TEST_DATA_DIR env var pointing to a directory containing a
2015_P16N/ subdirectory with:
  003DL000.000   — Downlooker PD0 binary (beam coordinates, EX byte 0x04)
  003_01.cnv     — CTD time-series (binary SBE)
  003.nc         — LDEO_IX processed reference output

Data source: NCEI archive 0221195 (2015_P16N GO-SHIP cruise).
"""
from __future__ import annotations

import os
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.ctd import assign_bin_depths, load_ctd
from ladcp.ingestion.rdi import load_rdi
from ladcp.solution.inverse import EnsembleData, InverseResult, compute_inverse, prepare_superensembles
from ladcp.transforms.beam2earth import beam2earth

THETA_DEG = 20.0  # RDI Workhorse 300 kHz beam angle


@pytest.fixture(scope="module")
def test_data_dir() -> Path:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        pytest.skip("TEST_DATA_DIR not set — see test_data/sources.md")
    path = Path(env) / "2015_P16N"
    if not path.exists():
        pytest.skip(f"2015_P16N directory not found at {path}")
    return path


@pytest.fixture(scope="module")
def dl_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "003DL000.000"
    if not p.exists():
        pytest.skip(f"DL PD0 file not found: {p}")
    return p


@pytest.fixture(scope="module")
def cnv_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "003_01.cnv"
    if not p.exists():
        pytest.skip(f"CTD file not found: {p}")
    return p


@pytest.fixture(scope="module")
def ref_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "003.nc"
    if not p.exists():
        pytest.skip(f"Reference NetCDF not found: {p}")
    return p


@pytest.fixture(scope="module")
def inverse_result(dl_path: Path, cnv_path: Path) -> InverseResult:
    """Run full pipeline on P16N cast 003 raw data."""
    rdi = load_rdi(dl_path)
    ctd = load_ctd(cnv_path)

    # beam2earth: file is in beam coordinates (EX byte 0x04), need explicit transform.
    # The rdi.u/v/w/e fields hold raw beam data for beam-coord files.
    u_earth, v_earth, w_earth = beam2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG,
        gimbaled=True,
    )

    # assign_bin_depths returns positive-down z_m (nens,) and izm (nbin, nens).
    z_m, izm_pos = assign_bin_depths(rdi, ctd, looker="down")

    # EnsembleData depth convention: negative = below surface.
    z_neg = -z_m           # (nens,) negative down
    izm_neg = -izm_pos     # (nbin, nens) negative down

    # Correlation-based weight: mean over 4 beams, normalised to 0–1.
    weight = np.nanmean(rdi.corr.astype(np.float64), axis=2) / 128.0

    # Bottom-track velocity in Earth frame (approximate: using raw beam-frame data
    # directly here since no beam→earth transform is applied to BT; bvels held fixed
    # at 2 cm/s nominal).  btrack_vel_ms is (4, nens), already in m/s.
    bt_u = rdi.btrack_vel_ms[0]  # beam 1 ≈ East proxy
    bt_v = rdi.btrack_vel_ms[1]  # beam 2 ≈ North proxy
    bt_w = rdi.btrack_vel_ms[2]  # beam 3 ≈ vertical proxy
    bvel = np.stack([bt_u, bt_v, bt_w], axis=1)    # (nens, 3)
    bvels = np.full_like(bvel, 0.02)                # 2 cm/s nominal std
    hbot = np.nanmean(rdi.btrack_range_m, axis=0)  # (nens,) mean of 4-beam ranges

    ens = EnsembleData(
        u=u_earth,
        v=v_earth,
        w=w_earth,
        weight=weight,
        izm=izm_neg,
        z=z_neg,
        time_jul=rdi.time_julian,  # already Julian days
        bvel=bvel,
        bvels=bvels,
        hbot=hbot,
        izd=np.arange(rdi.nbin),
        izu=np.array([], dtype=int),
        slat=np.full(rdi.nens, np.nan),
        slon=np.full(rdi.nens, np.nan),
    )

    se = prepare_superensembles(ens, dz=16.0)
    return compute_inverse(se)


@pytest.mark.integration
def test_inverse_result_is_inverse_result(inverse_result: InverseResult):
    """compute_inverse must return an InverseResult."""
    assert isinstance(inverse_result, InverseResult)


@pytest.mark.integration
def test_inverse_profile_has_z_and_velocity(inverse_result: InverseResult):
    """Profile must have depth and velocity arrays of equal length."""
    assert inverse_result.z.shape == inverse_result.u.shape == inverse_result.v.shape
    assert len(inverse_result.z) > 0


@pytest.mark.integration
def test_inverse_depth_range(inverse_result: InverseResult, ref_path: Path):
    """Computed profile must reach at least 80 % of the reference max depth."""
    ds = netCDF4.Dataset(ref_path)
    ref_z = np.array(ds.variables["z"][:])  # z coordinate in 003.nc (positive m)
    ref_u = np.array(ds.variables["u"][:])
    ds.close()
    ref_valid = np.isfinite(ref_u)
    ref_max_depth = float(ref_z[ref_valid].max())

    result_max_depth = float(inverse_result.z.max())
    assert result_max_depth > 0.8 * ref_max_depth, (
        f"Profile only reaches {result_max_depth:.0f} m vs ref {ref_max_depth:.0f} m"
    )


@pytest.mark.integration
@pytest.mark.xfail(strict=False, reason="Pipeline gaps: missing BT beam→earth rotation, no QC, no GPS/SADCP constraint; remove once pipeline complete")
def test_inverse_u_rmse(inverse_result: InverseResult, ref_path: Path):
    """RMS error in u vs LDEO_IX reference must be < 0.05 m/s (bins with nvel >= 3)."""
    ds = netCDF4.Dataset(ref_path)
    ref_z = np.array(ds.variables["z"][:])      # positive m, increasing downward
    ref_u = np.array(ds.variables["u"][:])      # m/s
    ref_nvel = np.array(ds.variables["nvel"][:])
    ds.close()

    result_u = np.interp(ref_z, inverse_result.z, inverse_result.u,
                         left=np.nan, right=np.nan)

    # Restrict comparison to well-observed reference bins
    valid = np.isfinite(ref_u) & np.isfinite(result_u) & (ref_nvel >= 3)
    assert valid.sum() > 10, (
        f"Too few overlapping depth bins to compare: {valid.sum()} bins"
    )

    rmse = float(np.sqrt(np.mean((result_u[valid] - ref_u[valid]) ** 2)))
    assert rmse < 0.05, f"u RMSE {rmse:.4f} m/s exceeds 0.05 m/s tolerance"


@pytest.mark.integration
@pytest.mark.xfail(strict=False, reason="Pipeline gaps: missing BT beam→earth rotation, no QC, no GPS/SADCP constraint; remove once pipeline complete")
def test_inverse_v_rmse(inverse_result: InverseResult, ref_path: Path):
    """RMS error in v vs LDEO_IX reference must be < 0.05 m/s (bins with nvel >= 3)."""
    ds = netCDF4.Dataset(ref_path)
    ref_z = np.array(ds.variables["z"][:])
    ref_v = np.array(ds.variables["v"][:])
    ref_nvel = np.array(ds.variables["nvel"][:])
    ds.close()

    result_v = np.interp(ref_z, inverse_result.z, inverse_result.v,
                         left=np.nan, right=np.nan)

    valid = np.isfinite(ref_v) & np.isfinite(result_v) & (ref_nvel >= 3)
    assert valid.sum() > 10, (
        f"Too few overlapping depth bins to compare: {valid.sum()} bins"
    )

    rmse = float(np.sqrt(np.mean((result_v[valid] - ref_v[valid]) ** 2)))
    assert rmse < 0.05, f"v RMSE {rmse:.4f} m/s exceeds 0.05 m/s tolerance"
