"""Unit tests for shear solution — src/ladcp/solution/shear.py."""
import numpy as np
import pytest
from ladcp.solution.shear import ShearProfile, compute_shear, _central_diff_shear, _bin_average_shear


def _constant_shear_inputs(
    nbin: int = 10,
    nens: int = 20,
    du_dz: float = 1e-3,
    dv_dz: float = -5e-4,
    dw_dz: float = 0.0,
    bin_spacing: float = 10.0,
    first_bin_depth: float = 10.0,
):
    """Synthetic field with exact linear shear du/dz everywhere."""
    izm = np.outer(
        np.arange(nbin) * bin_spacing + first_bin_depth,
        np.ones(nens),
    )  # (nbin, nens) positive depths
    u = izm * du_dz   # u = du_dz * z → shear = du_dz everywhere
    v = izm * dv_dz
    w = izm * dw_dz
    weight = np.ones((nbin, nens), dtype=np.float64)
    return u, v, w, izm, weight


def test_compute_shear_returns_shear_profile():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    assert isinstance(result, ShearProfile)


def test_shear_profile_fields_exist():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    for field in ("z", "u_shear", "v_shear", "w_shear",
                  "u_shear_err", "v_shear_err", "w_shear_err",
                  "n", "u_rel", "v_rel", "w_rel"):
        assert hasattr(result, field), f"ShearProfile missing field: {field}"


def test_shear_profile_arrays_same_length():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    nz = len(result.z)
    for arr in (result.u_shear, result.v_shear, result.w_shear,
                result.u_shear_err, result.v_shear_err, result.w_shear_err,
                result.n, result.u_rel, result.v_rel, result.w_rel):
        assert arr.shape == (nz,), f"Expected ({nz},), got {arr.shape}"


def test_z_axis_starts_at_half_dz():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    assert result.z[0] == pytest.approx(5.0)   # dz/2 = 5 m


def test_integrated_profile_is_zero_mean():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    assert np.mean(result.u_rel) == pytest.approx(0.0, abs=1e-12)
    assert np.mean(result.v_rel) == pytest.approx(0.0, abs=1e-12)


def test_central_diff_shear_exact_gradient():
    """Linear u = du_dz * z → shear = du_dz everywhere in interior."""
    du_dz = 1e-3
    nbin, nens = 6, 4
    bin_spacing = 8.0
    izm = np.outer(np.arange(nbin) * bin_spacing + 8.0, np.ones(nens))
    u = izm * du_dz
    v = np.zeros_like(u)
    w = np.zeros_like(u)
    weight_mask = np.ones((nbin, nens))
    su, sv, sw = _central_diff_shear(u, v, w, izm, weight_mask)
    # Interior bins (rows 1 to nbin-2) should equal du_dz exactly
    interior = su[1:-1, :]
    assert np.allclose(interior, du_dz, rtol=1e-10)


def test_central_diff_shear_boundary_nan():
    """First and last row must be NaN (no neighbour on one side)."""
    u, v, w, izm, weight = _constant_shear_inputs()
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    assert np.all(np.isnan(su[0, :]))
    assert np.all(np.isnan(su[-1, :]))


def test_central_diff_shear_weight_mask_applies():
    """Bins with weight ≤ weight_min become NaN shear (caller pre-masks)."""
    u, v, w, izm, _ = _constant_shear_inputs(nbin=6, nens=4)
    # weight_mask is 1 everywhere except column 2 → NaN
    mask = np.ones_like(u)
    mask[:, 2] = np.nan
    su, _, _ = _central_diff_shear(u, v, w, izm, mask)
    assert np.all(np.isnan(su[:, 2]))


def test_central_diff_shear_output_shape():
    u, v, w, izm, weight = _constant_shear_inputs(nbin=8, nens=5)
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    assert su.shape == (8, 5)
    assert sv.shape == (8, 5)
    assert sw.shape == (8, 5)


def test_bin_average_recovers_known_shear():
    """Constant shear field: bin mean should equal the true shear value."""
    du_dz = 2e-3
    u, v, w, izm, weight = _constant_shear_inputs(
        nbin=10, nens=30, du_dz=du_dz, dv_dz=0.0, dw_dz=0.0,
        bin_spacing=10.0, first_bin_depth=10.0,
    )
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, vsm, wsm, use, vse, wse, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    # Interior bins (have enough samples) should recover du_dz
    populated = nn > 2
    assert populated.sum() > 0, "No bins had more than 2 samples"
    assert np.allclose(usm[populated], du_dz, atol=1e-10)


def test_bin_average_outlier_rejection():
    """Inject a single large outlier; after 2σ editing it must not influence the mean."""
    nbin, nens = 6, 40
    du_dz = 1e-3
    u, v, w, izm, weight = _constant_shear_inputs(
        nbin=nbin, nens=nens, du_dz=du_dz, dv_dz=0.0, dw_dz=0.0,
        bin_spacing=10.0, first_bin_depth=10.0,
    )
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    # Inject spike in bin 3, first ensemble — 100× larger than normal
    su[3, 0] = 0.1
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, _, _, _, _, _, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    populated = nn > 2
    # After 2σ editing the mean in every bin should still be ≈ du_dz
    assert np.allclose(usm[populated], du_dz, atol=1e-4)


def test_bin_average_few_samples_returns_nan():
    """Bins with ≤ 2 valid samples return NaN (not enough to estimate std)."""
    nbin, nens = 4, 2   # only 2 ensembles → ≤ 2 finite shear estimates per bin
    u, v, w, izm, weight = _constant_shear_inputs(nbin=nbin, nens=nens)
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, _, _, _, _, _, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    assert np.all(np.isnan(usm[nn <= 2]))


def test_bin_average_output_shapes():
    u, v, w, izm, weight = _constant_shear_inputs(nbin=8, nens=20)
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, vsm, wsm, use, vse, wse, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    nz = len(z_bins)
    for arr in (usm, vsm, wsm, use, vse, wse, nn):
        assert arr.shape == (nz,)
