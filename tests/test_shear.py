"""Unit tests for shear solution — src/ladcp/solution/shear.py."""
import numpy as np
import pytest
from ladcp.solution.shear import ShearProfile, compute_shear


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
