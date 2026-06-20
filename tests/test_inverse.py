import numpy as np
import pytest
from ladcp.solution.inverse import (
    EnsembleData,
    SuperEnsemble,
    prepare_superensembles,
    _flatten_obs,
    _build_obs_matrix,
    _build_ctd_matrix,
    _apply_weights,
    _add_smoothness,
    _add_zero_mean,
)
import scipy.sparse


def _make_ens(*, n_bins: int = 5, n_ens: int = 40,
               u_val: float = 0.5, v_val: float = 0.1) -> EnsembleData:
    """Synthetic EnsembleData where CTD descends 5 m per ensemble."""
    z = np.linspace(-20.0, -20.0 - 5.0 * (n_ens - 1), n_ens)
    izm = np.zeros((n_bins, n_ens))
    for b in range(n_bins):
        izm[b] = z - (b + 0.5) * 8.0
    return EnsembleData(
        u=np.full((n_bins, n_ens), u_val),
        v=np.full((n_bins, n_ens), v_val),
        w=np.zeros((n_bins, n_ens)),
        weight=np.ones((n_bins, n_ens)),
        izm=izm,
        z=z,
        time_jul=np.linspace(2457000.0, 2457001.0, n_ens),
        bvel=np.full((n_ens, 3), np.nan),
        bvels=np.full((n_ens, 3), np.nan),
        hbot=np.full(n_ens, np.nan),
        izd=np.arange(n_bins),
        izu=np.array([], dtype=int),
        slat=np.full(n_ens, np.nan),
        slon=np.full(n_ens, np.nan),
    )


def test_prepare_superensembles_reduces_n_ens():
    """Super-ensemble averaging must collapse multiple raw ensembles."""
    ens = _make_ens(n_ens=40)
    se = prepare_superensembles(ens, dz=20.0)
    assert se.ru.shape[1] < 40


def test_prepare_superensembles_preserves_mean_velocity():
    """Constant velocity field should survive depth-window averaging unchanged."""
    ens = _make_ens(u_val=0.5, v_val=0.1)
    se = prepare_superensembles(ens, dz=20.0)
    assert np.allclose(np.nanmean(se.ru), 0.5, atol=0.02)
    assert np.allclose(np.nanmean(se.rv), 0.1, atol=0.02)


def test_prepare_superensembles_default_dz():
    """Default dz (inferred from izm[:, 0] bin spacing) must give the same shape as explicit dz."""
    ens = _make_ens(n_ens=40)
    se_default = prepare_superensembles(ens)
    # default dz = median(|diff(izm[:, 0])|) = 8.0 m (bin spacing in _make_ens)
    se_explicit = prepare_superensembles(ens, dz=8.0)
    assert se_default.ru.shape == se_explicit.ru.shape


def test_prepare_superensembles_bvel_preserved():
    """Bottom-track velocity should be averaged into super-ensembles."""
    ens = _make_ens(n_ens=40)
    ens.bvel[:] = [0.1, -0.2, -1.0]
    se = prepare_superensembles(ens, dz=20.0)
    assert np.allclose(np.nanmean(se.bvel[0]), 0.1, atol=0.02)


def test_build_obs_matrix_shape():
    """obs matrix rows = n_obs, cols = number of unique depth bins."""
    izv = np.array([10.0, 20.0, 30.0, 10.0])  # 3 unique bins
    A = _build_obs_matrix(izv, dz=10.0)
    assert A.shape == (4, 3)


def test_build_obs_matrix_one_nonzero_per_row():
    """Each observation maps to exactly one depth bin."""
    izv = np.array([5.0, 15.0, 25.0, 35.0])
    A = _build_obs_matrix(izv, dz=10.0)
    assert np.allclose(np.asarray(A.sum(axis=1)).ravel(), 1.0)


def test_build_ctd_matrix_shape():
    """ctd matrix rows = n_obs, cols = n_se."""
    jprof = np.array([0, 0, 1, 1, 2, 2], dtype=int)
    A = _build_ctd_matrix(jprof, n_se=3)
    assert A.shape == (6, 3)
    assert np.allclose(np.asarray(A.sum(axis=1)).ravel(), 1.0)


def test_flatten_obs_removes_nan_weight():
    """Observations with NaN or zero weight must be excluded."""
    ens = _make_ens(n_bins=3, n_ens=20)
    se = prepare_superensembles(ens, dz=10.0)
    # Force some NaN weights
    se.weight[:, :2] = np.nan
    d_u, d_v, izv, jprof, wm = _flatten_obs(se, velerr=0.05, weightmin=0.05)
    assert np.all(np.isfinite(d_u))
    assert np.all(np.isfinite(wm))
    assert np.all(wm >= 0.05)


def test_apply_weights_scales_rows():
    """_apply_weights must scale each row of A and d by the observation weight."""
    # Simple 4-observation system: 3 depth bins, 2 super-ensembles
    izv = np.array([10.0, 20.0, 30.0, 20.0])
    jprof = np.array([0, 0, 1, 1], dtype=int)
    A_ocean = _build_obs_matrix(izv, dz=10.0)
    A_ctd = _build_ctd_matrix(jprof, n_se=2)
    d = np.array([0.1, 0.2, 0.3, 0.4])
    wm = np.array([1.0, 2.0, 0.5, 3.0])

    A_o, A_c, d_w, idx_down, idx_up = _apply_weights(A_ocean, A_ctd, d, wm)

    # Row k of A_o == row k of A_ocean * wm[k]
    for k in range(4):
        assert np.allclose(A_o[k], np.asarray(A_ocean[k].todense()).ravel() * wm[k])
    assert np.allclose(d_w, d * wm)

    # idx_down and idx_up together cover all row indices
    all_covered = np.union1d(idx_down, idx_up)
    assert set(all_covered) == set(range(len(d)))


def test_add_smoothness_increases_rows():
    """Smoothness adds curvature rows for interior columns."""
    n_obs, n_zbins, n_se = 20, 8, 5
    A_o = np.random.rand(n_obs, n_zbins)
    A_c = np.random.rand(n_obs, n_se)
    d = np.random.rand(n_obs)
    A_o2, A_c2, d2 = _add_smoothness(A_o, A_c, d, smoofac=1.0)
    # Must add at least n_zbins - 2 curvature rows (interior bins)
    assert A_o2.shape[0] > n_obs
    assert A_c2.shape[0] == A_o2.shape[0]
    assert len(d2) == A_o2.shape[0]


def test_add_smoothness_zero_smoofac_still_runs():
    """smoofac=0 must not raise and must add at least boundary rows."""
    A_o = np.eye(6)
    A_c = np.zeros((6, 3))
    d = np.zeros(6)
    A_o2, A_c2, d2 = _add_smoothness(A_o, A_c, d, smoofac=0.0)
    assert A_o2.shape[0] >= 6


def test_add_zero_mean_appends_one_row():
    """Zero-mean adds exactly one constraint row."""
    A_o = np.eye(5)
    A_c = np.zeros((5, 3))
    d = np.ones(5)
    A_o2, A_c2, d2 = _add_zero_mean(A_o, A_c, d)
    assert A_o2.shape[0] == 6
    assert d2[-1] == 0.0  # RHS = 0 for zero-mean
