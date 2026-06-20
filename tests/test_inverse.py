import numpy as np
import pytest
from ladcp.solution.inverse import EnsembleData, SuperEnsemble, prepare_superensembles


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
