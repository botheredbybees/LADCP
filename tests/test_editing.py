"""Unit tests for ladcp.qa.editing.edit_sidelobes."""
from __future__ import annotations

import numpy as np

from ladcp.qa.editing import edit_sidelobes
from ladcp.solution.inverse import EnsembleData


def _make_ens(
    izm: np.ndarray,
    z: np.ndarray,
    weight: np.ndarray | None = None,
    hbot: np.ndarray | None = None,
) -> EnsembleData:
    """Minimal EnsembleData for sidelobe editing tests."""
    n_bins, n_ens = izm.shape
    if weight is None:
        weight = np.ones((n_bins, n_ens))
    if hbot is None:
        hbot = np.full(n_ens, np.nan)
    return EnsembleData(
        u=np.zeros((n_bins, n_ens)),
        v=np.zeros((n_bins, n_ens)),
        w=np.zeros((n_bins, n_ens)),
        weight=weight,
        izm=izm,
        z=z,
        time_jul=np.zeros(n_ens),
        bvel=np.zeros((n_ens, 3)),
        bvels=np.full((n_ens, 3), 0.02),
        hbot=hbot,
        izd=np.arange(n_bins),
        izu=np.array([], dtype=int),
        slat=np.full(n_ens, np.nan),
        slon=np.full(n_ens, np.nan),
    )


def test_surface_sidelobe_masks_shallow_bins():
    # z=-100, theta=20°, cell=8m:
    #   f = 1 - cos(20°) ≈ 0.06031
    #   margin = 1.5 * 8 = 12.0
    #   zlim_surface = 0.06031 * (-100) - 12 ≈ -18.03
    # bin at izm=-5  → -5 > -18.03  → contaminated (shallower than limit)
    # bin at izm=-30 → -30 < -18.03 → clean
    izm = np.array([[-5.0], [-30.0]])   # (2, 1)
    z   = np.array([-100.0])            # (1,)
    ens = _make_ens(izm, z)
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "shallow bin should be masked"
    assert result.weight[1, 0] == 1.0, "deep bin should be untouched"


def test_bottom_sidelobe_masks_deep_bins():
    # z=-100, zbottom=150 (explicit), theta=20°, cell=8m:
    #   hab = -100 + 150 = 50
    #   zlim_bot = -150 + 0.06031*50 + 12 ≈ -134.98
    #   zlim_surface ≈ -18.03  (neither test bin is shallower)
    # bin at izm=-120 → -120 > -134.98 → clean
    # bin at izm=-140 → -140 < -134.98 → bottom-contaminated
    izm = np.array([[-120.0], [-140.0]])
    z   = np.array([-100.0])
    ens = _make_ens(izm, z)
    result = edit_sidelobes(ens, zbottom=150.0, theta_deg=20.0, cell_size_m=8.0)
    assert result.weight[0, 0] == 1.0, "mid-column bin should be untouched"
    assert np.isnan(result.weight[1, 0]), "near-bottom bin should be masked"


def test_no_hbot_skips_bottom_edit():
    # hbot all NaN → auto-derived zbottom = nanmedian(NaN) = NaN → skip bottom mask
    # bin at -5 is surface-contaminated (z=-100, zlim_surface≈-18.03)
    # bin at -140 would be bottom-contaminated with zbottom=150, but is NOT masked here
    izm = np.array([[-5.0], [-140.0]])
    z   = np.array([-100.0])
    ens = _make_ens(izm, z, hbot=np.array([np.nan]))
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "surface bin still masked"
    assert result.weight[1, 0] == 1.0, "bottom mask skipped when no BT data"


def test_existing_nan_weights_preserved():
    # A bin well inside the safe zone has a pre-existing NaN weight.
    # After editing it should still be NaN (copy preserves it).
    # z=-100, zlim_surface≈-18.03: bin at -50 is not surface-contaminated.
    izm    = np.array([[-50.0]])
    z      = np.array([-100.0])
    weight = np.array([[np.nan]])
    ens    = _make_ens(izm, z, weight=weight)
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "pre-existing NaN must not be cleared"


def test_explicit_zbottom_overrides_auto():
    # hbot=[200] → auto zbottom = 50+200 = 250 → zlim_bot ≈ -225.94
    # explicit zbottom=100  → hab=50,          zlim_bot ≈  -84.98
    # bin at izm=-90:
    #   surface: zlim_surface = 0.06031*(-50)-12 ≈ -15.02; -90 < -15.02 → NOT surface-masked
    #   with explicit zbottom=100: -90 < -84.98 → contaminated
    #   with auto  zbottom=250: -90 > -225.94  → clean
    # Passing explicit zbottom=100 should mask the bin.
    izm = np.array([[-90.0]])
    z   = np.array([-50.0])
    ens = _make_ens(izm, z, hbot=np.array([200.0]))
    result = edit_sidelobes(ens, zbottom=100.0, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "explicit zbottom must be used, not auto-derived"


def test_returns_new_ensemble_not_mutated():
    # The function must return a new EnsembleData and not alter the input weight array.
    izm = np.array([[-5.0], [-50.0]])   # shallow bin will be surface-masked
    z   = np.array([-100.0])
    ens = _make_ens(izm, z)
    original_weight = ens.weight.copy()
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    np.testing.assert_array_equal(ens.weight, original_weight)  # not mutated
    assert result is not ens                                     # new object
    assert np.isnan(result.weight[0, 0])                        # sanity: mask applied


def test_surface_mask_is_per_ensemble():
    # Two ensembles at very different depths produce different zlim_surface values,
    # so the same bin depth maps to different mask outcomes — proving the mask
    # is computed per-ensemble via broadcasting, not from a scalar.
    # z[0]=-200: zlim_surface = 0.06031*(-200) - 12 ≈ -24.06
    # z[1]=-20:  zlim_surface = 0.06031*(-20)  - 12 ≈ -13.21
    # bin at izm=-18:
    #   ens 0: -18 > -24.06 → masked (shallower than limit for deep ADCP)
    #   ens 1: -18 < -13.21 → clean  (deeper than limit for shallow ADCP)
    z   = np.array([-200.0, -20.0])
    izm = np.array([[-18.0, -18.0]])   # (1, 2) same bin depth, two ensembles
    ens = _make_ens(izm, z)
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "deep-ADCP ensemble: shallow bin is masked"
    assert result.weight[0, 1] == 1.0, "shallow-ADCP ensemble: bin is below limit"
