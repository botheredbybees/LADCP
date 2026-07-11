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


# --- edit_outliers (loadrdi.m::outlier port) ---

from ladcp.qa.editing import edit_mask_bins, edit_outliers


def _make_ens_outlier(n_ul: int = 3, n_dl: int = 3, n_ens: int = 300, seed: int = 7):
    """Combined UL+DL ensemble with smooth fields for outlier tests.

    1 s pings -> loadrdi's 5-minute block = 300 ensembles = one block here.
    """
    rng = np.random.default_rng(seed)
    n_bins = n_ul + n_dl
    t = 2457000.0 + np.arange(n_ens) / 86400.0
    base = 0.3 + 0.05 * np.sin(np.arange(n_ens) / 40.0)
    u = base + rng.normal(0, 0.02, (n_bins, n_ens))
    v = -base + rng.normal(0, 0.02, (n_bins, n_ens))
    w = 1.0 + rng.normal(0, 0.02, (n_bins, n_ens))
    ens = EnsembleData(
        u=u, v=v, w=w,
        weight=np.ones((n_bins, n_ens)),
        izm=np.tile(-np.arange(n_bins, dtype=float)[:, None] * 8 - 50, (1, n_ens)),
        z=np.full(n_ens, -500.0),
        time_jul=t,
        bvel=rng.normal(0, 0.02, (n_ens, 3)),
        bvels=np.full((n_ens, 3), 0.02),
        hbot=np.full(n_ens, 100.0) + rng.normal(0, 1.0, n_ens),
        izd=np.arange(n_ul, n_ul + n_dl),
        izu=np.arange(n_ul - 1, -1, -1),
        slat=np.full(n_ens, np.nan),
        slon=np.full(n_ens, np.nan),
    )
    return ens


def test_edit_outliers_masks_dl_velocity_spike():
    ens = _make_ens_outlier()
    dl_row = ens.izd[1]
    ens.u[dl_row, 50] = 15.0  # physically impossible spike
    result = edit_outliers(ens)
    assert np.isnan(result.u[dl_row, 50])
    assert np.isnan(result.weight[dl_row, 50])
    # unlike the weight-only editors, outlier() NaNs all velocity components
    assert np.isnan(result.v[dl_row, 50])
    assert np.isnan(result.w[dl_row, 50])


def test_edit_outliers_masks_ul_block_independently():
    ens = _make_ens_outlier()
    ul_row = ens.izu[0]
    ens.w[ul_row, 120] = -9.0
    result = edit_outliers(ens)
    assert np.isnan(result.w[ul_row, 120])
    # clean cells survive both sweeps
    frac_masked = np.mean(~np.isfinite(result.u))
    assert frac_masked < 0.05


def test_edit_outliers_bottom_track_spike():
    ens = _make_ens_outlier()
    ens.bvel[70, 0] = 5.0
    result = edit_outliers(ens)
    assert np.all(np.isnan(result.bvel[70]))
    assert np.isnan(result.hbot[70])
    assert np.isfinite(result.bvel[71]).all()


def test_edit_outliers_does_not_mutate_input():
    ens = _make_ens_outlier()
    ens.u[ens.izd[0], 10] = 20.0
    u_orig = ens.u.copy()
    result = edit_outliers(ens)
    np.testing.assert_array_equal(ens.u, u_orig)
    assert result is not ens


# --- edit_mask_bins (edit_data.m bin masking, incl. zero-blanking bin 1) ---


def test_edit_mask_bins_masks_weight_rows_only():
    ens = _make_ens_outlier()
    # DL bin 0 and UL bin 0 (0-based instrument bins), as edit_data.m does
    # when the blanking distance is zero.
    result = edit_mask_bins(ens, dn_bins=[0], up_bins=[0])
    assert np.all(np.isnan(result.weight[ens.izd[0], :]))
    assert np.all(np.isnan(result.weight[ens.izu[0], :]))
    # velocities untouched (edit_data.m only NaNs d.weight)
    assert np.isfinite(result.u).all()
    # other rows untouched
    assert np.isfinite(result.weight[ens.izd[1], :]).all()


def test_edit_mask_bins_ignores_out_of_range():
    ens = _make_ens_outlier()
    result = edit_mask_bins(ens, dn_bins=[99], up_bins=[])
    assert np.isfinite(result.weight).all()


# --- build_ldeo_weights (loadrdi.m lines 408-533) ---

from ladcp.qa.editing import build_ldeo_weights


def _weight_inputs(nbin=6, nens=200, seed=13):
    rng = np.random.default_rng(seed)
    cm = rng.uniform(80, 128, (nbin, nens))
    ts = rng.uniform(60, 90, (nbin, nens))
    pitch = rng.normal(0, 0.5, nens)
    roll = rng.normal(0, 0.5, nens)
    v = rng.normal(0, 0.3, (nbin, nens))
    w = 1.0 + rng.normal(0, 0.1, (nbin, nens))
    izd = np.arange(3, 6)
    izu = np.array([2, 1, 0])
    return cm, ts, pitch, roll, v, w, izd, izu


def test_weights_normalized_by_median_of_max():
    cm, ts, pitch, roll, v, w, izd, izu = _weight_inputs()
    ts[:] = 0.0  # kill the echo penalty for this test
    wt = build_ldeo_weights(cm, ts, pitch, roll, v, w, izd, izu)
    col_max = np.nanmax(wt, axis=0)
    assert abs(np.nanmedian(col_max) - 1.0) < 1e-9


def test_weights_tilt_masks_whole_ensembles():
    cm, ts, pitch, roll, v, w, izd, izu = _weight_inputs()
    pitch[7] = 25.0    # tilt > 22 deg
    wt = build_ldeo_weights(cm, ts, pitch, roll, v, w, izd, izu)
    assert np.all(np.isnan(wt[:, 7]))
    # neighbors 6 and 8 are masked too (the 25-deg JUMP trips the
    # tilt-derivative check -- faithful to loadrdi.m); 9 is clean
    assert np.isfinite(wt[:, 9]).all()


def test_weights_tilt_derivative_masks_ensembles():
    cm, ts, pitch, roll, v, w, izd, izu = _weight_inputs()
    roll[50] = 8.0     # jump: neighbors get tiltd > 4
    wt = build_ldeo_weights(cm, ts, pitch, roll, v, w, izd, izu)
    assert np.all(np.isnan(wt[:, 50]))


def test_weights_echo_penalty_reduces_strong_echo():
    cm, ts, pitch, roll, v, w, izd, izu = _weight_inputs()
    cm[:] = 100.0
    ts[:] = 70.0
    ts[2, 30] = 90.0   # strong echo anomaly in row 2
    wt = build_ldeo_weights(cm, ts, pitch, roll, v, w, izd, izu)
    # the anomalous cell is the row max -> factor (1 - 1^1.5) = 0
    assert wt[2, 30] < 1e-12
    assert wt[2, 31] > 0.5


def test_weights_non_pinging_ensembles_masked_per_block():
    cm, ts, pitch, roll, v, w, izd, izu = _weight_inputs()
    ts[:] = 0.0
    # UL block (rows 0-2) flat w AND v at ensembles 10-12 -> non-pinging
    w[izu, 10:13] = 1.0
    v[izu, 10:13] = 0.1
    wt = build_ldeo_weights(cm, ts, pitch, roll, v, w, izd, izu)
    assert np.all(np.isnan(wt[izu][:, 10:13]))
    assert np.isfinite(wt[izd][:, 10:13]).all()  # DL block untouched


# --- edit_error_velocity (loadrdi.m elim, lines 1117-1134) ---

from ladcp.qa.editing import edit_error_velocity
from ladcp.transforms.beam2earth import janus_error_velocity


def test_error_velocity_edit_masks_large_e():
    u = np.ones((3, 4))
    v = np.ones((3, 4)) * 2
    w = np.ones((3, 4)) * 3
    e = np.zeros((3, 4))
    e[1, 2] = 0.7
    e[2, 0] = -0.9
    u2, v2, w2 = edit_error_velocity(u, v, w, e)
    for a in (u2, v2, w2):
        assert np.isnan(a[1, 2]) and np.isnan(a[2, 0])
        assert np.isfinite(a).sum() == 10
    # inputs untouched
    assert np.isfinite(u).all()


def test_error_velocity_edit_nan_e_kept():
    # missing-beam / 3-beam cells (e = NaN or 0) never trip the threshold
    u = np.ones((2, 2))
    e = np.array([[np.nan, 0.0], [0.49, 0.51]])
    u2, _, _ = edit_error_velocity(u, u, u, e)
    assert np.isfinite(u2[0, 0]) and np.isfinite(u2[0, 1])
    assert np.isfinite(u2[1, 0]) and np.isnan(u2[1, 1])


def test_janus_error_velocity_formula():
    theta = 20.0
    b1, b2, b3, b4 = 1.0, 2.0, 0.5, 0.25
    expected = (b1 + b2 - b3 - b4) / (4.0 * np.cos(np.radians(theta)))
    got = janus_error_velocity(
        np.array([b1]), np.array([b2]), np.array([b3]), np.array([b4]), theta)
    assert abs(got[0] - expected) < 1e-15
    # homogeneous flow across beams -> zero error velocity
    same = np.array([0.8])
    assert abs(janus_error_velocity(same, same, same, same, theta)[0]) < 1e-15


# --- edit_ppi (edit_data.m 208-262, previous-ping interference) ---

from ladcp.qa.editing import edit_ppi


def _make_ppi_ens(dt_s: float = 1.6, n_ens: int = 5, zbottom: float = 5000.0):
    """DL-only ensemble set: 3 bins per ensemble, one inside the PPI band.

    PPI band centre = -zbottom + 1500*dt/2*cos(20 deg) (~ -zbottom + 1128 m
    for dt=1.6 s); band half-width 90 m.
    """
    hab_centre = 1500.0 * dt_s / 2.0 * np.cos(np.radians(20.0))
    z_in = -zbottom + hab_centre           # inside the band
    z_above = -zbottom + hab_centre + 200  # above the band
    z_below = -zbottom + 50                # near bottom, below the band
    izm = np.tile(np.array([[z_above], [z_in], [z_below]]), (1, n_ens))
    z = np.full(n_ens, -3000.0)
    ens = _make_ens(izm, z, hbot=np.full(n_ens, zbottom - 3000.0))
    ens = ens.__class__(**{**ens.__dict__,
                           'time_jul': np.arange(n_ens) * dt_s / 86400.0})
    return ens


def test_ppi_masks_band_dl_only():
    # band centre ~1128 m above bottom needs max_hab > that (A16N used 1200)
    ens = _make_ppi_ens()
    out = edit_ppi(ens, npng=1, beam_angle_deg=20.0, max_hab_m=1200.0)
    # ensemble 0 has no dt -> untouched (edit_data.m starts at ensemble 2)
    assert out.weight[1, 0] == 1.0
    # band bin masked for ensembles 2..end
    assert np.all(np.isnan(out.weight[1, 1:]))
    # bins outside the band untouched
    assert np.all(out.weight[0, :] == 1.0)
    assert np.all(out.weight[2, :] == 1.0)
    # input not mutated
    assert np.all(ens.weight == 1.0)


def test_ppi_max_hab_disables_edit():
    # band centre ~1128 m above bottom > max_hab 1000 -> no edit at all
    ens = _make_ppi_ens(dt_s=1.6)
    out = edit_ppi(ens, npng=1, beam_angle_deg=20.0, max_hab_m=1000.0)
    assert np.all(np.isfinite(out.weight))


def test_ppi_multi_ping_dt_scaling():
    # npng=2 halves the effective ping interval -> band centre ~564 m
    # above bottom; the bin at ~1128 m is no longer masked
    ens = _make_ppi_ens(dt_s=1.6)
    out = edit_ppi(ens, npng=2, beam_angle_deg=20.0)
    assert np.all(np.isfinite(out.weight[1, :]))
