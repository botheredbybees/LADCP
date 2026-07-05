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
    _add_bottom_track,
    _add_sadcp,
    _add_barotropic,
    InverseParams,
    InverseResult,
    compute_inverse,
    rotup2down,
    offsetup2down,
    _medianan_na,
)
import scipy.sparse


def _make_ens(*, n_bins: int = 5, n_ens: int = 40,
               u_val: float = 0.5, v_val: float = 0.1,
               noise: float = 0.0) -> EnsembleData:
    """Synthetic EnsembleData where CTD descends 5 m per ensemble.

    Parameters
    ----------
    noise:
        Std of Gaussian noise added to u/v. A non-zero value prevents
        degenerate zero-variance windows (ruvs=0 → weight=NaN).
    """
    rng = np.random.default_rng(42)
    z = np.linspace(-20.0, -20.0 - 5.0 * (n_ens - 1), n_ens)
    izm = np.zeros((n_bins, n_ens))
    for b in range(n_bins):
        izm[b] = z - (b + 0.5) * 8.0
    u_arr = np.full((n_bins, n_ens), u_val)
    v_arr = np.full((n_bins, n_ens), v_val)
    if noise > 0.0:
        u_arr = u_arr + rng.normal(0.0, noise, u_arr.shape)
        v_arr = v_arr + rng.normal(0.0, noise, v_arr.shape)
    return EnsembleData(
        u=u_arr,
        v=v_arr,
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
    """smoofac=0 disables all smoothness constraints — no rows added."""
    A_o = np.eye(6)
    A_c = np.zeros((6, 3))
    d = np.zeros(6)
    A_o2, A_c2, d2 = _add_smoothness(A_o, A_c, d, smoofac=0.0)
    assert A_o2.shape[0] == 6  # no rows added when smoofac=0


def test_add_smoothness_two_column_matrix():
    """Boundary stencil must not crash with minimum viable 2-column matrix."""
    A_o = np.eye(2)
    A_c = np.zeros((2, 2))
    d = np.zeros(2)
    A_o2, A_c2, d2 = _add_smoothness(A_o, A_c, d, smoofac=1.0)
    assert A_o2.shape[0] >= 2


def test_add_zero_mean_appends_one_row():
    """Zero-mean adds exactly one constraint row."""
    A_o = np.eye(5)
    A_c = np.zeros((5, 3))
    d = np.ones(5)
    A_o2, A_c2, d2 = _add_zero_mean(A_o, A_c, d)
    assert A_o2.shape[0] == 6
    assert d2[-1] == 0.0  # RHS = 0 for zero-mean


def test_add_bottom_track_appends_rows():
    """One constraint row per ensemble with valid bottom track."""
    n_obs, n_zbins, n_se = 10, 5, 4
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.zeros((n_obs, n_se))
    d = np.zeros(n_obs)
    bvel = np.array([0.1, np.nan, -0.2, 0.0])   # 3 valid, 1 NaN
    bvels = np.array([0.01, 0.01, 0.01, 0.01])
    A_o2, A_c2, d2 = _add_bottom_track(A_o, A_c, d, bvel, bvels,
                                        botfac=1.0, velerr=0.05)
    # Should append one row per finite bvel
    n_finite = int(np.sum(np.isfinite(bvel)))
    assert A_c2.shape[0] == n_obs + n_finite


def test_add_barotropic_appends_one_row():
    """Barotropic constraint adds exactly one row to each component."""
    n_obs, n_zbins, n_se = 10, 5, 4
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.ones((n_obs, n_se))  # Use ones so col_scale is non-zero
    d_u = np.zeros(n_obs)
    d_v = np.zeros(n_obs)
    dt = np.ones(n_se) * 100.0
    A_ou, A_cu, du, A_ov, A_cv, dv = _add_barotropic(
        A_o, A_c, d_u, A_o.copy(), A_c.copy(), d_v,
        u_ship=0.5, v_ship=0.1, dt=dt, barofac=1.0,
    )
    assert A_cu.shape[0] == n_obs + 1
    assert du[-1] != 0.0   # RHS = -u_ship * weight


def test_add_barotropic_weight_formula():
    """Row weight matches MATLAB lainbaro: barvelerr=2*nav_error/T, fac=sqrt(sum|Ac|).

    For cast-003-like inputs: T=11157s, nav_error=30m, velerr=0.05
    barvelerr = 2*30/11157 ≈ 0.005378 m/s  (matches 003.nc log)
    fac_nav   = velerr/barvelerr ≈ 9.30
    fac       = sqrt(sum(|A_ctd|)) = sqrt(40) for ones(10,4)
    RHS       = -u_ship * barofac * fac_nav * fac
    """
    import math
    n_obs, n_zbins, n_se = 10, 5, 4
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.ones((n_obs, n_se))
    d_u = np.zeros(n_obs)
    d_v = np.zeros(n_obs)
    # Use T = 11157 s (11 SEs * ~1014 s each — cast-003 scale)
    T = 11157.0
    dt = np.full(n_se, T / n_se)

    nav_error = 30.0
    velerr = 0.05
    u_ship = 1.0  # 1 m/s for easy verification

    A_ou, A_cu, du, A_ov, A_cv, dv = _add_barotropic(
        A_o, A_c, d_u, A_o.copy(), A_c.copy(), d_v,
        u_ship=u_ship, v_ship=0.0, dt=dt, barofac=1.0,
        nav_error=nav_error, velerr=velerr,
    )

    barvelerr = 2.0 * nav_error / T          # ≈ 0.005378
    fac_nav = velerr / barvelerr             # ≈ 9.298
    fac = math.sqrt(float(np.abs(A_c).sum()))  # sqrt(40) ≈ 6.325
    expected_rhs = -u_ship * fac_nav * fac

    assert abs(barvelerr - 0.005378) < 1e-5, f"barvelerr={barvelerr}"
    assert abs(du[-1] - expected_rhs) < 1e-10, f"du[-1]={du[-1]}, expected={expected_rhs}"


def test_solve_lsq_identity():
    """Identity system must recover exact solution."""
    from ladcp.solution.inverse import _solve_lsq
    A = np.eye(4)
    d = np.array([1.0, 2.0, 3.0, 4.0])
    m, me = _solve_lsq(A, d)
    assert np.allclose(m, d)
    assert np.all(me >= 0)


def test_solve_lsq_overdetermined_consistent():
    """Consistent overdetermined system must find exact solution."""
    from ladcp.solution.inverse import _solve_lsq
    # 3 equations, 2 unknowns: x=1, y=2
    A = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    d = np.array([1.0, 2.0, 3.0])
    m, me = _solve_lsq(A, d)
    assert np.allclose(m, [1.0, 2.0], atol=1e-10)
    assert np.all(np.isfinite(me))


def test_solve_lsq_error_shape():
    """Error vector must have same length as solution vector."""
    from ladcp.solution.inverse import _solve_lsq
    A = np.random.rand(20, 5)
    d = np.random.rand(20)
    m, me = _solve_lsq(A, d)
    assert m.shape == me.shape == (5,)


def test_compute_inverse_returns_result():
    """compute_inverse must return an InverseResult with z, u, v arrays.

    Note: noise=0.02 is required — a perfectly constant velocity field yields
    ruvs=0 in every depth window, which propagates NaN into wm and leaves zero
    valid observations. The noise gives the window statistics a non-zero spread.
    """
    ens = _make_ens(n_bins=5, n_ens=60, noise=0.02)
    se = prepare_superensembles(ens, dz=10.0)
    result = compute_inverse(se)
    assert isinstance(result, InverseResult)
    assert result.z.shape == result.u.shape == result.v.shape
    assert len(result.z) > 0


def test_compute_inverse_zero_mean_no_constraint():
    """Without external constraints, result mean should be near zero.

    Note: noise=0.02 keeps ruvs non-NaN so observations survive the weight filter.
    The zero-mean fallback (_add_zero_mean) is triggered when botfac=barofac=0.
    """
    ens = _make_ens(n_bins=5, n_ens=60, u_val=0.0, v_val=0.0, noise=0.02)
    se = prepare_superensembles(ens, dz=10.0)
    params = InverseParams(barofac=0.0, botfac=0.0, sadcpfac=0.0)
    result = compute_inverse(se, params=params)
    assert abs(result.ubar) < 0.05
    assert abs(result.vbar) < 0.05


def test_compute_inverse_down_up_same_length():
    """Down-cast and up-cast profiles must have same length as full profile.

    Note: noise=0.02 is required for the same reason as test_compute_inverse_returns_result.
    """
    ens = _make_ens(n_bins=5, n_ens=60, noise=0.02)
    se = prepare_superensembles(ens, dz=10.0)
    result = compute_inverse(se)
    assert result.u_do.shape == result.u.shape
    assert result.u_up.shape == result.u.shape


def test_compute_inverse_with_bottom_track():
    """Bottom-track branch in compute_inverse must not crash and produce a result."""
    ens = _make_ens(n_bins=5, n_ens=60, noise=0.02)
    # Set valid bottom-track on all ensembles
    ens.bvel[:] = [0.05, -0.02, -0.5]
    ens.bvels[:] = [0.02, 0.02, 0.02]
    se = prepare_superensembles(ens, dz=10.0)
    params = InverseParams(botfac=1.0, barofac=0.0, sadcpfac=0.0)
    result = compute_inverse(se, params=params)
    assert isinstance(result, InverseResult)
    assert result.u.shape == result.z.shape


def test_compute_inverse_with_barotropic():
    """Barotropic branch in compute_inverse must not crash and produce a result."""
    ens = _make_ens(n_bins=5, n_ens=60, noise=0.02)
    se = prepare_superensembles(ens, dz=10.0)
    params = InverseParams(botfac=0.0, barofac=1.0, sadcpfac=0.0)
    result = compute_inverse(se, params=params, u_ship=0.1, v_ship=-0.05)
    assert isinstance(result, InverseResult)
    assert result.u.shape == result.z.shape


def test_barotropic_activates_for_stationary_ship():
    """u_ship=0.0 must still activate the barotropic constraint (GPS says ship stopped).

    The constraint sets mean(u_ctd) = u_ship, not u_ocean directly, so we can't
    assert ubar ≈ 0 without bottom-track data.  We verify it runs and produces a
    result that differs from the no-GPS zero-mean fallback.
    """
    ens = _make_ens(n_bins=5, n_ens=60, u_val=0.5, noise=0.02)
    se = prepare_superensembles(ens, dz=10.0)
    params = InverseParams(botfac=0.0, barofac=1.0, sadcpfac=0.0)
    result_gps = compute_inverse(se, params=params, u_ship=0.0, v_ship=0.0)
    result_no_gps = compute_inverse(se, params=params, u_ship=None, v_ship=None)
    assert isinstance(result_gps, InverseResult)
    # The two solutions must differ (GPS constraint changes the system)
    assert not np.allclose(result_gps.u, result_no_gps.u)


def test_barotropic_disabled_when_no_gps():
    """u_ship=None must disable the barotropic constraint entirely."""
    ens = _make_ens(n_bins=5, n_ens=60, u_val=1.0, v_val=0.5, noise=0.02)
    se = prepare_superensembles(ens, dz=10.0)
    params = InverseParams(botfac=0.0, barofac=1.0, sadcpfac=0.0)
    # No GPS → falls back to zero-mean; u_ship=None should not crash.
    result_no_gps = compute_inverse(se, params=params, u_ship=None, v_ship=None)
    assert isinstance(result_no_gps, InverseResult)
    # Same call with default (also None) should behave identically.
    result_default = compute_inverse(se, params=params)
    np.testing.assert_allclose(result_no_gps.u, result_default.u)


def test_compute_inverse_with_sadcp():
    """SADCP branch in compute_inverse must not crash and produce a result."""
    ens = _make_ens(n_bins=5, n_ens=60, noise=0.02)
    se = prepare_superensembles(ens, dz=10.0)
    params = InverseParams(botfac=0.0, barofac=0.0, sadcpfac=1.0)
    sadcp_z = np.array([50.0, 100.0, 150.0])
    sadcp_u = np.array([0.1, 0.05, -0.02])
    sadcp_v = np.array([-0.05, 0.02, 0.01])
    sadcp_err = np.array([0.03, 0.03, 0.03])
    result = compute_inverse(
        se, params=params,
        sadcp_z=sadcp_z, sadcp_u=sadcp_u, sadcp_v=sadcp_v, sadcp_err=sadcp_err,
    )
    assert isinstance(result, InverseResult)
    assert result.u.shape == result.z.shape


def test_add_sadcp_appends_rows():
    """Each finite SADCP measurement adds one row to A_ocean."""
    n_obs, n_zbins, n_se = 10, 8, 4
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.zeros((n_obs, n_se))
    d = np.zeros(n_obs)
    sadcp_z = np.array([15.0, 25.0, 35.0, np.nan])  # 3 valid
    sadcp_vel = np.array([0.1, 0.2, 0.3, np.nan])
    sadcp_err = np.array([0.02, 0.02, 0.02, 0.02])
    A_o2, A_c2, d2 = _add_sadcp(
        A_o, A_c, d,
        sadcp_z=sadcp_z, sadcp_vel=sadcp_vel, sadcp_err=sadcp_err,
        dz=10.0, sadcpfac=1.0, velerr=0.05,
    )
    assert A_o2.shape[0] == n_obs + 3  # 3 finite measurements


def test_add_sadcp_zeros_in_A_ctd():
    """SADCP constraint rows must have no A_ctd contribution."""
    n_obs, n_zbins, n_se = 5, 4, 3
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.zeros((n_obs, n_se))
    d = np.zeros(n_obs)
    A_o2, A_c2, d2 = _add_sadcp(
        A_o, A_c, d,
        sadcp_z=np.array([10.0]), sadcp_vel=np.array([0.5]),
        sadcp_err=np.array([0.02]), dz=10.0, sadcpfac=1.0, velerr=0.05,
    )
    assert A_c2[-1].sum() == 0.0  # last row of A_ctd = zeros


def _make_dual_ens(n_ul: int = 3, n_dl: int = 3, n_ens: int = 40) -> EnsembleData:
    """Synthetic dual-instrument EnsembleData with unit eastward flow."""
    n_bins = n_ul + n_dl
    z = np.linspace(-20.0, -220.0, n_ens)
    izm = np.tile(z, (n_bins, 1))
    return EnsembleData(
        u=np.ones((n_bins, n_ens)),
        v=np.zeros((n_bins, n_ens)),
        w=np.zeros((n_bins, n_ens)),
        weight=np.ones((n_bins, n_ens)),
        izm=izm,
        z=z,
        time_jul=np.linspace(2457000.0, 2457001.0, n_ens),
        bvel=np.full((n_ens, 3), np.nan),
        bvels=np.full((n_ens, 3), np.nan),
        hbot=np.full(n_ens, np.nan),
        izd=np.arange(n_ul, n_bins),
        izu=np.arange(n_ul),
        slat=np.full(n_ens, np.nan),
        slon=np.full(n_ens, np.nan),
    )


def test_rotup2down_constant_offset_is_noop():
    """A constant DL-UL compass offset is absorbed into hoff: hrot ~ 0."""
    ens = _make_dual_ens()
    n_ens = ens.z.size
    hdg_dl = np.full(n_ens, 100.0)
    hdg_ul = np.full(n_ens, 10.0)  # constant 90-deg mounting offset
    out, hrot = rotup2down(ens, hdg_dl, hdg_ul)
    np.testing.assert_allclose(hrot, 0.0, atol=1e-10)
    np.testing.assert_allclose(out.u, ens.u, atol=1e-12)
    np.testing.assert_allclose(out.v, ens.v, atol=1e-12)


def test_rotup2down_splits_fluctuation():
    """UL heading fluctuation +/-4 deg: DL and UL counter-rotate by 2 deg.

    If the UL compass reads 4 deg high, its Earth vector is 4 deg clockwise
    of the DL's; the correction rotates UL 2 deg CCW and DL 2 deg CW.
    """
    ens = _make_dual_ens()
    n_ens = ens.z.size
    d = np.where(np.arange(n_ens) % 2 == 0, 4.0, -4.0)  # zero-mean fluctuation
    hdg_dl = np.full(n_ens, 100.0)
    hdg_ul = 10.0 + d
    out, hrot = rotup2down(ens, hdg_dl, hdg_ul)
    np.testing.assert_allclose(hrot, -d, atol=1e-8)
    c2, s2 = np.cos(np.radians(2.0)), np.sin(np.radians(2.0))
    k = 0  # ensemble with d=+4 (hrot=-4): DL rotated 2 deg CW, UL 2 deg CCW
    assert abs(out.u[ens.izd[0], k] - c2) < 1e-8
    assert abs(out.v[ens.izd[0], k] - (-s2)) < 1e-8
    assert abs(out.u[ens.izu[0], k] - c2) < 1e-8
    assert abs(out.v[ens.izu[0], k] - (+s2)) < 1e-8


def test_rotup2down_returns_new_ensemble():
    """Input EnsembleData must not be mutated."""
    ens = _make_dual_ens()
    n_ens = ens.z.size
    u_before = ens.u.copy()
    hdg_dl = np.full(n_ens, 100.0)
    hdg_ul = 10.0 + np.where(np.arange(n_ens) % 2 == 0, 4.0, -4.0)
    out, _ = rotup2down(ens, hdg_dl, hdg_ul)
    np.testing.assert_array_equal(ens.u, u_before)
    assert out is not ens


def _make_biased_dual_ens(
    bias: float, n_ul: int = 3, n_dl: int = 3, n_ens: int = 10,
) -> EnsembleData:
    """Dual-instrument EnsembleData: UL reads 2*bias high in u, DL matches truth."""
    ens = _make_dual_ens(n_ul=n_ul, n_dl=n_dl, n_ens=n_ens)
    ens.izm[:, :] = -100.0  # constant positive depth 100 m for every bin/ensemble
    ens.u[ens.izu, :] = 1.0 + 2.0 * bias
    ens.u[ens.izd, :] = 1.0
    ens.bvel[:, 0] = 0.2  # arbitrary finite bottom-track u, to check it also shifts
    ens.bvel[:, 1] = 0.0
    ens.bvel[:, 2] = 0.0
    return ens


def test_medianan_na_matches_matlab_round_half_away_from_zero():
    """medianan.m's `round([-na:na]+n/2)` uses half-away-from-zero rounding.

    n=21, na=2 gives x/2=10.5 -> offsets [8.5,9.5,10.5,11.5,12.5] -> MATLAB
    round() (half away from zero) picks 1-based [9,10,11,12,13]. numpy's
    default round-half-to-even would instead pick [8,10,10,12,12] --
    duplicating two indices and dropping one. Distinct, non-uniform values
    are required to expose this (a uniform column can't distinguish which
    indices were averaged).
    """
    n = 21
    col = np.arange(n, dtype=float)[:, np.newaxis]  # column 0..20, already sorted
    y = _medianan_na(col, na=2)
    # MATLAB 1-based [9,10,11,12,13] -> 0-based [8,9,10,11,12] -> mean = 10.0
    assert abs(y[0] - 10.0) < 1e-10


def test_offsetup2down_splits_residual_toward_consensus():
    """UL reads 2*bias high vs a first guess matching DL (u=1): both converge to 1+bias.

    prepinv.m's offsetup2down shifts UL/DL by half the UL-DL residual, in
    opposite directions, so neither instrument is treated as "more correct".
    """
    bias = 0.6
    ens = _make_biased_dual_ens(bias)
    dr_z = np.array([0.0, 200.0])
    dr_u = np.array([1.0, 1.0])
    dr_v = np.array([0.0, 0.0])
    out, uoff = offsetup2down(ens, dr_z, dr_u, dr_v)
    np.testing.assert_allclose(out.u[ens.izu, :], 1.0 + bias, atol=1e-10)
    np.testing.assert_allclose(out.u[ens.izd, :], 1.0 + bias, atol=1e-10)
    np.testing.assert_allclose(out.v, 0.0, atol=1e-10)
    # Bottom track shifts with the DL correction (-uoff/2), matching rotup2down.
    np.testing.assert_allclose(out.bvel[:, 0], 0.2 + bias, atol=1e-10)


def test_offsetup2down_zero_factor_is_noop():
    """factor=0 disables the correction entirely (matches p.offsetup2down=0)."""
    ens = _make_biased_dual_ens(bias=0.6)
    dr_z = np.array([0.0, 200.0])
    dr_u = np.array([1.0, 1.0])
    dr_v = np.array([0.0, 0.0])
    out, uoff = offsetup2down(ens, dr_z, dr_u, dr_v, factor=0.0)
    np.testing.assert_allclose(out.u, ens.u, atol=1e-12)
    np.testing.assert_allclose(out.v, ens.v, atol=1e-12)
    np.testing.assert_allclose(out.bvel, ens.bvel, atol=1e-12, equal_nan=True)


def test_offsetup2down_returns_new_ensemble():
    """Input EnsembleData must not be mutated."""
    ens = _make_biased_dual_ens(bias=0.6)
    u_before = ens.u.copy()
    dr_z = np.array([0.0, 200.0])
    dr_u = np.array([1.0, 1.0])
    dr_v = np.array([0.0, 0.0])
    out, _ = offsetup2down(ens, dr_z, dr_u, dr_v)
    np.testing.assert_array_equal(ens.u, u_before)
    assert out is not ens
