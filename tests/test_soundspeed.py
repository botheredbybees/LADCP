from __future__ import annotations

import numpy as np

from ladcp.solution.inverse import EnsembleData
from ladcp.transforms.soundspeed import (
    apply_sound_speed_correction,
    depth_to_pressure,
    sound_speed,
)

# --- sound_speed (sounds.m / Chen & Millero 1977, UNESCO 44) ---


def test_sound_speed_matches_sounds_m():
    # Parity target is sounds.m as executed, not its comment: run under
    # Octave 9.2, sounds(10000, 40, 40) = 1732.139394 (measured 2026-07-10).
    # The in-file UNESCO check-value comment (1731.995) does NOT hold for
    # the code as written -- a ~0.14 m/s (1e-4 relative) discrepancy
    # inherited from the FORTRAN translation. Since the correction only
    # uses the RATIO ss/sv, this has no practical effect, and matching
    # LDEO's actual arithmetic is what validation requires.
    c = sound_speed(np.array([10000.0]), np.array([40.0]), 40.0)
    assert abs(c[0] - 1732.139394) < 1e-4


def test_sound_speed_surface_freshwater_range():
    # ~1449 m/s at S=35, T=10, P=0 (canonical seawater value ballpark)
    c = sound_speed(np.array([0.0]), np.array([10.0]), 35.0)
    assert 1440.0 < c[0] < 1500.0


# --- depth_to_pressure (press.m / GEOSECS) ---


def test_depth_to_pressure_zero_at_surface():
    assert abs(depth_to_pressure(np.array([0.0]))[0]) < 0.05


def test_depth_to_pressure_monotonic_and_plausible():
    p = depth_to_pressure(np.array([1000.0, 4000.0]))
    assert 1000.0 < p[0] < 1020.0     # ~1% over depth in dbar
    assert 4020.0 < p[1] < 4120.0
    assert p[1] > p[0]


# --- apply_sound_speed_correction (getdpthi.m lines 182-207 + 428-441) ---


def _make_ens(n_ul: int = 2, n_dl: int = 2, n_ens: int = 4) -> EnsembleData:
    n_bins = n_ul + n_dl
    z = np.full(n_ens, -1000.0)
    # izm = z + signed bin offsets: UL rows above instrument, DL rows below
    # rows: UL far, UL near, DL near, DL far
    offsets = np.array([16.0, 8.0, -8.0, -16.0])
    izm = z[np.newaxis, :] + offsets[:, np.newaxis]
    return EnsembleData(
        u=np.full((n_bins, n_ens), 0.5),
        v=np.full((n_bins, n_ens), -0.5),
        w=np.full((n_bins, n_ens), 1.0),
        weight=np.ones((n_bins, n_ens)),
        izm=izm,
        z=z,
        time_jul=np.arange(n_ens, dtype=float),
        bvel=np.full((n_ens, 3), 0.2),
        bvels=np.full((n_ens, 3), 0.02),
        hbot=np.full(n_ens, 150.0),
        izd=np.array([2, 3]),
        izu=np.array([1, 0]),
        slat=np.full(n_ens, np.nan),
        slon=np.full(n_ens, np.nan),
    )


def test_velocities_scaled_per_instrument():
    ens = _make_ens()
    n_ens = 4
    ss = np.full(n_ens, 1515.0)
    sv_dl = np.full(n_ens, 1500.0)   # sc_dl = 1.01
    sv_ul = np.full(n_ens, 1530.0)   # sc_ul = 1515/1530
    out = apply_sound_speed_correction(ens, ss=ss, sv_dl=sv_dl, sv_ul=sv_ul)
    np.testing.assert_allclose(out.u[ens.izd, :], 0.5 * 1.01)
    np.testing.assert_allclose(out.u[ens.izu, :], 0.5 * 1515.0 / 1530.0)
    np.testing.assert_allclose(out.w[ens.izd, :], 1.0 * 1.01)
    # weight untouched
    np.testing.assert_array_equal(out.weight, ens.weight)


def test_bottom_track_scaled_with_dl_ratio():
    ens = _make_ens()
    ss = np.full(4, 1515.0)
    out = apply_sound_speed_correction(
        ens, ss=ss, sv_dl=np.full(4, 1500.0), sv_ul=np.full(4, 1500.0)
    )
    np.testing.assert_allclose(out.bvel, 0.2 * 1.01)
    np.testing.assert_allclose(out.hbot, 150.0 * 1.01)


def test_izm_bin_offsets_scaled_but_not_z():
    ens = _make_ens()
    ss = np.full(4, 1515.0)
    out = apply_sound_speed_correction(
        ens, ss=ss, sv_dl=np.full(4, 1500.0), sv_ul=np.full(4, 1530.0)
    )
    # instrument depth itself must NOT be scaled (it comes from CTD pressure)
    np.testing.assert_array_equal(out.z, ens.z)
    # DL far bin: offset -16 -> -16*1.01; izm = -1000 - 16.16
    np.testing.assert_allclose(out.izm[3, :], -1000.0 - 16.0 * 1.01)
    # UL far bin: offset +16 -> * (1515/1530)
    np.testing.assert_allclose(out.izm[0, :], -1000.0 + 16.0 * 1515.0 / 1530.0)


def test_input_not_mutated():
    ens = _make_ens()
    u_orig = ens.u.copy()
    izm_orig = ens.izm.copy()
    out = apply_sound_speed_correction(
        ens, ss=np.full(4, 1515.0), sv_dl=np.full(4, 1500.0), sv_ul=np.full(4, 1500.0)
    )
    np.testing.assert_array_equal(ens.u, u_orig)
    np.testing.assert_array_equal(ens.izm, izm_orig)
    assert out is not ens
