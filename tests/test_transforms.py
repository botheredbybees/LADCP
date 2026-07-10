"""Unit tests for src/ladcp/transforms/beam2earth.py."""

import numpy as np

from ladcp.transforms.beam2earth import (
    beam2earth,
    beam2xyz,
    reconstruct_3beam,
    uvrot,
)

THETA = 20.0  # RDI Workhorse 300 kHz


class TestBeam2xyz:
    def test_zero_beams_give_zero(self):
        b = np.zeros((5, 10))
        Vx, Vy, Vz = beam2xyz(b, b, b, b, THETA)
        assert np.all(Vx == 0.0)
        assert np.all(Vy == 0.0)
        assert np.all(Vz == 0.0)

    def test_output_shape_matches_input(self):
        rng = np.random.default_rng(42)
        b = rng.standard_normal((25, 100))
        Vx, Vy, Vz = beam2xyz(b, b, b, b, THETA)
        assert Vx.shape == (25, 100)
        assert Vy.shape == (25, 100)
        assert Vz.shape == (25, 100)

    def test_nan_propagation(self):
        """NaN propagates per axis: Vx←b1,b2; Vy←b3,b4; Vz←all four."""
        b0 = np.ones((3, 4))
        b_nan = b0.copy()
        b_nan[1, 2] = np.nan
        # NaN in b1 → Vx NaN, Vz NaN, Vy unaffected (depends only on b3, b4)
        Vx, Vy, Vz = beam2xyz(b_nan, b0, b0, b0, THETA)
        assert np.isnan(Vx[1, 2])
        assert np.isfinite(Vy[1, 2])  # b3, b4 are clean
        assert np.isnan(Vz[1, 2])
        # NaN in b3 → Vy NaN, Vz NaN, Vx unaffected
        Vx2, Vy2, Vz2 = beam2xyz(b0, b0, b_nan, b0, THETA)
        assert np.isfinite(Vx2[1, 2])
        assert np.isnan(Vy2[1, 2])
        assert np.isnan(Vz2[1, 2])
        # Other cells not affected
        assert np.isfinite(Vx[0, 0])

    def test_vx_from_b1_b2_only(self):
        """Vx depends only on b1 and b2; Vy depends only on b3 and b4.

        Down-looking (default) matrix per loadrdi.m b2earth: VX = +b1 - b2.
        """
        theta = np.radians(THETA)
        uvfac = 1.0 / (2.0 * np.sin(theta))
        b = np.zeros((1, 1))
        b1 = np.full((1, 1), 0.5)
        b2 = np.full((1, 1), -0.3)
        Vx, Vy, Vz = beam2xyz(b1, b2, b, b, THETA)
        expected_Vx = uvfac * (0.5 - (-0.3))
        assert abs(float(Vx[0, 0]) - expected_Vx) < 1e-10
        assert float(Vy[0, 0]) == 0.0

    def test_vx_up_looking_matrix(self):
        """Up-looking matrix per loadrdi.m b2earth: VX = -b1 + b2."""
        theta = np.radians(THETA)
        uvfac = 1.0 / (2.0 * np.sin(theta))
        b = np.zeros((1, 1))
        b1 = np.full((1, 1), 0.5)
        b2 = np.full((1, 1), -0.3)
        Vx, Vy, Vz = beam2xyz(b1, b2, b, b, THETA, beams_up=True)
        expected_Vx = uvfac * (-0.5 + (-0.3))
        assert abs(float(Vx[0, 0]) - expected_Vx) < 1e-10
        assert float(Vy[0, 0]) == 0.0

    def test_vy_from_b3_b4_only(self):
        theta = np.radians(THETA)
        uvfac = 1.0 / (2.0 * np.sin(theta))
        b = np.zeros((1, 1))
        b3 = np.full((1, 1), 0.2)
        b4 = np.full((1, 1), 0.6)
        Vx, Vy, Vz = beam2xyz(b, b, b3, b4, THETA)
        expected_Vy = uvfac * (-0.2 + 0.6)
        assert abs(float(Vy[0, 0]) - expected_Vy) < 1e-10
        assert float(Vx[0, 0]) == 0.0

    def test_vz_from_all_beams(self):
        """Down-looking: Vz = wfac * (+b1 + b2 + b3 + b4).

        Positive beam velocity = water toward the (down-facing) transducer
        = upward water motion, so Vz must come out positive.
        """
        theta = np.radians(THETA)
        wfac = 1.0 / (4.0 * np.cos(theta))
        b = np.full((1, 1), 0.1)
        Vx, Vy, Vz = beam2xyz(b, b, b, b, THETA)
        expected_Vz = wfac * (+0.4)
        assert abs(float(Vz[0, 0]) - expected_Vz) < 1e-10

    def test_vz_up_looking(self):
        """Up-looking: Vz = wfac * (-b1 - b2 - b3 - b4)."""
        theta = np.radians(THETA)
        wfac = 1.0 / (4.0 * np.cos(theta))
        b = np.full((1, 1), 0.1)
        Vx, Vy, Vz = beam2xyz(b, b, b, b, THETA, beams_up=True)
        expected_Vz = wfac * (-0.4)
        assert abs(float(Vz[0, 0]) - expected_Vz) < 1e-10

    def test_scalar_theta(self):
        """theta_deg can be a float scalar."""
        b = np.ones((2, 3))
        Vx, Vy, Vz = beam2xyz(b, -b, b, -b, 20.0)
        assert Vx.shape == (2, 3)


class TestBeam2earth:
    def test_identity_rotation_zero_heading_pitch_roll(self):
        """heading=pitch=roll=0 → rotation matrix is identity: u=Vx, v=Vy, w=Vz."""
        theta_deg = 20.0
        theta = np.radians(theta_deg)
        uvfac = 1.0 / (2.0 * np.sin(theta))
        wfac = 1.0 / (4.0 * np.cos(theta))

        nbin, nens = 5, 8
        b1 = np.random.default_rng(0).standard_normal((nbin, nens)) * 0.1
        b2 = np.random.default_rng(1).standard_normal((nbin, nens)) * 0.1
        b3 = np.random.default_rng(2).standard_normal((nbin, nens)) * 0.1
        b4 = np.random.default_rng(3).standard_normal((nbin, nens)) * 0.1

        heading = np.zeros(nens)
        pitch = np.zeros(nens)
        roll = np.zeros(nens)

        u, v, w = beam2earth(b1, b2, b3, b4, heading, pitch, roll, theta_deg)

        Vx_expected = uvfac * (+b1 - b2)
        Vy_expected = uvfac * (-b3 + b4)
        Vz_expected = wfac * (+b1 + b2 + b3 + b4)

        np.testing.assert_allclose(u, Vx_expected, rtol=1e-10)
        np.testing.assert_allclose(v, Vy_expected, rtol=1e-10)
        np.testing.assert_allclose(w, Vz_expected, rtol=1e-10)

    def test_output_shape(self):
        nbin, nens = 25, 100
        b = np.random.default_rng(42).standard_normal((nbin, nens))
        heading = np.random.default_rng(0).uniform(0, 360, nens)
        pitch = np.random.default_rng(1).uniform(-5, 5, nens)
        roll = np.random.default_rng(2).uniform(-5, 5, nens)
        u, v, w = beam2earth(b, b, b, b, heading, pitch, roll, 20.0)
        assert u.shape == (nbin, nens)
        assert v.shape == (nbin, nens)
        assert w.shape == (nbin, nens)

    def test_nan_propagation(self):
        """NaN beam → NaN Earth velocity for that cell."""
        nbin, nens = 3, 4
        b = np.ones((nbin, nens))
        b_nan = b.copy()
        b_nan[1, 2] = np.nan
        heading = np.full(nens, 90.0)
        pitch = np.zeros(nens)
        roll = np.zeros(nens)
        u, v, w = beam2earth(b_nan, b, b, b, heading, pitch, roll, 20.0)
        assert np.isnan(u[1, 2])
        assert np.isnan(v[1, 2])
        assert np.isnan(w[1, 2])
        assert np.isfinite(u[0, 0])

    def test_gimbaled_no_effect_when_pitch_roll_zero(self):
        """Gimbaled correction is zero when pitch=roll=0 (arcsin(0/sqrt(cos²+0))=0)."""
        nbin, nens = 3, 5
        b = np.random.default_rng(7).standard_normal((nbin, nens))
        heading = np.full(nens, 45.0)
        pitch = np.zeros(nens)
        roll = np.zeros(nens)
        u_g, v_g, w_g = beam2earth(
            b, b, b, b, heading, pitch, roll, 20.0, gimbaled=True
        )
        u_ng, v_ng, w_ng = beam2earth(
            b, b, b, b, heading, pitch, roll, 20.0, gimbaled=False
        )
        np.testing.assert_allclose(u_g, u_ng, rtol=1e-10)
        np.testing.assert_allclose(v_g, v_ng, rtol=1e-10)
        np.testing.assert_allclose(w_g, w_ng, rtol=1e-10)

    def test_all_nan_beams_give_nan_output(self):
        b = np.full((3, 4), np.nan)
        heading = np.zeros(4)
        pitch = np.zeros(4)
        roll = np.zeros(4)
        u, v, w = beam2earth(b, b, b, b, heading, pitch, roll, 20.0)
        assert np.all(np.isnan(u))
        assert np.all(np.isnan(v))
        assert np.all(np.isnan(w))


class TestUvrot:
    def test_zero_rotation_is_identity(self):
        u = np.array([1.0, 0.0, -0.5])
        v = np.array([0.0, 1.0,  0.3])
        ur, vr = uvrot(u, v, 0.0)
        np.testing.assert_allclose(ur, u, atol=1e-14)
        np.testing.assert_allclose(vr, v, atol=1e-14)

    def test_90_ccw_rotates_east_to_north(self):
        """Unit East vector (1, 0) rotated 90° CCW becomes North (0, 1)."""
        ur, vr = uvrot(np.array([1.0]), np.array([0.0]), 90.0)
        assert abs(float(ur[0])) < 1e-14
        assert abs(float(vr[0]) - 1.0) < 1e-14

    def test_180_rotates_east_to_west(self):
        ur, vr = uvrot(np.array([1.0]), np.array([0.0]), 180.0)
        assert abs(float(ur[0]) + 1.0) < 1e-14
        assert abs(float(vr[0])) < 1e-14

    def test_negative_90_cw(self):
        """Negative angle rotates clockwise: (1, 0) → (0, -1)."""
        ur, vr = uvrot(np.array([1.0]), np.array([0.0]), -90.0)
        assert abs(float(ur[0])) < 1e-14
        assert abs(float(vr[0]) + 1.0) < 1e-14

    def test_preserves_magnitude(self):
        """Rotation must not change vector magnitude."""
        rng = np.random.default_rng(99)
        u = rng.standard_normal(100)
        v = rng.standard_normal(100)
        ur, vr = uvrot(u, v, 17.5)
        mag_before = np.sqrt(u**2 + v**2)
        mag_after = np.sqrt(ur**2 + vr**2)
        np.testing.assert_allclose(mag_after, mag_before, rtol=1e-12)

    def test_nan_propagates(self):
        u = np.array([1.0, np.nan, 0.5])
        v = np.array([0.0, 1.0,   0.5])
        ur, vr = uvrot(u, v, 12.3)
        assert np.isnan(ur[1])
        assert np.isnan(vr[1])
        assert np.isfinite(ur[0])
        assert np.isfinite(ur[2])

    def test_2d_array_shape_preserved(self):
        u = np.ones((5, 10))
        v = np.zeros((5, 10))
        ur, vr = uvrot(u, v, 45.0)
        assert ur.shape == (5, 10)
        assert vr.shape == (5, 10)


# --- reconstruct_3beam (loadrdi.m b2earth lines 1713-1726) ---


class TestReconstruct3Beam:
    def _beams(self):
        # b1+b2-b3-b4 = 0 by construction (zero error velocity), so the
        # reconstruction is exact for these cells.
        b1 = np.array([[0.30, 0.10]])
        b2 = np.array([[0.20, 0.40]])
        b3 = np.array([[0.25, 0.35]])
        b4 = np.array([[0.25, 0.15]])
        return b1, b2, b3, b4

    def test_each_missing_beam_reconstructed_exactly(self):
        for missing in range(4):
            beams = [b.copy() for b in self._beams()]
            truth = beams[missing][0, 0]
            beams[missing][0, 0] = np.nan
            r, n3 = reconstruct_3beam(*beams)
            assert n3 == 1
            assert abs(r[missing][0, 0] - truth) < 1e-12
            # untouched cell intact
            assert r[missing][0, 1] == self._beams()[missing][0, 1]

    def test_two_missing_beams_stay_nan(self):
        b1, b2, b3, b4 = (b.copy() for b in self._beams())
        b1[0, 0] = np.nan
        b2[0, 0] = np.nan
        r, n3 = reconstruct_3beam(b1, b2, b3, b4)
        assert n3 == 0
        assert np.isnan(r[0][0, 0]) and np.isnan(r[1][0, 0])

    def test_all_finite_untouched(self):
        beams = self._beams()
        r, n3 = reconstruct_3beam(*beams)
        assert n3 == 0
        for orig, out in zip(beams, r):
            np.testing.assert_array_equal(orig, out)

    def test_inputs_not_mutated(self):
        beams = [b.copy() for b in self._beams()]
        beams[0][0, 0] = np.nan
        b1_orig = beams[0].copy()
        reconstruct_3beam(*beams)
        np.testing.assert_array_equal(beams[0], b1_orig)

    def test_reconstructed_cell_has_zero_error_velocity(self):
        # loadrdi.m asserts |VE| < 1e-9 for reconstructed cells.
        beams = [b.copy() for b in self._beams()]
        beams[2][0, 1] = np.nan
        r, _ = reconstruct_3beam(*beams)
        ve = r[0][0, 1] + r[1][0, 1] - r[2][0, 1] - r[3][0, 1]
        assert abs(ve) < 1e-12


class TestBeam2EarthAllow3Beam:
    def test_allow_3beam_fills_single_missing_beam(self):
        nbin, nens = 3, 5
        rng = np.random.default_rng(3)
        b1, b2, b3 = (rng.normal(0, 0.3, (nbin, nens)) for _ in range(3))
        b4 = b1 + b2 - b3  # zero error velocity everywhere
        h = np.zeros(nens)
        p = np.zeros(nens)
        r = np.zeros(nens)
        u_full, v_full, w_full = beam2earth(
            b1, b2, b3, b4, h, p, r, 20.0, gimbaled=False
        )
        b3_gap = b3.copy()
        b3_gap[1, 2] = np.nan
        u3, v3, w3 = beam2earth(
            b1, b2, b3_gap, b4, h, p, r, 20.0, gimbaled=False, allow_3beam=True
        )
        np.testing.assert_allclose(u3, u_full, atol=1e-12)
        np.testing.assert_allclose(w3, w_full, atol=1e-12)
        # without the flag the cell is NaN (existing behavior preserved)
        u0, _, _ = beam2earth(
            b1, b2, b3_gap, b4, h, p, r, 20.0, gimbaled=False
        )
        assert np.isnan(u0[1, 2])
