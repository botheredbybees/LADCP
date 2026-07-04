"""E1 diagnostic: fit the rotation angle mapping UL (u,v) onto DL (u,v).

Loads P16N cast 003 DL and UL PD0 files, runs each through the same
ingestion + beam2earth path as tests/integration/test_inverse_p16n_cast003.py
(each instrument with its OWN heading; UL with the current negate-pitch
convention), then compares the near-package water velocity seen by the two
instruments per ensemble (prepinv.m uses bins 1:4 of each instrument for the
same purpose, line ~329: iz=1:4).

Per ensemble: w = u + i*v averaged over the 4 near-package bins of each
instrument; the best-fit rotation mapping UL onto DL is
    theta = angle(w_dl * conj(w_ul))   [CCW positive]
Reported: circular mean, circular std, and correlation of theta against
package heading and cast phase (downcast/upcast).

Interpretation guide (VALIDATION_PLAN.md Phase 1 E1):
  ~0 deg           -> transform fine
  constant ~87 deg -> UL heading effectively unused/cancelled
  heading-dependent-> mirrored-frame (handedness) composition bug

Also runs a faithful port of loadrdi.m::b2earth (the LDEO reference transform,
UL beam-matrix sign flip + fixed-sensor pitch correction) on the same data as
a cross-check of what the correct convention produces.

Usage:  TEST_DATA_DIR=test_data uv run python scripts/diag_ul_dl_rotation.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ladcp.ingestion.ctd import assign_bin_depths, load_ctd  # noqa: E402
from ladcp.ingestion.rdi import load_rdi  # noqa: E402
from ladcp.transforms.beam2earth import beam2earth  # noqa: E402

THETA_DEG = 20.0
NREF_BINS = 4  # prepinv.m iz=1:4
MIN_SPEED = 0.05  # m/s floor: angle of a near-zero vector is noise


def ldeo_b2earth(b1, b2, b3, b4, heading, pitch, roll, theta_deg, beams_up):
    """Faithful vectorized port of loadrdi.m::b2earth (lines 1571-1767).

    Fixed-sensor attitude convention (Martini/Pluddeman):
        RR = roll;  PP = asin(sin(pitch)*cos(roll)/KA),
        KA = sqrt(1 - (sin(pitch)*sin(roll))^2)
    Beam->instrument matrix depends on beams_up (loadrdi lines 1736-1748).
    No 3-beam solutions (NaN in any beam -> NaN out), matching our Python path.
    """
    d2r = np.pi / 180.0
    RR = roll * d2r
    KA = np.sqrt(1.0 - (np.sin(pitch * d2r) * np.sin(roll * d2r)) ** 2)
    PP = np.arcsin(np.sin(pitch * d2r) * np.cos(roll * d2r) / KA)
    HH = heading * d2r

    CP, SP = np.cos(PP), np.sin(PP)
    CR, SR = np.cos(RR), np.sin(RR)
    CH, SH = np.cos(HH), np.sin(HH)

    S = np.sin(theta_deg * d2r)
    C = np.cos(theta_deg * d2r)
    VXS = VYS = 1.0 / (2.0 * S)
    VZS = 1.0 / (4.0 * C)

    if beams_up:
        VX = VXS * (-b1 + b2)
        VY = VYS * (-b3 + b4)
        VZ = VZS * (-b1 - b2 - b3 - b4)
    else:
        VX = VXS * (+b1 - b2)
        VY = VYS * (-b3 + b4)
        VZ = VZS * (+b1 + b2 + b3 + b4)

    u = VX * (CH * CR + SH * SR * SP) + VY * SH * CP + VZ * (CH * SR - SH * CR * SP)
    v = -VX * (SH * CR - CH * SR * SP) + VY * CH * CP - VZ * (SH * SR + CH * SP * CR)
    w = -VX * (SR * CP) + VY * SP + VZ * (CP * CR)
    return u, v, w


def circ_stats(theta_deg):
    """Circular mean and std (degrees) of an angle sample."""
    z = np.exp(1j * np.radians(theta_deg))
    zm = np.nanmean(z)
    mean = np.degrees(np.angle(zm))
    R = np.abs(zm)
    std = np.degrees(np.sqrt(-2.0 * np.log(max(R, 1e-12))))
    return mean, std


def near_package(u, v):
    """Mean complex velocity over the NREF_BINS bins nearest the transducer."""
    with np.errstate(invalid="ignore"):
        return np.nanmean(u[:NREF_BINS] + 1j * v[:NREF_BINS], axis=0)


def fit_rotation(u_dl, v_dl, u_ul_a, v_ul_a):
    """Per-ensemble rotation angle (deg, CCW) mapping UL (u,v) onto DL (u,v).

    Uses the mean over the NREF_BINS bins closest to the package of each
    instrument (rows 0..3 of each raw array = nearest the transducer).
    """
    wd = near_package(u_dl, v_dl)
    wu = near_package(u_ul_a, v_ul_a)
    ok = (
        np.isfinite(wd)
        & np.isfinite(wu)
        & (np.abs(wd) > MIN_SPEED)
        & (np.abs(wu) > MIN_SPEED)
    )
    theta = np.full(wd.shape, np.nan)
    theta[ok] = np.degrees(np.angle(wd[ok] * np.conj(wu[ok])))
    return theta, ok


def model_fits(a, b, label_a="UL", label_b="DL"):
    """Fit b ~ rotation(a) vs b ~ reflection(conj(a)); report which wins.

    Rotation model:   b = e^{i theta} a      -> theta = angle(sum b conj(a))
    Reflection model: b = e^{i alpha} conj(a)-> alpha = angle(sum b a)
    rho = |sum| / sqrt(sum|a|^2 sum|b|^2) is the model correlation (1 = exact).
    Magnitude-weighted, so noise from slow ensembles barely contributes.
    """
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    norm = np.sqrt(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2))
    s_rot = np.sum(b * np.conj(a))
    s_ref = np.sum(b * a)
    rho_rot = np.abs(s_rot) / norm
    rho_ref = np.abs(s_ref) / norm
    th_rot = np.degrees(np.angle(s_rot))
    th_ref = np.degrees(np.angle(s_ref))
    print(
        f"  model fit {label_a}->{label_b} (n={ok.sum()}): "
        f"ROTATION rho={rho_rot:.3f} theta={th_rot:+7.2f} deg | "
        f"REFLECTION rho={rho_ref:.3f} alpha={th_ref:+7.2f} deg"
    )
    return rho_rot, th_rot, rho_ref, th_ref


def report(tag, theta, ok, hdg_dl, downcast):
    mean, std = circ_stats(theta[ok])
    print(f"\n--- {tag} ---")
    print(f"  n ensembles used: {ok.sum()} / {len(ok)}")
    print(f"  rotation UL->DL: circ mean = {mean:+8.2f} deg, circ std = {std:6.2f} deg")
    for phase, sel in (("downcast", downcast & ok), ("upcast", (~downcast) & ok)):
        if sel.sum() > 10:
            m, s = circ_stats(theta[sel])
            print(f"    {phase:9s}: mean {m:+8.2f}, std {s:6.2f} deg (n={sel.sum()})")
    # correlation vs heading: center angle near its circular mean to unwrap
    th = (theta[ok] - mean + 180.0) % 360.0 - 180.0
    hd = hdg_dl[ok]
    r_h = np.corrcoef(th, hd)[0, 1]
    r_sin = np.corrcoef(th, np.sin(np.radians(hd)))[0, 1]
    r_cos = np.corrcoef(th, np.cos(np.radians(hd)))[0, 1]
    r_ph = np.corrcoef(th, downcast[ok].astype(float))[0, 1]
    print(f"  corr(theta, DL heading)      r = {r_h:+.3f}")
    print(f"  corr(theta, sin/cos heading) r = {r_sin:+.3f} / {r_cos:+.3f}")
    print(f"  corr(theta, cast phase)      r = {r_ph:+.3f}  (phase: downcast=1)")


def main():
    data_dir = Path(os.environ.get("TEST_DATA_DIR", "test_data")) / "2015_P16N"
    rdi = load_rdi(data_dir / "003DL000.000")
    rdi_ul = load_rdi(data_dir / "003UL000.000")
    ctd = load_ctd(data_dir / "003_01.cnv")

    print("=== Fixed-leader configuration ===")
    for name, r in (("DL", rdi), ("UL", rdi_ul)):
        print(
            f"  {name}: sysconfig=0x{r.sysconfig:04x}  beams_up={r.beams_up}"
            f"  beam_angle={r.beam_angle_deg} deg"
            f"  EA(hdg_align)={r.hdg_align_deg:+.2f} deg"
            f"  EB(hdg_bias)={r.hdg_bias_deg:+.2f} deg"
            f"  EX(coord)=0x{r.coord_transform:02x}"
        )

    print("\n=== Raw attitude statistics ===")
    for name, r in (("DL", rdi), ("UL", rdi_ul)):
        hm, hs = circ_stats(r.heading[np.isfinite(r.heading)])
        print(
            f"  {name}: heading circ-mean {hm:+7.2f} (std {hs:5.2f}); "
            f"pitch mean {np.nanmean(r.pitch):+6.2f}; "
            f"roll mean {np.nanmean(r.roll):+6.2f}; "
            f"roll median {np.nanmedian(r.roll):+6.2f}"
        )

    # Package depth and cast phase from CTD (as the integration test does)
    z_m, _ = assign_bin_depths(rdi, ctd, looker="down")
    i_deep = int(np.nanargmax(z_m))
    downcast = np.arange(rdi.nens) < i_deep
    print(f"\n  max package depth {np.nanmax(z_m):.0f} m at ens {i_deep}/{rdi.nens}")

    # Time-align UL to DL (same as integration test)
    ul_idx = np.argmin(
        np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0
    )

    # DL-UL raw compass offset
    dh = np.degrees(
        np.angle(np.exp(1j * np.radians(rdi.heading - rdi_ul.heading[ul_idx])))
    )
    m, s = circ_stats(dh[np.isfinite(dh)])
    md, _ = circ_stats(dh[downcast & np.isfinite(dh)])
    mu, _ = circ_stats(dh[(~downcast) & np.isfinite(dh)])
    print(
        f"  raw heading offset DL-UL: circ mean {m:+.2f} deg (std {s:.2f}); "
        f"downcast {md:+.2f}, upcast {mu:+.2f}"
    )

    # --- Path A: production beam2earth as called by the integration test ---
    # (post-fix convention: loadrdi fixed-sensor case, per-instrument beams_up,
    #  raw attitude for both instruments)
    u_dl, v_dl, _ = beam2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG, gimbaled=False, beams_up=False,
    )
    u_ul, v_ul, _ = beam2earth(
        rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
        rdi_ul.heading, rdi_ul.pitch, rdi_ul.roll,
        THETA_DEG, gimbaled=False, beams_up=True,
    )
    theta, ok = fit_rotation(u_dl, v_dl, u_ul[:, ul_idx], v_ul[:, ul_idx])
    report("Path A: production beam2earth (gimbaled=False, beams_up per inst)",
           theta, ok, rdi.heading, downcast)
    wu_a_pkg = near_package(u_ul[:, ul_idx], v_ul[:, ul_idx])
    wd_a_pkg = near_package(u_dl, v_dl)
    model_fits(wu_a_pkg, wd_a_pkg, "UL(A)", "DL(A)")
    model_fits(wu_a_pkg[downcast], wd_a_pkg[downcast], "UL(A) down", "DL(A) down")
    model_fits(wu_a_pkg[~downcast], wd_a_pkg[~downcast], "UL(A) up  ", "DL(A) up  ")

    # --- Path B: faithful loadrdi.m b2earth port ---
    u_dl2, v_dl2, _ = ldeo_b2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG, beams_up=False,
    )
    u_ul2, v_ul2, _ = ldeo_b2earth(
        rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
        rdi_ul.heading, rdi_ul.pitch, rdi_ul.roll,
        THETA_DEG, beams_up=True,
    )
    theta2, ok2 = fit_rotation(u_dl2, v_dl2, u_ul2[:, ul_idx], v_ul2[:, ul_idx])
    report("Path B: loadrdi.m b2earth port (DL down-matrix, UL up-matrix)",
           theta2, ok2, rdi.heading, downcast)
    wu_b_pkg = near_package(u_ul2[:, ul_idx], v_ul2[:, ul_idx])
    wd_b_pkg = near_package(u_dl2, v_dl2)
    model_fits(wu_b_pkg, wd_b_pkg, "UL(B)", "DL(B)")
    model_fits(wu_b_pkg[downcast], wd_b_pkg[downcast], "UL(B) down", "DL(B) down")
    model_fits(wu_b_pkg[~downcast], wd_b_pkg[~downcast], "UL(B) up  ", "DL(B) up  ")

    # --- Cross-checks: same instrument, Python path vs loadrdi port ---
    print("\n--- Cross-check: our DL vs loadrdi DL (same beams, same attitude) ---")
    model_fits(wd_a_pkg, wd_b_pkg, "DL(py)", "DL(ldeo)")
    print("--- Cross-check: our UL vs loadrdi UL (same beams, same attitude) ---")
    model_fits(wu_a_pkg, wu_b_pkg, "UL(py)", "UL(ldeo)")

    # --- w sanity: during downcast, water moves up relative to the package ---
    _, _, w_dl_py = beam2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG, gimbaled=False, beams_up=False,
    )
    _, _, w_dl_ld = ldeo_b2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG, beams_up=False,
    )
    wm_py = np.nanmean(w_dl_py[:NREF_BINS][:, downcast])
    wm_ld = np.nanmean(w_dl_ld[:NREF_BINS][:, downcast])
    print("\n--- w sanity check (downcast mean near-package w; expect > 0:")
    print("    package sinks ~1 m/s, so water moves UP relative to package) ---")
    print(f"  our beam2earth DL w = {wm_py:+.3f}; loadrdi port DL w = {wm_ld:+.3f} m/s")


if __name__ == "__main__":
    main()
