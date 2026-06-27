"""Diagnose n_se discrepancy: 524 (Python) vs 827 (MATLAB).

Investigates two hypotheses:
  A. dz mismatch: test uses dz=16.0, MATLAB uses median bin spacing (~8m)
  B. depth variable: Python uses ens.z (CTD depth), MATLAB uses d.izm(1,:) (top bin depth)

Run with:
    TEST_DATA_DIR=test_data python scripts/diagnose_n_se.py
"""
from __future__ import annotations
import os, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ladcp.ingestion.rdi import load_rdi
from ladcp.ingestion.ctd import assign_bin_depths, load_ctd
from ladcp.transforms.beam2earth import beam2earth, uvrot
from ladcp.qa.editing import edit_sidelobes, edit_large_velocities, edit_w_outliers
from ladcp.solution.inverse import EnsembleData, prepare_superensembles, _window_boundaries

TEST_DATA_DIR = os.environ.get("TEST_DATA_DIR", "test_data")
P16N = Path(TEST_DATA_DIR) / "2015_P16N"
THETA_DEG = 20.0
DROT_DEG = 12.318441

def main():
    dl_path = P16N / "003DL000.000"
    ul_path = P16N / "003UL000.000"
    cnv_path = P16N / "003_01.cnv"

    if not dl_path.exists():
        sys.exit(f"DL file not found: {dl_path}")

    print("=== Loading DL PD0 ===")
    rdi = load_rdi(dl_path)
    print(f"  n_ens={rdi.nens}  n_bins={rdi.nbin}")
    print(f"  blen_m (bin size) = {rdi.blen_m:.2f} m   <- MATLAB avdz default")
    print(f"  blnk_m (blanking) = {rdi.blnk_m:.2f} m")

    ctd = load_ctd(cnv_path)
    print(f"\n=== CTD ===")
    print(f"  n_samples={len(ctd.time_julian)}")

    # Build combined DL+UL array (same as test fixture)
    u_dl, v_dl, w_dl = beam2earth(rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll, THETA_DEG, gimbaled=True)
    u_dl, v_dl = uvrot(u_dl, v_dl, -DROT_DEG)
    z_m, izm_dl_pos = assign_bin_depths(rdi, ctd, looker="down")
    z_neg = -z_m
    weight_dl = np.nanmean(rdi.corr.astype(np.float64), axis=2) / 128.0

    rdi_ul = load_rdi(ul_path)
    u_ul, v_ul, w_ul = beam2earth(rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
        rdi_ul.heading, -rdi_ul.pitch, rdi_ul.roll, THETA_DEG, gimbaled=True)
    u_ul, v_ul = uvrot(u_ul, v_ul, -DROT_DEG)
    _, izm_ul_pos = assign_bin_depths(rdi_ul, ctd, looker="up")
    weight_ul = np.nanmean(rdi_ul.corr.astype(np.float64), axis=2) / 128.0

    ul_idx = np.argmin(
        np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0)
    u_ul_a = u_ul[:, ul_idx]
    v_ul_a = v_ul[:, ul_idx]
    w_ul_a = w_ul[:, ul_idx]
    weight_ul_a = weight_ul[:, ul_idx]
    izm_ul_neg_a = -izm_ul_pos[:, ul_idx]

    n_ul = rdi_ul.nbin
    n_dl = rdi.nbin
    u_comb   = np.vstack([u_ul_a[::-1, :],    u_dl])
    v_comb   = np.vstack([v_ul_a[::-1, :],    v_dl])
    w_comb   = np.vstack([w_ul_a[::-1, :],    w_dl])
    weight_c = np.vstack([weight_ul_a[::-1,:], weight_dl])
    izm_comb = np.vstack([izm_ul_neg_a[::-1,:],-izm_dl_pos])

    izu = np.arange(n_ul - 1, -1, -1, dtype=int)
    izd = np.arange(n_ul, n_ul + n_dl, dtype=int)

    bt_u_e, bt_v_e, bt_w_e = beam2earth(
        rdi.btrack_vel_ms[0], rdi.btrack_vel_ms[1],
        rdi.btrack_vel_ms[2], rdi.btrack_vel_ms[3],
        rdi.heading, rdi.pitch, rdi.roll, THETA_DEG, gimbaled=True)
    bt_u_e, bt_v_e = uvrot(bt_u_e, bt_v_e, -DROT_DEG)
    bvel = np.stack([bt_u_e, bt_v_e, bt_w_e], axis=1)
    bvels = np.full_like(bvel, 0.02)
    hbot = np.nanmean(rdi.btrack_range_m, axis=0)

    slat = np.full(rdi.nens, np.nan)
    slon = np.full(rdi.nens, np.nan)

    ens = EnsembleData(
        u=u_comb, v=v_comb, w=w_comb, weight=weight_c,
        izm=izm_comb, z=z_neg,
        time_jul=rdi.time_julian,
        bvel=bvel, bvels=bvels, hbot=hbot,
        izd=izd, izu=izu,
        slat=slat, slon=slon,
    )

    ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)
    ens = edit_large_velocities(ens)
    ens = edit_w_outliers(ens)

    print(f"\n=== Depth variable comparison ===")
    print(f"  ens.z  (CTD depth):      min={ens.z.min():.1f}  max={ens.z.max():.1f}  n={len(ens.z)}")
    top_bin_depth = ens.izm[0, :]  # shallowest combined bin (= MATLAB d.izm(1,:))
    print(f"  izm[0] (top bin depth):  min={np.nanmin(top_bin_depth):.1f}  max={np.nanmax(top_bin_depth):.1f}")
    print(f"  Rate of change per ens:")
    print(f"    d(ens.z)/dt:  median abs diff = {np.nanmedian(np.abs(np.diff(ens.z))):.4f} m/ens")
    print(f"    d(izm[0])/dt: median abs diff = {np.nanmedian(np.abs(np.diff(top_bin_depth))):.4f} m/ens")

    print(f"\n=== n_se at different dz values ===")
    print(f"  MATLAB default avdz = blen_m = {rdi.blen_m:.1f} m")
    for dz_test in [rdi.blen_m, rdi.blen_m * 2, 8.0, 10.0, 16.0, None]:
        # Hypothesis A: using ens.z with different dz
        if dz_test is None:
            dz_val = float(np.nanmedian(np.abs(np.diff(ens.izm[:, 0]))))
            label = f"None → {dz_val:.1f}m (auto from izm[:,0])"
        else:
            dz_val = float(dz_test)
            label = f"{dz_val:.1f}m"
        wins_z = _window_boundaries(ens.z, dz_val)
        print(f"  dz={label}  ->  n_se (using ens.z)  = {len(wins_z)}")

    print(f"\n  --- Hypothesis B: use izm[0,:] instead of ens.z ---")
    for dz_test in [rdi.blen_m, rdi.blen_m * 2, 8.0, 16.0]:
        dz_val = float(dz_test)
        # Replace NaN in top_bin_depth with interpolated values
        izm0 = top_bin_depth.copy()
        valid = np.isfinite(izm0)
        if valid.any():
            izm0 = np.interp(np.arange(len(izm0)), np.where(valid)[0], izm0[valid])
        wins_izm = _window_boundaries(izm0, dz_val)
        print(f"  dz={dz_val:.1f}m  ->  n_se (using izm[0,:]) = {len(wins_izm)}")

    print(f"\n  MATLAB reference: n_se = 827")
    print(f"  Python current:   n_se = 524 (dz=16.0, ens.z)")

    # Show actual auto-computed dz from izm[:,0]
    auto_dz = float(np.nanmedian(np.abs(np.diff(ens.izm[:, 0]))))
    print(f"\n  Auto dz from izm[:,0] = {auto_dz:.3f} m  (bin spacing in first ensemble)")
    auto_dz_z = float(np.nanmedian(np.abs(np.diff(ens.izm[:, 0]))))  # same
    print(f"  rdi.blen_m = {rdi.blen_m:.3f} m  (from PD0 fixed leader)")

    print(f"\n=== Best match attempt ===")
    se_best = prepare_superensembles(ens, dz=rdi.blen_m)
    print(f"  prepare_superensembles(ens, dz=blen_m={rdi.blen_m:.1f}) → n_se={se_best.izm.shape[1]}")


if __name__ == "__main__":
    main()
