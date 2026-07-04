"""Diagnostic: compare su-ensemble ru values against reference u profile.

If ru[k, m] = u_ocean(z_k) - u_ctd(m), and u_ctd is small (~0.05 m/s),
then ru should approximate u_ocean.  Anti-correlation in ru vs reference
means the bug is in prepare_superensembles; otherwise it's in compute_inverse.

Run from repo root:
    TEST_DATA_DIR=test_data uv run python scripts/diag_ru_vs_ref.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ladcp.ingestion.ctd import assign_bin_depths, load_ctd
from ladcp.ingestion.rdi import load_rdi
from ladcp.qa.editing import edit_large_velocities, edit_sidelobes, edit_w_outliers
from ladcp.solution.inverse import (
    EnsembleData,
    InverseParams,
    _build_ctd_matrix,
    _build_obs_matrix,
    _flatten_obs,
    compute_inverse,
    prepare_superensembles,
)
from ladcp.transforms.beam2earth import beam2earth, uvrot

THETA_DEG = 20.0
DROT_DEG = 12.318441


def load_data():
    env = os.environ.get("TEST_DATA_DIR", "")
    data_dir = Path(env) / "2015_P16N"

    rdi = load_rdi(data_dir / "003DL000.000")
    ctd = load_ctd(data_dir / "003_01.cnv")
    rdi_ul = load_rdi(data_dir / "003UL000.000")

    ref_ds = netCDF4.Dataset(data_dir / "003.nc")
    ref_z = np.array(ref_ds.variables["z"][:])
    ref_u = np.array(ref_ds.variables["u"][:])
    ref_nvel = np.array(ref_ds.variables["nvel"][:])
    u_ship = float(ref_ds.uship)
    v_ship = float(ref_ds.vship)
    ref_ds.close()

    u_dl, v_dl, w_dl = beam2earth(rdi.u, rdi.v, rdi.w, rdi.e,
                                   rdi.heading, rdi.pitch, rdi.roll,
                                   THETA_DEG, gimbaled=True)
    u_dl, v_dl = uvrot(u_dl, v_dl, -DROT_DEG)
    z_m, izm_dl_pos = assign_bin_depths(rdi, ctd, looker="down")
    z_neg = -z_m
    weight_dl = np.nanmean(rdi.corr.astype(np.float64), axis=2) / 128.0

    u_ul, v_ul, w_ul = beam2earth(rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
                                   rdi_ul.heading, -rdi_ul.pitch, rdi_ul.roll,
                                   THETA_DEG, gimbaled=True)
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
    u_comb = np.vstack([u_ul_a[::-1, :], u_dl])
    v_comb = np.vstack([v_ul_a[::-1, :], v_dl])
    w_comb = np.vstack([w_ul_a[::-1, :], w_dl])
    weight_comb = np.vstack([weight_ul_a[::-1, :], weight_dl])
    izm_comb = np.vstack([izm_ul_neg_a[::-1, :], -izm_dl_pos])

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

    sadcp_path = data_dir / "sadcp_003.npz"
    npz = np.load(sadcp_path)
    sadcp_z, sadcp_u, sadcp_v, sadcp_err = npz["z"], npz["u"], npz["v"], npz["err"]

    ens = EnsembleData(
        u=u_comb, v=v_comb, w=w_comb, weight=weight_comb,
        izm=izm_comb, z=z_neg,
        time_jul=rdi.time_julian,
        bvel=bvel, bvels=bvels, hbot=hbot,
        izd=izd, izu=izu,
        slat=slat, slon=slon,
    )
    ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)
    ens = edit_large_velocities(ens)
    ens = edit_w_outliers(ens)

    return ens, ref_z, ref_u, ref_nvel, u_ship, v_ship, sadcp_z, sadcp_u, sadcp_v, sadcp_err


def run():
    print("Loading data...")
    ens, ref_z, ref_u, ref_nvel, u_ship, v_ship, sadcp_z, sadcp_u, sadcp_v, sadcp_err = load_data()

    print("Building superensembles...")
    se = prepare_superensembles(ens)
    print(f"  n_se={se.ru.shape[1]}, n_bins={se.ru.shape[0]}")
    print(f"  SE depth range: {-se.z.max():.0f} to {-se.z.min():.0f} m (CTD depth)")

    # === Strategy 1: look at mean ru at each depth ===
    # For each SE m, mean ru across DL bins (rows n_ul..n_ul+n_dl-1)
    # gives the "average ADCP velocity in the water column below the instrument"
    # which should loosely approximate u_ocean(z_dl_mean) - u_ctd(m)
    print("\n=== Strategy 1: ru averaged across DL bins vs reference u ===")
    print("  (high u_ctd means ru can differ from u_ocean; focus on SIGN pattern)")

    n_bins, n_se = se.ru.shape
    n_ul = ens.izu.max() + 1
    n_dl_bins = n_bins - n_ul

    # Downcast SEs: z_ctd increasing (going deeper)
    # Find SEs in key depth ranges
    for band_lo, band_hi in [(900, 1100), (1000, 1200), (1200, 1500), (1500, 2000)]:
        # SE depth: se.z is NEGATIVE (CTD convention), so deeper = more negative
        in_band = (-se.z >= band_lo) & (-se.z < band_hi)
        if in_band.sum() < 3:
            continue
        band_idx = np.where(in_band)[0]
        # Mean ru across DL bins (approximate u_ocean at DL mean depth - u_ctd)
        ru_dl_mean = np.nanmean(se.ru[n_ul:, in_band], axis=0)
        izm_dl_mean = np.nanmean(-se.izm[n_ul:, in_band], axis=0)
        # Reference u at those depths
        ref_u_interp = np.interp(izm_dl_mean, ref_z, ref_u,
                                  left=np.nan, right=np.nan)
        valid = np.isfinite(ru_dl_mean) & np.isfinite(ref_u_interp)
        if valid.sum() < 2:
            continue
        corr = np.corrcoef(ru_dl_mean[valid], ref_u_interp[valid])[0, 1]
        print(f"  CTD {band_lo}-{band_hi}m: n_se={in_band.sum()}, "
              f"ru_dl_mean={np.nanmean(ru_dl_mean):.4f} m/s, "
              f"ref_u_at_dlmean={np.nanmean(ref_u_interp):.4f} m/s, "
              f"corr(ru, ref_u)={corr:.3f}")

    # === Strategy 2: per-bin ru depth profile averaged over 1000-2000m SEs ===
    print("\n=== Strategy 2: bin-by-bin ru profile for SEs near 1500m CTD depth ===")
    # Find SEs where CTD is near 1500m (should have DL bins at 1500-1700m)
    target_z = 1500.0
    se_target = np.argmin(np.abs(-se.z - target_z))
    print(f"  Closest SE to z=1500m: SE{se_target}, CTD z={-se.z[se_target]:.1f}m")
    print(f"  {'bin':>4}  {'izm[m]':>8}  {'ru[m/s]':>10}  {'ref_u[m/s]':>12}  {'diff':>8}")
    for k in range(n_bins):
        iz = float(-se.izm[k, se_target])
        ru_k = float(se.ru[k, se_target])
        ref_u_k = float(np.interp(iz, ref_z, ref_u)) if np.isfinite(iz) else np.nan
        if np.isfinite(ru_k) and iz > 0:
            flag = "DL" if k >= n_ul else "UL"
            print(f"  {k:3d}{flag}  {iz:8.1f}  {ru_k:+10.4f}  {ref_u_k:+12.4f}  {ru_k-ref_u_k:+8.4f}")

    # === Strategy 3: Full inversion and correlation by band ===
    print("\n=== Strategy 3: Full inversion vs reference correlation by depth band ===")
    result = compute_inverse(se, u_ship=u_ship, v_ship=v_ship,
                              sadcp_z=sadcp_z, sadcp_u=sadcp_u,
                              sadcp_v=sadcp_v, sadcp_err=sadcp_err)
    for band_lo, band_hi in [(0,500),(500,1000),(1000,1500),(1500,2000),(2000,2500),(2500,3000),(3000,4400)]:
        mask = (ref_z >= band_lo) & (ref_z < band_hi) & np.isfinite(ref_u) & (ref_nvel >= 3)
        if mask.sum() < 3:
            continue
        py_u = np.interp(ref_z[mask], result.z, result.u)
        valid = np.isfinite(py_u)
        if valid.sum() < 3:
            continue
        corr = np.corrcoef(ref_u[mask][valid], py_u[valid])[0, 1]
        rmse = np.sqrt(np.mean((py_u[valid] - ref_u[mask][valid])**2))
        print(f"  {band_lo}-{band_hi}m: corr={corr:+.3f}, rmse={rmse:.4f} m/s")

    # === Strategy 4: Check if the vertical GRADIENT of ru is right ===
    # Compare ru shear vs reference u shear
    print("\n=== Strategy 4: ru shear vs reference u shear at 1500m SE ===")
    print("  Shear = ru[k+1] - ru[k] per 8m bin")
    print(f"  {'z_mid':>8}  {'ru_shear':>12}  {'ref_shear':>12}")
    for k in range(1, n_bins - 1):
        iz0 = float(-se.izm[k-1, se_target])
        iz1 = float(-se.izm[k, se_target])
        if not (np.isfinite(iz0) and np.isfinite(iz1) and iz0 > 0 and iz1 > 0):
            continue
        ru_sh = float(se.ru[k, se_target] - se.ru[k-1, se_target])
        ref0 = float(np.interp(iz0, ref_z, ref_u))
        ref1 = float(np.interp(iz1, ref_z, ref_u))
        ref_sh = ref1 - ref0
        if np.isfinite(ru_sh) and np.isfinite(ref_sh):
            flag = "DL" if k >= n_ul else "UL"
            print(f"  {(iz0+iz1)/2:8.1f}m{flag}  {ru_sh:+12.5f}  {ref_sh:+12.5f}")

    print("\nDone.")


if __name__ == "__main__":
    run()
