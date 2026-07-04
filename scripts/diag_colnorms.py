"""Diagnostic: check A_ocean column norms to identify ill-constrained depth bins.

Run from repo root:
    TEST_DATA_DIR=test_data uv run python scripts/diag_colnorms.py

Checks:
1. Which depth bins are "ill-constrained" (sqrt(col_sum) < 0.3 * median) — these are
   the bins where MATLAB's lainsmoo adds regularization even at smoofac=0.
2. Where those bins fall in the depth profile relative to the 1000-2000m anti-correlation.
3. Shape of anti-correlation: mirror vs wiggly.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ladcp.ingestion.ctd import assign_bin_depths, compute_ship_velocity, load_ctd
from ladcp.ingestion.rdi import load_rdi
from ladcp.qa.editing import edit_large_velocities, edit_sidelobes, edit_w_outliers
from ladcp.solution.inverse import (
    EnsembleData,
    InverseParams,
    _apply_weights,
    _build_ctd_matrix,
    _build_obs_matrix,
    _flatten_obs,
    compute_inverse,
    prepare_superensembles,
)
from ladcp.transforms.beam2earth import beam2earth, uvrot

THETA_DEG = 20.0
DROT_DEG = 12.318441

def run():
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        print("ERROR: TEST_DATA_DIR not set")
        sys.exit(1)
    data_dir = Path(env) / "2015_P16N"

    # --- Load raw data (same as integration test) ---
    rdi = load_rdi(data_dir / "003DL000.000")
    ctd = load_ctd(data_dir / "003_01.cnv")
    rdi_ul = load_rdi(data_dir / "003UL000.000")
    ref_ds = netCDF4.Dataset(data_dir / "003.nc")
    ref_z   = np.array(ref_ds.variables["z"][:])
    ref_u   = np.array(ref_ds.variables["u"][:])
    ref_nvel = np.array(ref_ds.variables["nvel"][:])
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
    try:
        _ds = netCDF4.Dataset(data_dir / "003.nc")
        u_ship = float(_ds.uship)
        v_ship = float(_ds.vship)
        _ds.close()
    except Exception:
        u_ship, v_ship = None, None

    sadcp_path = data_dir / "sadcp_003.npz"
    if sadcp_path.exists():
        npz = np.load(sadcp_path)
        sadcp_z, sadcp_u, sadcp_v, sadcp_err = npz["z"], npz["u"], npz["v"], npz["err"]
    else:
        sadcp_z = sadcp_u = sadcp_v = sadcp_err = None

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

    se = prepare_superensembles(ens)

    # --- Full inversion for profile comparison ---
    result = compute_inverse(
        se,
        u_ship=u_ship, v_ship=v_ship,
        sadcp_z=sadcp_z, sadcp_u=sadcp_u,
        sadcp_v=sadcp_v, sadcp_err=sadcp_err,
    )

    # --- Shape of anti-correlation ---
    print("\n=== Profile comparison at 1000-2000m ===")
    for band_lo, band_hi in [(0, 500), (500, 1000), (1000, 1500), (1500, 2000)]:
        mask = (ref_z >= band_lo) & (ref_z < band_hi) & np.isfinite(ref_u) & (ref_nvel >= 3)
        if mask.sum() < 3:
            continue
        py_u = np.interp(ref_z[mask], result.z, result.u)
        valid = np.isfinite(py_u)
        if valid.sum() < 3:
            continue
        corr = np.corrcoef(ref_u[mask][valid], py_u[valid])[0, 1]
        rmse = np.sqrt(np.mean((py_u[valid] - ref_u[mask][valid]) ** 2))
        py_mean = np.mean(py_u[valid])
        ref_mean = np.mean(ref_u[mask][valid])
        print(f"  {band_lo}-{band_hi}m: corr={corr:.3f}, rmse={rmse:.4f}, py_mean={py_mean:.3f}, ref_mean={ref_mean:.3f}")

    # --- Quick visual of profile shape at depth ---
    print("\n=== Python vs Reference u [m/s] at selected depths ===")
    print(f"{'depth':>6}  {'python':>8}  {'ref':>8}  {'diff':>8}")
    for target_z in range(1000, 2100, 100):
        py_u = float(np.interp(target_z, result.z, result.u))
        ref_idx = np.argmin(np.abs(ref_z - target_z))
        ru = ref_u[ref_idx]
        print(f"  {target_z:4d}m  {py_u:+.4f}  {ru:+.4f}  {(py_u-ru):+.4f}")

    # --- A_ocean column norm diagnostics ---
    print("\n=== A_ocean column norm analysis ===")
    params = InverseParams()  # dz=10.0 default
    d_u, d_v, izv, jprof, wm = _flatten_obs(se, params.velerr, params.weightmin)
    n_se = se.izm.shape[1]
    A_ocean_sp = _build_obs_matrix(izv, params.dz)
    A_ctd_sp = _build_ctd_matrix(jprof, n_se)
    n_zbins = A_ocean_sp.shape[1]
    A_o_u, A_c_u, dw_u, idx_down, idx_up = _apply_weights(A_ocean_sp, A_ctd_sp, d_u, wm)

    print(f"  dz used: {params.dz} m")
    print(f"  n_zbins (A_ocean cols): {n_zbins}")
    print(f"  max depth of A_ocean: {n_zbins * params.dz:.0f} m")

    # Column norms (sqrt of column sum, like MATLAB fs=sqrt(full(sum(A))))
    col_norms = np.sqrt(np.abs(A_o_u).sum(axis=0))
    median_norm = max(float(np.median(col_norms[col_norms > 0])), 0.01)
    ill_mask = col_norms < 0.3 * median_norm

    print(f"  median col_norm: {median_norm:.4f}")
    print(f"  threshold (0.3*median): {0.3 * median_norm:.4f}")
    print(f"  n ill-constrained bins: {ill_mask.sum()} / {n_zbins}")

    if ill_mask.any():
        ill_depths = (np.where(ill_mask)[0] + 1) * params.dz
        print(f"  ill-constrained depth range: {ill_depths.min():.0f} – {ill_depths.max():.0f} m")
        print(f"  ill-constrained depths (first 20): {ill_depths[:20]}")
    else:
        print("  No ill-constrained bins found (all within 0.3*median threshold)")

    # Also show column norm vs depth profile in key ranges
    print("\n  col_norm by depth range:")
    for lo, hi in [(0, 500), (500, 1000), (1000, 1500), (1500, 2000)]:
        lo_idx = int(lo / params.dz)
        hi_idx = min(int(hi / params.dz), n_zbins)
        if hi_idx > lo_idx:
            norms_in_range = col_norms[lo_idx:hi_idx]
            n_zero = (norms_in_range == 0).sum()
            n_ill = (norms_in_range < 0.3 * median_norm).sum()
            print(f"    {lo}-{hi}m: mean_norm={np.mean(norms_in_range):.4f}, n_zero={n_zero}, n_ill={n_ill}/{len(norms_in_range)}")

    # Also check with dz=8m (what MATLAB uses)
    print("\n=== With dz=8m (MATLAB equivalent) ===")
    dz_matlab = float(np.nanmedian(np.abs(np.diff(se.izm[:, 0]))))
    print(f"  Computed dz from se.izm: {dz_matlab:.2f} m")
    A_ocean_8 = _build_obs_matrix(izv, dz_matlab)
    n_zbins_8 = A_ocean_8.shape[1]
    A_o_8, _, _, _, _ = _apply_weights(A_ocean_8, A_ctd_sp, d_u, wm)
    col_norms_8 = np.sqrt(np.abs(A_o_8).sum(axis=0))
    median_norm_8 = max(float(np.median(col_norms_8[col_norms_8 > 0])), 0.01)
    ill_mask_8 = col_norms_8 < 0.3 * median_norm_8
    print(f"  n_zbins (8m): {n_zbins_8}")
    print(f"  n ill-constrained (8m): {ill_mask_8.sum()} / {n_zbins_8}")
    if ill_mask_8.any():
        ill_depths_8 = (np.where(ill_mask_8)[0] + 1) * dz_matlab
        print(f"  ill-constrained depth range (8m): {ill_depths_8.min():.0f} – {ill_depths_8.max():.0f} m")

    print("\nDone.")


if __name__ == "__main__":
    run()
