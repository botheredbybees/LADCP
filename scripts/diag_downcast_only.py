"""Diagnostic: downcast-only vs full-cast inversion to isolate upcast contamination."""
from __future__ import annotations
import sys, dataclasses
from pathlib import Path
import netCDF4, numpy as np
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from ladcp.ingestion.ctd import assign_bin_depths, load_ctd
from ladcp.ingestion.rdi import load_rdi
from ladcp.qa.editing import edit_large_velocities, edit_sidelobes, edit_w_outliers
from ladcp.solution.inverse import EnsembleData, compute_inverse, prepare_superensembles
from ladcp.transforms.beam2earth import beam2earth, uvrot

THETA_DEG = 20.0; DROT_DEG = 12.318441
data_dir = Path("test_data/2015_P16N")
rdi = load_rdi(data_dir / "003DL000.000"); ctd = load_ctd(data_dir / "003_01.cnv")
rdi_ul = load_rdi(data_dir / "003UL000.000")
ref_ds = netCDF4.Dataset(data_dir / "003.nc")
ref_z = np.array(ref_ds.variables["z"][:]); ref_u = np.array(ref_ds.variables["u"][:])
ref_v = np.array(ref_ds.variables["v"][:]); ref_nvel = np.array(ref_ds.variables["nvel"][:])
ref_ds.close()
npz = np.load(data_dir / "sadcp_003.npz")

u_dl, v_dl, w_dl = beam2earth(rdi.u, rdi.v, rdi.w, rdi.e, rdi.heading, rdi.pitch, rdi.roll, THETA_DEG, gimbaled=True)
u_dl, v_dl = uvrot(u_dl, v_dl, -DROT_DEG)
z_m, izm_dl_pos = assign_bin_depths(rdi, ctd, looker="down")
z_neg = -z_m; weight_dl = np.nanmean(rdi.corr.astype(np.float64), axis=2) / 128.0
u_ul, v_ul, w_ul = beam2earth(rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e, rdi_ul.heading, -rdi_ul.pitch, rdi_ul.roll, THETA_DEG, gimbaled=True)
u_ul, v_ul = uvrot(u_ul, v_ul, -DROT_DEG)
_, izm_ul_pos = assign_bin_depths(rdi_ul, ctd, looker="up")
weight_ul = np.nanmean(rdi_ul.corr.astype(np.float64), axis=2) / 128.0
ul_idx = np.argmin(np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0)
u_ul_a = u_ul[:, ul_idx]; v_ul_a = v_ul[:, ul_idx]; w_ul_a = w_ul[:, ul_idx]
weight_ul_a = weight_ul[:, ul_idx]; izm_ul_neg_a = -izm_ul_pos[:, ul_idx]
n_ul = rdi_ul.nbin; n_dl = rdi.nbin
u_comb = np.vstack([u_ul_a[::-1, :], u_dl]); v_comb = np.vstack([v_ul_a[::-1, :], v_dl])
w_comb = np.vstack([w_ul_a[::-1, :], w_dl]); weight_comb = np.vstack([weight_ul_a[::-1, :], weight_dl])
izm_comb = np.vstack([izm_ul_neg_a[::-1, :], -izm_dl_pos])
izu = np.arange(n_ul - 1, -1, -1, dtype=int); izd = np.arange(n_ul, n_ul + n_dl, dtype=int)
bt_u_e, bt_v_e, bt_w_e = beam2earth(rdi.btrack_vel_ms[0], rdi.btrack_vel_ms[1], rdi.btrack_vel_ms[2], rdi.btrack_vel_ms[3], rdi.heading, rdi.pitch, rdi.roll, THETA_DEG, gimbaled=True)
bt_u_e, bt_v_e = uvrot(bt_u_e, bt_v_e, -DROT_DEG)
bvel = np.stack([bt_u_e, bt_v_e, bt_w_e], axis=1); bvels = np.full_like(bvel, 0.02)
hbot = np.nanmean(rdi.btrack_range_m, axis=0)
ens = EnsembleData(u=u_comb, v=v_comb, w=w_comb, weight=weight_comb, izm=izm_comb, z=z_neg,
    time_jul=rdi.time_julian, bvel=bvel, bvels=bvels, hbot=hbot, izd=izd, izu=izu,
    slat=np.full(rdi.nens, np.nan), slon=np.full(rdi.nens, np.nan))
ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)
ens = edit_large_velocities(ens); ens = edit_w_outliers(ens)
ens.weight[n_ul - 1, :] = np.nan; ens.weight[n_ul, :] = np.nan; ens.weight[n_ul + 6:, :] = np.nan

se = prepare_superensembles(ens)
n_se = se.ru.shape[1]
turnaround = int(np.argmin(se.z))
print(f"n_se={n_se}, turnaround SE={turnaround}, CTD depth={-se.z[turnaround]:.0f}m")

def run_inv(se_mod):
    return compute_inverse(se_mod, u_ship=0.0017, v_ship=-0.0019,
        sadcp_z=npz["z"], sadcp_u=npz["u"], sadcp_v=npz["v"], sadcp_err=npz["err"])

result_full = run_inv(se)

# Downcast-only (blank upcast SE weights)
w_down = se.weight.copy(); w_down[:, turnaround:] = np.nan
se_down = dataclasses.replace(se, weight=w_down)
result_down = run_inv(se_down)

def rmse(result):
    valid = np.isfinite(ref_u) & (ref_nvel >= 3)
    py_u = np.interp(ref_z[valid], result.z, result.u)
    return float(np.sqrt(np.mean((py_u - ref_u[valid]) ** 2)))

print(f"Full RMSE={rmse(result_full):.4f}, Downcast-only RMSE={rmse(result_down):.4f}")
print(f"\n{'depth':>6}  {'u_full':>8}  {'u_down':>8}  {'ref_u':>8}  {'bias_full':>10}  {'bias_down':>10}")
for z_w in range(0, 3500, 250):
    uf = float(np.interp(z_w, result_full.z, result_full.u))
    ud = float(np.interp(z_w, result_down.z, result_down.u))
    ri = np.argmin(np.abs(ref_z - z_w)); ru = ref_u[ri]
    if np.isfinite(ru):
        print(f"  {z_w:4d}m  {uf:+8.4f}  {ud:+8.4f}  {ru:+8.4f}  {uf-ru:+10.4f}  {ud-ru:+10.4f}")

print()
for band_lo, band_hi in [(0,500),(500,1000),(1000,2000),(2000,3000),(3000,4400)]:
    mask = (ref_z >= band_lo) & (ref_z < band_hi) & np.isfinite(ref_u) & (ref_nvel >= 3)
    if mask.sum() < 3: continue
    pf = np.interp(ref_z[mask], result_full.z, result_full.u)
    pd = np.interp(ref_z[mask], result_down.z, result_down.u)
    valid = np.isfinite(pf) & np.isfinite(pd)
    if valid.sum() < 3: continue
    rf = float(np.sqrt(np.mean((pf[valid] - ref_u[mask][valid]) ** 2)))
    rd = float(np.sqrt(np.mean((pd[valid] - ref_u[mask][valid]) ** 2)))
    cf = float(np.corrcoef(ref_u[mask][valid], pf[valid])[0, 1])
    cd = float(np.corrcoef(ref_u[mask][valid], pd[valid])[0, 1])
    print(f"{band_lo}-{band_hi}m: full rmse={rf:.4f} corr={cf:+.3f} | down rmse={rd:.4f} corr={cd:+.3f}")
