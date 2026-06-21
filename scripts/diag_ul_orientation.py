"""Diagnostic: sweep UL orientation variants to find the best beam2earth convention.

Runs the full inverse solver for 5 variants and prints per-band correlation and RMSE
against the 003.nc LDEO_IX reference. Mirrors the method used to pin down magnetic
declination sign.

Usage:
    TEST_DATA_DIR=test_data/2015_P16N python scripts/diag_ul_orientation.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import netCDF4
import numpy as np

# Add src to path so we can import ladcp without installing
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ladcp.ingestion.ctd import assign_bin_depths, compute_ship_velocity, load_ctd
from ladcp.ingestion.rdi import load_rdi
from ladcp.qa.editing import edit_sidelobes
from ladcp.solution.inverse import EnsembleData, compute_inverse, prepare_superensembles
from ladcp.transforms.beam2earth import beam2earth, uvrot

THETA_DEG = 20.0
DROT_DEG = 12.318441


def load_dl(dl_path: Path, cnv_path: Path, sadcp_path: Path):
    """Load DL data, apply beam2earth + declination, assign depths. Returns dict."""
    rdi = load_rdi(dl_path)
    ctd = load_ctd(cnv_path)

    u_earth, v_earth, w_earth = beam2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG, gimbaled=True,
    )
    u_earth, v_earth = uvrot(u_earth, v_earth, -DROT_DEG)

    z_m, izm_pos = assign_bin_depths(rdi, ctd, looker="down")

    weight = np.nanmean(rdi.corr.astype(np.float64), axis=2) / 128.0

    bt_u_e, bt_v_e, bt_w_e = beam2earth(
        rdi.btrack_vel_ms[0], rdi.btrack_vel_ms[1],
        rdi.btrack_vel_ms[2], rdi.btrack_vel_ms[3],
        rdi.heading, rdi.pitch, rdi.roll,
        THETA_DEG, gimbaled=True,
    )
    bt_u_e, bt_v_e = uvrot(bt_u_e, bt_v_e, -DROT_DEG)
    bvel = np.stack([bt_u_e, bt_v_e, bt_w_e], axis=1)
    bvels = np.full_like(bvel, 0.02)
    hbot = np.nanmean(rdi.btrack_range_m, axis=0)

    if ctd.lat is not None:
        slat = np.interp(rdi.time_julian, ctd.time_julian, ctd.lat, left=np.nan, right=np.nan)
        slon = np.interp(rdi.time_julian, ctd.time_julian, ctd.lon, left=np.nan, right=np.nan)
        u_ship, v_ship = compute_ship_velocity(ctd.lat, ctd.lon, ctd.time_julian)
    else:
        slat = np.full(rdi.nens, np.nan)
        slon = np.full(rdi.nens, np.nan)
        u_ship, v_ship = 0.0, 0.0

    sadcp_z = sadcp_u = sadcp_v = sadcp_err = None
    if sadcp_path.exists():
        npz = np.load(sadcp_path)
        sadcp_z, sadcp_u, sadcp_v, sadcp_err = npz["z"], npz["u"], npz["v"], npz["err"]

    return dict(
        rdi=rdi, ctd=ctd,
        u=u_earth, v=v_earth, w=w_earth, weight=weight,
        z_neg=-z_m, izm_neg=-izm_pos,
        bvel=bvel, bvels=bvels, hbot=hbot,
        slat=slat, slon=slon,
        u_ship=u_ship, v_ship=v_ship,
        sadcp_z=sadcp_z, sadcp_u=sadcp_u, sadcp_v=sadcp_v, sadcp_err=sadcp_err,
    )


def score_vs_ref(result, ref_path: Path) -> dict:
    """Compare u, v against 003.nc by depth band. Returns RMSE and corr per band."""
    ds = netCDF4.Dataset(ref_path)
    ref_z = np.array(ds.variables["z"][:])
    ref_u = np.array(ds.variables["u"][:])
    ref_v = np.array(ds.variables["v"][:])
    ref_nvel = np.array(ds.variables["nvel"][:])
    ds.close()

    res_u = np.interp(ref_z, result.z, result.u, left=np.nan, right=np.nan)
    res_v = np.interp(ref_z, result.z, result.v, left=np.nan, right=np.nan)

    bands = [(0, 800), (800, 4500), (0, 4500)]
    scores = {}
    for lo, hi in bands:
        mask = (ref_z >= lo) & (ref_z <= hi) & np.isfinite(ref_u) & np.isfinite(res_u) & (ref_nvel >= 3)
        if mask.sum() < 5:
            scores[f"{lo}-{hi}"] = dict(n=0, u_rmse=np.nan, v_rmse=np.nan, u_corr=np.nan, v_corr=np.nan)
            continue
        du = result.u[np.searchsorted(result.z, ref_z[mask])] if False else res_u[mask]
        dv = res_v[mask]
        ru, rv = ref_u[mask], ref_v[mask]
        u_rmse = float(np.sqrt(np.mean((du - ru) ** 2)))
        v_rmse = float(np.sqrt(np.mean((dv - rv) ** 2)))
        u_corr = float(np.corrcoef(du, ru)[0, 1]) if len(du) > 1 else np.nan
        v_corr = float(np.corrcoef(dv, rv)[0, 1]) if len(dv) > 1 else np.nan
        scores[f"{lo}-{hi}"] = dict(n=int(mask.sum()), u_rmse=u_rmse, v_rmse=v_rmse,
                                    u_corr=u_corr, v_corr=v_corr)
    return scores


def run_variant(dl: dict, ul_data: dict | None, label: str):
    """Build EnsembleData from DL + optional UL, run inverse, return result."""
    rdi = dl["rdi"]
    n_dl = rdi.nbin

    if ul_data is None:
        # DL-only baseline
        ens = EnsembleData(
            u=dl["u"], v=dl["v"], w=dl["w"], weight=dl["weight"],
            izm=dl["izm_neg"], z=dl["z_neg"],
            time_jul=rdi.time_julian,
            bvel=dl["bvel"], bvels=dl["bvels"], hbot=dl["hbot"],
            izd=np.arange(n_dl, dtype=int),
            izu=np.array([], dtype=int),
            slat=dl["slat"], slon=dl["slon"],
        )
    else:
        rdi_ul = ul_data["rdi_ul"]
        n_ul = rdi_ul.nbin

        # Time-align UL to DL: for each DL ensemble, find nearest UL ensemble
        ul_indices = np.argmin(
            np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0
        )  # shape (n_dl_ens,)

        u_ul_a = ul_data["u_ul"][:, ul_indices]  # (n_ul_bins, n_dl_ens)
        v_ul_a = ul_data["v_ul"][:, ul_indices]
        w_ul_a = ul_data["w_ul"][:, ul_indices]
        w_ul_a_q = ul_data["weight_ul"][:, ul_indices]
        izm_ul_a = ul_data["izm_ul"][:, ul_indices]

        # Merge: reversed UL bins on top of DL bins (shallow→deep)
        u_comb = np.vstack([u_ul_a[::-1, :], dl["u"]])
        v_comb = np.vstack([v_ul_a[::-1, :], dl["v"]])
        w_comb = np.vstack([w_ul_a[::-1, :], dl["w"]])
        wt_comb = np.vstack([w_ul_a_q[::-1, :], dl["weight"]])
        izm_comb = np.vstack([izm_ul_a[::-1, :], dl["izm_neg"]])

        # izu: each element is the combined-array row index for that UL bin
        #   UL bin 0 (deepest) → combined row n_ul-1
        #   UL bin n_ul-1 (shallowest) → combined row 0
        izu = np.arange(n_ul - 1, -1, -1, dtype=int)
        izd = np.arange(n_ul, n_ul + n_dl, dtype=int)

        ens = EnsembleData(
            u=u_comb, v=v_comb, w=w_comb, weight=wt_comb,
            izm=izm_comb, z=dl["z_neg"],
            time_jul=rdi.time_julian,
            bvel=dl["bvel"], bvels=dl["bvels"], hbot=dl["hbot"],
            izd=izd, izu=izu,
            slat=dl["slat"], slon=dl["slon"],
        )

    ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)
    se = prepare_superensembles(ens, dz=16.0)
    result = compute_inverse(
        se,
        u_ship=dl["u_ship"], v_ship=dl["v_ship"],
        sadcp_z=dl["sadcp_z"], sadcp_u=dl["sadcp_u"],
        sadcp_v=dl["sadcp_v"], sadcp_err=dl["sadcp_err"],
    )
    return result


def make_ul_variant(rdi_ul, ctd, variant: str) -> dict:
    """Apply beam2earth to UL with given orientation variant."""
    h = rdi_ul.heading.copy()
    p = rdi_ul.pitch.copy()
    r = rdi_ul.roll.copy()

    if variant == "heading_180":
        h = (h + 180.0) % 360.0
    elif variant == "pitch_flip":
        p = -p
    elif variant == "roll_flip":
        r = -r
    elif variant == "pitch_roll_flip":
        p = -p
        r = -r
    # else: "no_flip" — use as-is

    u_ul, v_ul, w_ul = beam2earth(
        rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
        h, p, r,
        THETA_DEG, gimbaled=True,
    )
    u_ul, v_ul = uvrot(u_ul, v_ul, -DROT_DEG)

    _, izm_pos = assign_bin_depths(rdi_ul, ctd, looker="up")
    weight_ul = np.nanmean(rdi_ul.corr.astype(np.float64), axis=2) / 128.0

    return dict(
        rdi_ul=rdi_ul,
        u_ul=u_ul, v_ul=v_ul, w_ul=w_ul,
        weight_ul=weight_ul,
        izm_ul=-izm_pos,
    )


def main():
    env = os.environ.get("TEST_DATA_DIR", "test_data/2015_P16N")
    base = Path(env)

    dl_path = base / "003DL000.000"
    ul_path = base / "003UL000.000"
    cnv_path = base / "003_01.cnv"
    ref_path = base / "003.nc"
    sadcp_path = base / "sadcp_003.npz"

    for p in [dl_path, ul_path, cnv_path, ref_path]:
        if not p.exists():
            print(f"MISSING: {p}")
            sys.exit(1)

    print("Loading DL + CTD...")
    dl = load_dl(dl_path, cnv_path, sadcp_path)
    rdi_ul = load_rdi(ul_path)
    ctd = dl["ctd"]

    variants = [
        ("DL only (baseline)", None),
        ("DL + UL, no_flip", "no_flip"),
        ("DL + UL, pitch_flip", "pitch_flip"),
        ("DL + UL, roll_flip", "roll_flip"),
        ("DL + UL, pitch_roll_flip", "pitch_roll_flip"),
        ("DL + UL, heading_180", "heading_180"),
    ]

    print(f"\n{'Variant':<35} {'Band':>12} {'n':>5} {'u_rmse':>8} {'v_rmse':>8} {'u_corr':>8} {'v_corr':>8}")
    print("-" * 95)

    for label, vname in variants:
        if vname is None:
            ul_data = None
        else:
            ul_data = make_ul_variant(rdi_ul, ctd, vname)

        result = run_variant(dl, ul_data, label)
        scores = score_vs_ref(result, ref_path)

        first = True
        for band, s in scores.items():
            row_label = label if first else ""
            first = False
            print(f"{row_label:<35} {band:>12} {s['n']:>5} "
                  f"{s['u_rmse']:>8.4f} {s['v_rmse']:>8.4f} "
                  f"{s['u_corr']:>8.4f} {s['v_corr']:>8.4f}")
        print()


if __name__ == "__main__":
    main()
