"""Depth-stratified RMSE vs LDEO reference, before/after the beams_up fix.

Runs the exact pipeline of tests/integration/test_inverse_p16n_cast003.py in
two conventions:
  new    - loadrdi convention (gimbaled=False, beams_up per instrument, raw
           UL attitude): the fixed production path.
  legacy - the pre-fix path, reproduced exactly: up-looking beam matrix for
           BOTH instruments (beams_up=True), gimbaled=True, UL pitch negated.

Usage:  TEST_DATA_DIR=test_data uv run python scripts/diag_rmse_strata.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ladcp.ingestion.ctd import (  # noqa: E402
    assign_bin_depths,
    estimate_ctd_adcp_lag,
    load_ctd,
)
from ladcp.ingestion.rdi import best_ul_shift, load_rdi  # noqa: E402
from ladcp.qa.editing import (  # noqa: E402
    build_ldeo_weights,
    edit_large_velocities,
    edit_mask_bins,
    edit_outliers,
    edit_sidelobes,
    edit_w_outliers,
    tilt_from_pitch_roll,
)
from ladcp.solution.inverse import (  # noqa: E402
    EnsembleData,
    compute_inverse,
    offsetup2down,
    prepare_superensembles,
    rotup2down,
)
from ladcp.transforms.beam2earth import beam2earth, uvrot  # noqa: E402
from ladcp.transforms.soundspeed import (  # noqa: E402
    apply_sound_speed_correction,
    depth_to_pressure,
    sound_speed,
)

THETA_DEG = 20.0
DROT_DEG = 12.318441
STRATA = [(0, 1000), (1000, 2000), (2000, 3000), (3000, 4500)]


def run_pipeline(data_dir: Path, legacy: bool, rot: bool = False, offset: bool = False, stages: dict | None = None):
    """stages, if given a dict, is populated with intermediate pipeline
    outputs (post-transform combined arrays, post-edit ensemble, super-
    ensembles, final result) for octave_harness/diff_stages.py (M3)."""
    rdi = load_rdi(data_dir / "003DL000.000")
    rdi_ul = load_rdi(data_dir / "003UL000.000")
    ctd = load_ctd(data_dir / "003_01.cnv")

    ds = netCDF4.Dataset(data_dir / "003.nc")
    lat_deg = float(ds.variables["lat"][:])
    u_ship, v_ship = float(ds.uship), float(ds.vship)
    ref_z = np.array(ds.variables["z"][:])
    ref_u = np.array(ds.variables["u"][:])
    ref_v = np.array(ds.variables["v"][:])
    ref_nvel = np.array(ds.variables["nvel"][:])
    ds.close()

    if legacy:
        dl_kw = dict(gimbaled=True, beams_up=True)  # old up-matrix for DL too
        ul_kw = dict(gimbaled=True, beams_up=True)
        ul_pitch = -rdi_ul.pitch
    else:
        # allow_3beam: loadrdi.m p.allow_3beam_solutions defaults ON (this
        # cast: 14422 DL / 8473 UL 3-beam solutions in LDEO's log).
        dl_kw = dict(gimbaled=False, beams_up=False, allow_3beam=True)
        ul_kw = dict(gimbaled=False, beams_up=True, allow_3beam=True)
        ul_pitch = rdi_ul.pitch

    u_dl, v_dl, w_dl = beam2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll, THETA_DEG, **dl_kw,
    )
    u_dl, v_dl = uvrot(u_dl, v_dl, -DROT_DEG)
    # CTD-ADCP clock offset (loadctd.m besttlag): correlate the CTD's
    # pressure-derived sinking rate against ADCP earth-frame w, then sample
    # pressure at the corrected time and shift the ADCP time labels to the
    # CTD time frame (loadctd.m lines 410-443).
    lag_scans, lagdt_days, lag_corr = estimate_ctd_adcp_lag(
        rdi.time_julian, np.nanmedian(w_dl, axis=0), ctd, lat_deg=lat_deg,
    )
    z_m, izm_dl_pos = assign_bin_depths(
        rdi, ctd, looker="down", lat_deg=lat_deg, time_offset_days=lagdt_days,
    )
    cm_dl = np.median(rdi.corr.astype(np.float64), axis=2)
    ts_dl = np.median(rdi.echo.astype(np.float64), axis=2)

    u_ul, v_ul, w_ul = beam2earth(
        rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
        rdi_ul.heading, ul_pitch, rdi_ul.roll, THETA_DEG, **ul_kw,
    )
    u_ul, v_ul = uvrot(u_ul, v_ul, -DROT_DEG)
    _, izm_ul_pos = assign_bin_depths(
        rdi_ul, ctd, looker="up", lat_deg=lat_deg, time_offset_days=lagdt_days,
    )
    cm_ul = np.median(rdi_ul.corr.astype(np.float64), axis=2)
    ts_ul = np.median(rdi_ul.echo.astype(np.float64), axis=2)

    ul_idx = np.argmin(
        np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0
    )
    # loadrdi.m refines the nearest-time pairing by w cross-correlation
    # ("shift ADCP timeseries by 1 ensembles"); the UL clock is ~0.6 s off
    # the DL clock, so nearest-recorded-time picks the wrong neighbor.
    # SEQUENCE shift (iu = iu(iiu)), not a value shift -- see best_ul_shift.
    ul_shift, _ = best_ul_shift(w_dl, w_ul, ul_idx)
    ul_idx = ul_idx[np.clip(np.arange(len(ul_idx)) + ul_shift, 0, len(ul_idx) - 1)]
    n_ul, n_dl = rdi_ul.nbin, rdi.nbin
    u_comb = np.vstack([u_ul[:, ul_idx][::-1, :], u_dl])
    v_comb = np.vstack([v_ul[:, ul_idx][::-1, :], v_dl])
    w_comb = np.vstack([w_ul[:, ul_idx][::-1, :], w_dl])
    izm_comb = np.vstack([-izm_ul_pos[:, ul_idx][::-1, :], -izm_dl_pos])
    izu = np.arange(n_ul - 1, -1, -1, dtype=int)
    izd = np.arange(n_ul, n_ul + n_dl, dtype=int)
    # loadrdi.m weight: median-over-beams correlation, normalized, with
    # tilt/echo/non-pinging modifiers (build_ldeo_weights docstring).
    cm_comb = np.vstack([cm_ul[:, ul_idx][::-1, :], cm_dl])
    ts_comb = np.vstack([ts_ul[:, ul_idx][::-1, :], ts_dl])
    weight_comb = build_ldeo_weights(
        cm_comb, ts_comb, rdi.pitch, rdi.roll, v_comb, w_comb, izd, izu,
    )

    bt_u, bt_v, bt_w = beam2earth(
        rdi.btrack_vel_ms[0], rdi.btrack_vel_ms[1],
        rdi.btrack_vel_ms[2], rdi.btrack_vel_ms[3],
        rdi.heading, rdi.pitch, rdi.roll, THETA_DEG, **dl_kw,
    )
    bt_u, bt_v = uvrot(bt_u, bt_v, -DROT_DEG)
    bvel = np.stack([bt_u, bt_v, bt_w], axis=1)

    sadcp = data_dir / "sadcp_003.npz"
    npz = np.load(sadcp) if sadcp.exists() else None

    ens = EnsembleData(
        u=u_comb, v=v_comb, w=w_comb, weight=weight_comb,
        izm=izm_comb, z=-z_m, time_jul=rdi.time_julian + lagdt_days,
        bvel=bvel, bvels=np.full_like(bvel, 0.02),
        hbot=np.nanmean(rdi.btrack_range_m, axis=0),
        izd=izd, izu=izu,
        slat=np.full(rdi.nens, np.nan), slon=np.full(rdi.nens, np.nan),
    )
    if stages is not None:
        stages["post_transform"] = ens
    # loadrdi.m runs outlier() at ingestion, before edit_data.m's masking.
    ens = edit_outliers(ens)
    # getdpthi.m sound-speed correction (loadrdi.m:346 hardcodes soundc=0,
    # so LDEO always applies it): true sound speed at the instrument from
    # GEOSECS pressure + CTD temperature at ADCP time, salinity 34.5.
    temp_i = np.interp(
        rdi.time_julian + lagdt_days, ctd.time_julian, ctd.temp_c
    )
    ss = sound_speed(depth_to_pressure(z_m), temp_i, 34.5)
    ens = apply_sound_speed_correction(
        ens, ss=ss, sv_dl=rdi.sound_vel_ms, sv_ul=rdi_ul.sound_vel_ms[ul_idx],
    )
    ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)
    ens = edit_large_velocities(ens)
    ens = edit_w_outliers(ens)
    # edit_data.m: mask bin 1 of any instrument with zero blanking distance.
    ens = edit_mask_bins(
        ens,
        dn_bins=[0] if rdi.blnk_m == 0 else [],
        up_bins=[0] if rdi_ul.blnk_m == 0 else [],
    )
    if stages is not None:
        stages["post_edit"] = ens
    if rot or offset:
        # offsetup2down (process_cast.m step 12) is always paired with
        # rotup2down (step 10) in LDEO's default pipeline.
        ens, _ = rotup2down(ens, rdi.heading, rdi_ul.heading[ul_idx])

    # prepinv.m parity parameters: superens_std_min is the instrument-
    # derived Single_Ping_Err/sqrt(Pings_per_Ensemble) (0.083833 for the
    # P16N WH300s, quoted by the LDEO log); outlier_nblock is LDEO's
    # p.outlier_n, set once at loadrdi from the RAW ping rate.
    import math as _math
    _nblock = int(_math.ceil(
        5.0 / (float(np.nanmean(np.diff(rdi.time_julian))) * 24.0 * 60.0)))
    _tilt = tilt_from_pitch_roll(rdi.pitch, rdi.roll)

    def _solve(ens: EnsembleData):
        se = prepare_superensembles(
            ens, superens_std_min=0.083833,
            outlier_nblock=_nblock, tilt_deg=_tilt,
        )
        if stages is not None:
            stages["superensembles"] = se
        return compute_inverse(
            se, u_ship=u_ship, v_ship=v_ship,
            sadcp_z=npz["z"] if npz is not None else None,
            sadcp_u=npz["u"] if npz is not None else None,
            sadcp_v=npz["v"] if npz is not None else None,
            sadcp_err=npz["err"] if npz is not None else None,
        )

    if offset:
        # process_cast.m steps 10-12: a first solve (rotup2down only) supplies
        # the first-guess profile `dr` that offsetup2down needs, then the
        # ensembles are re-formed with the offset correction applied and
        # solved again. We skip LDEO's step-11 outlier-trimming (lanarrow) —
        # a documented simplification, not a silent shortcut.
        first_guess = _solve(ens)
        ens, _ = offsetup2down(ens, first_guess.z, first_guess.u, first_guess.v)

    res = _solve(ens)
    if stages is not None:
        stages["result"] = res
    return res, (ref_z, ref_u, ref_v, ref_nvel)


def stats(res, ref):
    ref_z, ref_u, ref_v, ref_nvel = ref
    ru = np.interp(ref_z, res.z, res.u, left=np.nan, right=np.nan)
    rv = np.interp(ref_z, res.z, res.v, left=np.nan, right=np.nan)
    valid = np.isfinite(ref_u) & np.isfinite(ru) & (ref_nvel >= 3)
    out = []
    tot_u = float(np.sqrt(np.mean((ru[valid] - ref_u[valid]) ** 2)))
    tot_v = float(np.sqrt(np.mean((rv[valid] - ref_v[valid]) ** 2)))
    out.append(("TOTAL", valid.sum(), tot_u, tot_v,
                float(np.corrcoef(ru[valid], ref_u[valid])[0, 1])))
    for z0, z1 in STRATA:
        s = valid & (ref_z >= z0) & (ref_z < z1)
        if s.sum() > 5:
            out.append((
                f"{z0}-{z1}m", s.sum(),
                float(np.sqrt(np.mean((ru[s] - ref_u[s]) ** 2))),
                float(np.sqrt(np.mean((rv[s] - ref_v[s]) ** 2))),
                float(np.corrcoef(ru[s], ref_u[s])[0, 1]),
            ))
    return out


def main():
    data_dir = Path(os.environ.get("TEST_DATA_DIR", "test_data")) / "2015_P16N"
    configs = [
        ("LEGACY (pre-fix)", True, False, False),
        ("NEW (loadrdi convention)", False, False, False),
        ("NEW + rotup2down only", False, True, False),
        ("NEW + rotup2down + offsetup2down (iterative)", False, True, True),
    ]
    for tag, legacy, rot, offset in configs:
        res, ref = run_pipeline(data_dir, legacy, rot, offset)
        print(f"\n=== {tag} ===")
        hdr = f"  {'stratum':>12s} {'n':>4s} {'u RMSE':>8s} {'v RMSE':>8s} {'r(u)':>6s}"
        print(hdr)
        for name, n, eu, ev, r in stats(res, ref):
            print(f"  {name:>12s} {n:4d} {eu:8.4f} {ev:8.4f} {r:+6.3f}")


if __name__ == "__main__":
    main()
