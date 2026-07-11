"""End-to-end LADCP processing pipeline (LDEO_IX-equivalent).

process_cast() chains the validated stages -- ingestion, beam->earth with
3-beam solutions, CTD depth registration (Saunders + clock-lag), LDEO
editing and weights, sound-speed correction, UL/DL merge with w-lag
pairing, rotup2down/offsetup2down, super-ensemble formation, inverse
solution -- exactly as validated against LDEO_IX on GO-SHIP P16N 2015
cast 003 (u RMSE 0.045 / v 0.033 vs the archived LDEO answer; Stage A
and formation match the MATLAB code under Octave to machine precision;
see octave_harness/REPORT.md).

Every cast-specific quantity is a parameter; CastParams.from_ldeo_nc()
reads them from an LDEO_IX output file's global attributes when
validating against archived results.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ladcp.ingestion.ctd import (
    CTDTimeSeries,
    assign_bin_depths,
    ctd_in_water_window,
    estimate_ctd_adcp_lag,
    load_ctd,
)
from ladcp.ingestion.rdi import best_ul_shift, cut_ensembles, load_rdi
from ladcp.qa.editing import (
    build_ldeo_weights,
    edit_error_velocity,
    edit_large_velocities,
    edit_mask_bins,
    edit_outliers,
    edit_ppi,
    edit_sidelobes,
    edit_w_outliers,
    tilt_from_pitch_roll,
)
from ladcp.solution.inverse import (
    EnsembleData,
    InverseResult,
    compute_inverse,
    offsetup2down,
    prepare_superensembles,
    rotup2down,
)
from ladcp.transforms.beam2earth import beam2earth, janus_error_velocity, uvrot
from ladcp.transforms.soundspeed import (
    apply_sound_speed_correction,
    depth_to_pressure,
    sound_speed,
)


@dataclass
class CastParams:
    """Cast- and instrument-specific processing parameters."""
    lat_deg: float                      # for Saunders pressure->depth
    drot_deg: float                     # magnetic declination East
    superens_std_min: float = 0.1      # Single_Ping_Err / sqrt(npng)
    theta_deg: float = 20.0            # beam angle (WH300: 20)
    u_ship: float | None = None        # GPS-derived mean ship velocity
    v_ship: float | None = None
    sadcp_z: np.ndarray | None = None  # optional SADCP constraint
    sadcp_u: np.ndarray | None = None
    sadcp_v: np.ndarray | None = None
    sadcp_err: np.ndarray | None = None
    edit_ppi: bool = False             # previous-ping-interference edit
    edit_ppi_layer_thickness: float = 180.0   # (edit_data.m defaults)
    edit_ppi_max_hab: float = 1000.0

    @classmethod
    def from_ldeo_nc(cls, path: Path | str, **overrides) -> CastParams:
        """Read cast parameters from an LDEO_IX output NetCDF's attributes.

        Uses GEN_Magnetic_deviation_deg (or drot), lat, uship/vship, and
        LADCP_dn_conf_single_ping_acc / _number_pings (superens_std_min =
        single-ping accuracy / sqrt(pings), prepinv.m line 41). When the
        archived run used a SADCP constraint, ladcp2cdf embeds the exact
        profile it consumed (z_sadcp/u_sadcp/v_sadcp/uerr_sadcp variables);
        those are read back so a validation rerun applies the same
        constraint. Without it, the inverse has large-vertical-scale
        freedom the archive did not (measured on I7N 003: a u tilt of
        +0.21 m/s above 1000 m vs -0.11 below 3000 m with the depth mean
        pinned by uship/vship).
        """
        import netCDF4

        ds = netCDF4.Dataset(str(path))
        try:
            drot_attr = getattr(ds, "GEN_Magnetic_deviation_deg", None)
            if drot_attr is None:
                drot_attr = ds.drot
            drot = float(drot_attr)
            lat = float(np.asarray(ds.variables["lat"][:]).ravel()[0]) \
                if "lat" in ds.variables else float(ds.lat)
            spe = getattr(ds, "LADCP_dn_conf_single_ping_acc", None)
            npng = getattr(ds, "LADCP_dn_conf_number_pings", 1)
            std_min = (float(spe) / math.sqrt(max(float(npng), 1.0))
                       if spe is not None and np.isfinite(float(spe))
                       else 0.1)
            kw = dict(
                lat_deg=lat, drot_deg=drot, superens_std_min=std_min,
                u_ship=float(ds.uship) if hasattr(ds, "uship") else None,
                v_ship=float(ds.vship) if hasattr(ds, "vship") else None,
            )
            # Per-cast PPI editing, recorded verbatim in some LDEO_IX
            # outputs (A16N 2013 IX_8: edit_PPI/edit_PPI_layer_thickness/
            # edit_PPI_max_hab global attributes).
            ppi = getattr(ds, "edit_PPI", None)
            if ppi is not None:
                kw["edit_ppi"] = bool(float(ppi))
                thick = getattr(ds, "edit_PPI_layer_thickness", None)
                if thick is not None:
                    kw["edit_ppi_layer_thickness"] = float(thick)
                hab = getattr(ds, "edit_PPI_max_hab", None)
                if hab is not None:
                    kw["edit_ppi_max_hab"] = float(hab)
            if "u_sadcp" in ds.variables:
                def _var(name: str) -> np.ndarray:
                    return np.ma.filled(
                        ds.variables[name][:].astype(np.float64), np.nan)
                sz, su, sv = _var("z_sadcp"), _var("u_sadcp"), _var("v_sadcp")
                serr = (_var("uerr_sadcp") if "uerr_sadcp" in ds.variables
                        else np.full_like(sz, 0.05))
                ok = np.isfinite(sz) & np.isfinite(su) & np.isfinite(sv)
                if ok.sum() >= 3:
                    kw.update(
                        sadcp_z=sz[ok], sadcp_u=su[ok], sadcp_v=sv[ok],
                        sadcp_err=np.where(np.isfinite(serr[ok]),
                                           serr[ok], 0.05),
                    )
        finally:
            ds.close()
        kw.update(overrides)
        return cls(**kw)


def process_cast(
    dl_path: Path | str,
    ul_path: Path | str,
    ctd: Path | str | CTDTimeSeries,
    params: CastParams,
    *,
    rot: bool = False,
    offset: bool = False,
    stages: dict | None = None,
) -> InverseResult:
    """Run the full validated pipeline on one cast.

    rot/offset enable prepinv.m's rotup2down/offsetup2down passes. Both
    default ON in LDEO_IX, but OFF here: measured on P16N 003 (2026-07-11,
    post formation parity), rotup2down is neutral vs the archived answer
    (u 0.0454 vs 0.0450) and offsetup2down WORSENS it (u 0.0611) --
    plausibly because LDEO's step-11 lanarrow trim between the two solves
    is not ported. The <0.05 validation targets are met with both off.
    stages, if a dict, receives intermediate outputs (post_transform,
    post_edit, superensembles, result).
    """
    theta = params.theta_deg
    rdi = load_rdi(Path(dl_path))
    rdi_ul = load_rdi(Path(ul_path))
    if not isinstance(ctd, CTDTimeSeries):
        ctd = load_ctd(Path(ctd))

    dl_kw = dict(gimbaled=False, beams_up=False, allow_3beam=True)
    ul_kw = dict(gimbaled=False, beams_up=True, allow_3beam=True)

    u_dl, v_dl, w_dl = beam2earth(
        rdi.u, rdi.v, rdi.w, rdi.e,
        rdi.heading, rdi.pitch, rdi.roll, theta, **dl_kw,
    )
    u_dl, v_dl, w_dl = edit_error_velocity(
        u_dl, v_dl, w_dl,
        janus_error_velocity(rdi.u, rdi.v, rdi.w, rdi.e, theta),
    )
    u_dl, v_dl = uvrot(u_dl, v_dl, -params.drot_deg)

    # CTD-ADCP clock offset (loadctd.m besttlag equivalent).
    _, lagdt_days, _ = estimate_ctd_adcp_lag(
        rdi.time_julian, np.nanmedian(w_dl, axis=0), ctd,
        lat_deg=params.lat_deg,
    )
    # Discard ensembles outside the CTD in-water window (loadctd.m:517),
    # e.g. pre-deployment / post-recovery pinging on deck.
    t_in0, t_in1 = ctd_in_water_window(ctd)
    keep = ((rdi.time_julian + lagdt_days >= t_in0)
            & (rdi.time_julian + lagdt_days <= t_in1))
    if not keep.all():
        rdi = cut_ensembles(rdi, keep)
        u_dl, v_dl, w_dl = u_dl[:, keep], v_dl[:, keep], w_dl[:, keep]
    z_m, izm_dl_pos = assign_bin_depths(
        rdi, ctd, looker="down", lat_deg=params.lat_deg,
        time_offset_days=lagdt_days,
    )
    cm_dl = np.median(rdi.corr.astype(np.float64), axis=2)
    ts_dl = np.median(rdi.echo.astype(np.float64), axis=2)

    u_ul, v_ul, w_ul = beam2earth(
        rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
        rdi_ul.heading, rdi_ul.pitch, rdi_ul.roll, theta, **ul_kw,
    )
    u_ul, v_ul, w_ul = edit_error_velocity(
        u_ul, v_ul, w_ul,
        janus_error_velocity(rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e, theta),
    )
    u_ul, v_ul = uvrot(u_ul, v_ul, -params.drot_deg)
    _, izm_ul_pos = assign_bin_depths(
        rdi_ul, ctd, looker="up", lat_deg=params.lat_deg,
        time_offset_days=lagdt_days,
    )
    cm_ul = np.median(rdi_ul.corr.astype(np.float64), axis=2)
    ts_ul = np.median(rdi_ul.echo.astype(np.float64), axis=2)

    # UL->DL pairing: nearest time refined by w cross-correlation
    # (loadrdi.m merge bestlag; SEQUENCE shift).
    ul_idx = np.argmin(
        np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0
    )
    ul_shift, _ = best_ul_shift(w_dl, w_ul, ul_idx)
    ul_idx = ul_idx[np.clip(np.arange(len(ul_idx)) + ul_shift,
                            0, len(ul_idx) - 1)]

    n_ul, n_dl = rdi_ul.nbin, rdi.nbin
    u_comb = np.vstack([u_ul[:, ul_idx][::-1, :], u_dl])
    v_comb = np.vstack([v_ul[:, ul_idx][::-1, :], v_dl])
    w_comb = np.vstack([w_ul[:, ul_idx][::-1, :], w_dl])
    izm_comb = np.vstack([-izm_ul_pos[:, ul_idx][::-1, :], -izm_dl_pos])
    izu = np.arange(n_ul - 1, -1, -1, dtype=int)
    izd = np.arange(n_ul, n_ul + n_dl, dtype=int)
    cm_comb = np.vstack([cm_ul[:, ul_idx][::-1, :], cm_dl])
    ts_comb = np.vstack([ts_ul[:, ul_idx][::-1, :], ts_dl])
    weight_comb = build_ldeo_weights(
        cm_comb, ts_comb, rdi.pitch, rdi.roll, v_comb, w_comb, izd, izu,
    )

    bt_u, bt_v, bt_w = beam2earth(
        rdi.btrack_vel_ms[0], rdi.btrack_vel_ms[1],
        rdi.btrack_vel_ms[2], rdi.btrack_vel_ms[3],
        rdi.heading, rdi.pitch, rdi.roll, theta, **dl_kw,
    )
    bt_u, bt_v, bt_w = edit_error_velocity(
        bt_u, bt_v, bt_w,
        janus_error_velocity(rdi.btrack_vel_ms[0], rdi.btrack_vel_ms[1],
                             rdi.btrack_vel_ms[2], rdi.btrack_vel_ms[3],
                             theta),
    )
    bt_u, bt_v = uvrot(bt_u, bt_v, -params.drot_deg)
    bvel = np.stack([bt_u, bt_v, bt_w], axis=1)

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

    ens = edit_outliers(ens)
    temp_i = np.interp(
        rdi.time_julian + lagdt_days, ctd.time_julian, ctd.temp_c
    )
    ss = sound_speed(depth_to_pressure(z_m), temp_i, 34.5)
    ens = apply_sound_speed_correction(
        ens, ss=ss, sv_dl=rdi.sound_vel_ms, sv_ul=rdi_ul.sound_vel_ms[ul_idx],
    )
    ens = edit_sidelobes(ens, theta_deg=theta, cell_size_m=rdi.blen_m)
    if params.edit_ppi:
        ens = edit_ppi(
            ens, npng=rdi.npng, beam_angle_deg=theta, ss=ss,
            layer_thickness_m=params.edit_ppi_layer_thickness,
            max_hab_m=params.edit_ppi_max_hab,
        )
    ens = edit_large_velocities(ens)
    ens = edit_w_outliers(ens)
    ens = edit_mask_bins(
        ens,
        dn_bins=[0] if rdi.blnk_m == 0 else [],
        up_bins=[0] if rdi_ul.blnk_m == 0 else [],
    )
    if stages is not None:
        stages["post_edit"] = ens

    if rot or offset:
        ens, _ = rotup2down(ens, rdi.heading, rdi_ul.heading[ul_idx])

    nblock = int(math.ceil(
        5.0 / (float(np.nanmean(np.diff(rdi.time_julian))) * 24.0 * 60.0)))
    tilt = tilt_from_pitch_roll(rdi.pitch, rdi.roll)

    def _solve(e: EnsembleData) -> InverseResult:
        se = prepare_superensembles(
            e, superens_std_min=params.superens_std_min,
            outlier_nblock=nblock, tilt_deg=tilt,
        )
        if stages is not None:
            stages["superensembles"] = se
        return compute_inverse(
            se, u_ship=params.u_ship, v_ship=params.v_ship,
            sadcp_z=params.sadcp_z, sadcp_u=params.sadcp_u,
            sadcp_v=params.sadcp_v, sadcp_err=params.sadcp_err,
        )

    if offset:
        first_guess = _solve(ens)
        ens, _ = offsetup2down(ens, first_guess.z, first_guess.u,
                               first_guess.v)
    res = _solve(ens)
    if stages is not None:
        stages["result"] = res
    return res
