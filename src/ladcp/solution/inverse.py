"""LADCP inverse velocity solver.

Replicates prepinv.m (super-ensemble formation) and getinv.m (constrained
least-squares inversion) from the LDEO_IX MATLAB reference implementation.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass


@dataclass
class EnsembleData:
    """Earth-frame ADCP data aligned with CTD depths, input to prepare_superensembles().

    Depth convention: negative = below surface. izd / izu are 0-indexed row indices
    into the first dimension of u/v/w/weight/izm.
    """
    u: np.ndarray         # (n_bins, n_ens) eastward velocity m/s
    v: np.ndarray         # (n_bins, n_ens) northward velocity m/s
    w: np.ndarray         # (n_bins, n_ens) vertical velocity m/s
    weight: np.ndarray    # (n_bins, n_ens) quality weight 0–1
    izm: np.ndarray       # (n_bins, n_ens) bin depth m (≤ 0)
    z: np.ndarray         # (n_ens,) CTD depth m (≤ 0)
    time_jul: np.ndarray  # (n_ens,) Julian day
    bvel: np.ndarray      # (n_ens, 3) bottom track u, v, w m/s (NaN = absent)
    bvels: np.ndarray     # (n_ens, 3) bottom track std m/s
    hbot: np.ndarray      # (n_ens,) height above bottom m (NaN = absent)
    izd: np.ndarray       # (n_dl_bins,) downlooker bin row indices (int)
    izu: np.ndarray       # (n_ul_bins,) uplooker bin row indices (int)
    slat: np.ndarray      # (n_ens,) latitude (NaN if unavailable)
    slon: np.ndarray      # (n_ens,) longitude (NaN if unavailable)


@dataclass
class SuperEnsemble:
    """Depth-window-averaged ADCP data, output of prepare_superensembles().

    Shape conventions: 2-D fields are (n_bins, n_se); 1-D fields are (n_se,).
    bvel / bvels are (3, n_se) — transposed from EnsembleData for column access.
    """
    ru: np.ndarray        # (n_bins, n_se) eastward velocity m/s
    rv: np.ndarray        # (n_bins, n_se) northward velocity m/s
    rw: np.ndarray        # (n_bins, n_se) vertical velocity m/s
    ruvs: np.ndarray      # (n_bins, n_se) combined U+V std (velocity uncertainty)
    weight: np.ndarray    # (n_bins, n_se) quality weight 0–1
    izm: np.ndarray       # (n_bins, n_se) bin depth m (≤ 0)
    z: np.ndarray         # (n_se,) CTD depth m (≤ 0)
    dt: np.ndarray        # (n_se,) time interval s
    time_jul: np.ndarray  # (n_se,) Julian day
    bvel: np.ndarray      # (3, n_se) bottom track u, v, w m/s
    bvels: np.ndarray     # (3, n_se) bottom track std m/s
    hbot: np.ndarray      # (n_se,) height above bottom m
    slat: np.ndarray      # (n_se,) latitude
    slon: np.ndarray      # (n_se,) longitude
    izd: np.ndarray       # downlooker bin row indices (unchanged from input)
    izu: np.ndarray       # uplooker bin row indices (unchanged from input)


def _window_boundaries(depth0: np.ndarray, avdz: float) -> list[np.ndarray]:
    """Partition ensemble indices into depth-triggered windows.

    Replicates the while-loop in prepinv.m lines 499–605. Each window spans
    consecutive ensembles until |depth0[t] - depth0[window_start]| > avdz.
    """
    n = len(depth0)
    windows: list[np.ndarray] = []
    ilast = 0
    while ilast < n:
        if ilast + 1 >= n:
            windows.append(np.array([ilast]))
            break
        depth_change = np.abs(depth0[ilast + 1:] - depth0[ilast])
        ii = np.where(depth_change > avdz)[0]
        end = int(ii[0]) + 1 if len(ii) > 0 else n - ilast
        i1 = np.arange(ilast, min(ilast + end, n))
        windows.append(i1)
        ilast = int(i1[-1]) + 1
    return windows


def prepare_superensembles(
    ens: EnsembleData,
    *,
    dz: float | None = None,
) -> SuperEnsemble:
    """Form super-ensembles by depth-window averaging (replicates prepinv.m).

    Parameters
    ----------
    ens:
        Earth-frame ADCP data from beam2earth + assign_bin_depths.
    dz:
        Depth window size m. Defaults to ``median(|diff(izm[:, 0])|)``
        (all bins of the first ensemble), matching MATLAB's
        ``medianan(abs(diff(d.izm(:,1))))``.
    """
    if dz is None:
        dz = float(np.nanmedian(np.abs(np.diff(ens.izm[:, 0]))))

    # Reference bins: 2nd and 3rd downlooker bins (0-indexed: offset 1, 2 from min(izd))
    # Replicates: izr = [min(d.izd)+1, min(d.izd)+2]  (MATLAB 1-indexed → same offsets)
    izr: np.ndarray
    if len(ens.izd) > 2:
        base = int(ens.izd.min())
        izr = np.array([base + 1, base + 2], dtype=int)
    else:
        izr = ens.izd.copy()
    if len(ens.izu) > 2:
        ul_top = int(ens.izu.max())
        izr = np.concatenate([izr, [ul_top - 1, ul_top - 2]])

    windows = _window_boundaries(ens.z, dz)
    n_se = len(windows)
    n_bins = ens.u.shape[0]

    ru = np.full((n_bins, n_se), np.nan)
    rv = np.full((n_bins, n_se), np.nan)
    rw = np.full((n_bins, n_se), np.nan)
    ruvs = np.full((n_bins, n_se), np.nan)
    weight_se = np.full((n_bins, n_se), np.nan)
    izm_se = np.full((n_bins, n_se), np.nan)
    z_se = np.full(n_se, np.nan)
    time_se = np.full(n_se, np.nan)
    bvel_se = np.full((3, n_se), np.nan)
    bvels_se = np.full((3, n_se), np.nan)
    hbot_se = np.full(n_se, np.nan)
    slat_se = np.full(n_se, np.nan)
    slon_se = np.full(n_se, np.nan)

    for im, i1 in enumerate(windows):
        u_win = ens.u[:, i1]   # (n_bins, n_win) — no weight masking here
        v_win = ens.v[:, i1]   # (MATLAB uses w=d.weight*0+1 = all ones)
        w_win = ens.w[:, i1]

        # Per-ensemble reference velocity: median of reference bins
        # MATLAB: ur = medianan(d.ru(izr, i1)) — plain median, no weight applied
        izr_valid = izr[izr < n_bins]
        ur_t = np.nanmedian(u_win[izr_valid], axis=0)  # (n_win,)
        vr_t = np.nanmedian(v_win[izr_valid], axis=0)

        # Remove per-ensemble reference, take median, add back mean reference
        # MATLAB: di.ru(:,im) = medianan(d.ru - ur_broadcast)' + mean(ur)
        u_deref = u_win - ur_t[np.newaxis, :]
        v_deref = v_win - vr_t[np.newaxis, :]

        ru[:, im] = np.nanmedian(u_deref, axis=1) + np.nanmean(ur_t)
        rv[:, im] = np.nanmedian(v_deref, axis=1) + np.nanmean(vr_t)
        rw[:, im] = np.nanmedian(w_win, axis=1)

        # Velocity uncertainty: combined U+V std over window
        # Fix I-4: single-sample window — use ddof=0 to return 0 not NaN (matches stdnan())
        n_win = len(i1)
        ddof = 1 if n_win > 1 else 0
        ruvs[:, im] = np.sqrt(
            np.nanstd(u_win, axis=1, ddof=ddof) ** 2
            + np.nanstd(v_win, axis=1, ddof=ddof) ** 2
        )

        weight_se[:, im] = np.nanmean(ens.weight[:, i1], axis=1)
        izm_se[:, im] = np.nanmean(ens.izm[:, i1], axis=1)
        z_se[im] = np.nanmean(ens.z[i1])
        time_se[im] = np.nanmean(ens.time_jul[i1])

        bvel_se[:, im] = np.nanmean(ens.bvel[i1], axis=0)
        # Fix I-3: detrend w-component by reference vertical velocity before computing std
        # Replicates prepinv.m lines 565-568: bvel(:,3) = bvel(:,3) - wr(1,:)'
        bvel_win = ens.bvel[i1].copy()
        wr_ref = float(np.nanmean(rw[:, im]))
        bvel_win[:, 2] = bvel_win[:, 2] - wr_ref
        bvels_se[:, im] = np.nanstd(bvel_win, axis=0, ddof=1)
        hbot_se[im] = np.nanmean(ens.hbot[i1])
        slat_se[im] = np.nanmedian(ens.slat[i1])
        slon_se[im] = np.nanmedian(ens.slon[i1])

    # Time interval between super-ensembles (seconds); mirror edge values
    # Fix M-1: guard against n_se==1 where dt_mid is empty
    dt_mid = np.diff(time_se) * 24.0 * 3600.0
    if len(dt_mid) == 0:
        dt = np.zeros(1)
    else:
        dt = np.concatenate([[dt_mid[0]], (dt_mid[:-1] + dt_mid[1:]) / 2.0, [dt_mid[-1]]])

    # Floor std at single_ping_err (≈0.01 m/s); propagate NaN from weight
    # Replicates prepinv.m's superens_std_min logic
    single_ping_err = 0.01
    zero_mask = ruvs == 0
    weight_se[zero_mask] = np.nan
    ruvs[ruvs < single_ping_err] = single_ping_err
    ruvs = ruvs + weight_se * 0  # NaN-propagate

    return SuperEnsemble(
        ru=ru, rv=rv, rw=rw, ruvs=ruvs, weight=weight_se, izm=izm_se,
        z=z_se, dt=dt, time_jul=time_se,
        bvel=bvel_se, bvels=bvels_se, hbot=hbot_se,
        slat=slat_se, slon=slon_se,
        izd=ens.izd, izu=ens.izu,
    )
