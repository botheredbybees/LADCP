"""LADCP inverse velocity solver.

Replicates prepinv.m (super-ensemble formation) and getinv.m (constrained
least-squares inversion) from the LDEO_IX MATLAB reference implementation.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg
from dataclasses import dataclass
from scipy.sparse import csr_matrix


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


def _ref_medianan(x: np.ndarray) -> np.ndarray:
    """Replicate MATLAB medianan(x, na=0) column-wise on a (n_bins, n_ens) array.

    For each column (ensemble), sorts valid (finite) rows, then returns the
    element at 1-based index round(n_valid / 2).  With 4 reference bins whose
    DL and UL values form separate clusters, this reliably picks a DL-side
    value rather than the midpoint returned by np.nanmedian, replicating the
    MATLAB prepinv.m reference-velocity extraction behaviour.
    """
    n_rows, n_cols = x.shape
    y = np.full(n_cols, np.nan)
    for j in range(n_cols):
        col = x[:, j]
        valid = col[np.isfinite(col)]
        n = len(valid)
        if n == 0:
            continue
        valid = np.sort(valid)
        idx_1based = round(n / 2)          # MATLAB round(n/2)
        idx = max(0, min(n - 1, idx_1based - 1))  # convert to 0-based
        y[j] = valid[idx]
    return y


def _window_boundaries(depth0: np.ndarray, avdz: float, oversample: float = 1.0) -> list[np.ndarray]:
    """Partition ensemble indices into depth-triggered windows.

    Replicates the while-loop in prepinv.m lines 499–605.

    MATLAB oversample expansion (default 1.0):
      - Initial range i1 spans from ilast+1 to the first ensemble where
        |depth[t] - depth[ilast]| > avdz.
      - Then i1 is expanded symmetrically around its mean to width
        length(i1) * oversample, making adjacent windows overlap.
      - ilast advances to max(i1) of the expanded window.
    Windows may therefore contain overlapping ensembles (smoothing effect).
    """
    n = len(depth0)
    windows: list[np.ndarray] = []
    # MATLAB uses 1-based indexing; ilast=1 means the comparison starts from
    # index 2 in the first iteration, effectively skipping ensemble 0 from
    # inclusion as the comparison anchor only. We use 0-based with the same
    # semantic: ilast is the anchor, i1 starts at ilast+1.
    ilast = 0
    while ilast < n - 1:
        depth_change = np.abs(depth0[ilast + 1:] - depth0[ilast])
        ii = np.where(depth_change > avdz)[0]
        if len(ii) == 0:
            end = n - ilast - 1
        else:
            end = int(ii[0]) + 1
        # Initial range: ilast+1 .. ilast+end (1-indexed: ilast+[1:end])
        i1_init = np.arange(ilast + 1, min(ilast + 1 + end, n))
        if len(i1_init) == 0:
            break
        # Oversample expansion: expand symmetrically around mean of i1_init
        i1l = len(i1_init) / 2.0 * oversample
        center = float(np.mean(i1_init))
        lo = int(round(center - i1l))
        hi = int(round(center + i1l))
        i1 = np.arange(max(0, lo), min(n, hi + 1))
        if len(i1) == 0:
            i1 = i1_init
        # Ensure minimum window size of 2
        if len(i1) < 2:
            i1 = np.array([i1[0], min(i1[0] + 1, n - 1)])
        windows.append(i1)
        ilast = int(i1[-1])
        if ilast >= n - 1:
            break
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
        # MATLAB prepinv.m computes w = d.weight*0+1, which in MATLAB evaluates to 
        # 1 for finite weights and NaN for NaN weights (since NaN * 0 = NaN).
        # Therefore, we MUST mask the velocities here where weight is NaN!
        wt_win = ens.weight[:, i1]
        nan_mask = np.isnan(wt_win)
        u_win = np.where(nan_mask, np.nan, ens.u[:, i1])
        v_win = np.where(nan_mask, np.nan, ens.v[:, i1])
        w_win = np.where(nan_mask, np.nan, ens.w[:, i1])
        izr_valid = izr[izr < n_bins]
        # Replicate MATLAB medianan(x, na=0): pick round(n/2)-th sorted value
        # per ensemble rather than the average of the two middle values.
        # When DL and UL velocities differ (e.g. due to UL compass offset),
        # this consistently selects a DL-cluster value as the reference.
        ur_t = _ref_medianan(u_win[izr_valid])  # (n_win,)
        vr_t = _ref_medianan(v_win[izr_valid])
        wr_t = _ref_medianan(w_win[izr_valid])

        # MATLAB prepinv.m computes the mean of the valid reference velocities,
        # then sets any NaN reference velocities to 0 before dereferencing!
        ruav = np.nanmean(ur_t)
        rvav = np.nanmean(vr_t)
        rwav = np.nanmean(wr_t)
        
        ur_t[np.isnan(ur_t)] = 0.0
        vr_t[np.isnan(vr_t)] = 0.0
        wr_t[np.isnan(wr_t)] = 0.0

        # Remove per-ensemble reference, average, add back mean reference
        u_deref = u_win - ur_t[np.newaxis, :]
        v_deref = v_win - vr_t[np.newaxis, :]
        w_deref = w_win - wr_t[np.newaxis, :]

        def _medianan(x: np.ndarray, n_avg: int) -> np.ndarray:
            # Replicates MATLAB's medianan(x, n) behavior
            x = np.moveaxis(x, 1, -1)
            y = np.full(x.shape[0], np.nan)
            for i in range(x.shape[0]):
                valid = x[i][np.isfinite(x[i])]
                li = len(valid)
                if li > 0:
                    if n_avg > 1:
                        if li > n_avg:
                            valid.sort()
                            # MATLAB: i1 = max([1, round(li/2 - n/2)])
                            i1_matlab = max(1, round(li / 2 - n_avg / 2))
                            i1 = i1_matlab - 1  # 0-indexed
                            y[i] = np.sum(valid[i1 : i1 + n_avg]) / n_avg
                        else:
                            y[i] = np.mean(valid)
                    else:
                        y[i] = valid[0]
            return y

        n_win = len(i1)
        iav = round(n_win / 2)

        ru[:, im] = _medianan(u_deref, iav) + ruav
        rv[:, im] = _medianan(v_deref, iav) + rvav
        rw[:, im] = _medianan(w_deref, iav) + rwav

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
        bvel_win[:, 2] = bvel_win[:, 2] - wr_t
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


def _flatten_obs(
    se: SuperEnsemble,
    velerr: float = 0.05,
    weightmin: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Flatten super-ensemble data into observation vectors for matrix construction.

    Column-major (Fortran) order matches MATLAB reshape(x, nbin*nt, 1).

    Returns
    -------
    d_u, d_v : (n_obs,) observed velocities m/s
    izv      : (n_obs,) positive bin depths m  (= -izm, column-major flattened)
    jprof    : (n_obs,) super-ensemble index 0..n_se-1
    wm       : (n_obs,) data weights (velerr / ruvs, NaN-propagated)
    """
    n_bins, n_se = se.izm.shape

    # Observation depth: positive (izv = -izm), column-major flatten
    izv_full = (-se.izm).ravel(order="F")   # (n_bins * n_se,)

    # Profile index per observation: ensemble 0 repeated n_bins times, etc.
    jprof_full = np.repeat(np.arange(n_se), n_bins)

    d_u_full = se.ru.ravel(order="F")
    d_v_full = se.rv.ravel(order="F")

    # Data weight: velerr / std  (std-based weighting, MATLAB std_weight=1)
    # NaN in weight propagates to wm, excluding those observations
    wm_full = velerr / se.ruvs + se.weight * 0  # NaN-propagate from weight
    wm_full = wm_full.ravel(order="F")

    # Keep only valid, well-weighted observations
    valid = (
        np.isfinite(d_u_full)
        & np.isfinite(d_v_full)
        & np.isfinite(wm_full)
        & (wm_full >= weightmin)
        & (izv_full > 0)
    )

    return (
        d_u_full[valid],
        d_v_full[valid],
        izv_full[valid],
        jprof_full[valid],
        wm_full[valid],
    )


def _build_obs_matrix(izv: np.ndarray, dz: float) -> csr_matrix:
    """Build A_ocean: maps each observation to a depth bin.

    Column j = round(izv[k] / dz) - 1 (0-indexed).
    Replicates lainseta(izv, dz) from getinv.m.
    """
    n_obs = len(izv)
    j = np.round(izv / dz).astype(int) - 1  # 0-indexed depth bin
    j = np.clip(j, 0, None)
    n_zbins = int(j.max()) + 1
    i = np.arange(n_obs)
    return csr_matrix((np.ones(n_obs), (i, j)), shape=(n_obs, n_zbins))


def _build_ctd_matrix(jprof: np.ndarray, n_se: int) -> csr_matrix:
    """Build A_ctd: maps each observation to its super-ensemble (time bin).

    Replicates lainseta(jprof, 1) from getinv.m.
    """
    n_obs = len(jprof)
    i = np.arange(n_obs)
    j = jprof.astype(int)
    return csr_matrix((np.ones(n_obs), (i, j)), shape=(n_obs, n_se))


def _apply_weights(
    A_ocean: csr_matrix,
    A_ctd: csr_matrix,
    d: np.ndarray,
    wm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Apply data weights to observation system; return dense arrays + split indices.

    Replicates lainweig() from getinv.m. Returns dense arrays (not sparse) because
    constraints are appended row-by-row in subsequent tasks.

    Returns
    -------
    A_o, A_c : dense float64 arrays with weights applied
    d_w      : weighted observation vector
    idx_down : row indices belonging to the downcast
    idx_up   : row indices belonging to the upcast
    """
    A_o = A_ocean.toarray() * wm[:, np.newaxis]
    A_c = A_ctd.toarray() * wm[:, np.newaxis]
    d_w = d * wm

    # Split down/up cast at the deepest observation depth.
    # The last column of A_ocean is always the deepest depth bin (n_zbins - 1)
    # because _build_obs_matrix sets n_zbins = j.max() + 1.
    deepest_col = A_ocean.shape[1] - 1
    rows_at_bottom = np.where(A_ocean.getcol(deepest_col).toarray().ravel() > 0)[0]
    if len(rows_at_bottom) > 0:
        split = int(np.median(rows_at_bottom))
    else:
        split = len(d) // 2
    idx_down = np.arange(0, split + 1)
    idx_up = np.arange(split, len(d))

    return A_o, A_c, d_w, idx_down, idx_up


def _add_smoothness(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    smoofac: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Append curvature-penalty rows to the ocean-velocity block (lainsmoo).

    For each interior column j (1..n_cols-2), adds one row with stencil
    [-1, 2, -1] scaled by smoofac * (median_norm / col_norm[j]).
    Also smooths the CTD block symmetrically (MATLAB calls lainsmoo twice).
    Boundary columns get first-derivative (slope) rows: [2,-2] and [-2,2].
    """
    def _smoo_one(A_target: np.ndarray, A_other: np.ndarray,
                  d_in: np.ndarray, fs0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n_rows, n_cols = A_target.shape
        col_norms = np.sqrt(np.abs(A_target).sum(axis=0))
        pos = col_norms > 0
        if not pos.any():
            return A_target, A_other, d_in
        median_norm = max(float(np.median(col_norms[pos])), 0.01)
        clipped = np.maximum(col_norms, median_norm * 0.1)
        fs = (median_norm / clipped) * fs0

        # Interior: curvature stencil [-1, 2, -1]
        cur = np.array([-1.0, 2.0, -1.0])
        smoo_rows_t, smoo_rows_o = [], []
        for j in range(1, n_cols - 1):
            if fs[j] > 0:
                row = np.zeros(n_cols)
                row[j - 1 : j + 2] = cur * fs[j]
                smoo_rows_t.append(row)
                smoo_rows_o.append(np.zeros(A_other.shape[1]))

        # Boundaries: slope constraint [2,-2] / [-2,2]
        if n_cols >= 2:
            if fs[0] > 0:
                row = np.zeros(n_cols)
                row[0:2] = np.array([2.0, -2.0]) * fs[0]
                smoo_rows_t.append(row)
                smoo_rows_o.append(np.zeros(A_other.shape[1]))
            if fs[-1] > 0:
                row = np.zeros(n_cols)
                row[-2:] = np.array([-2.0, 2.0]) * fs[-1]
                smoo_rows_t.append(row)
                smoo_rows_o.append(np.zeros(A_other.shape[1]))

        if not smoo_rows_t:
            return A_target, A_other, d_in

        block_t = np.array(smoo_rows_t)
        block_o = np.array(smoo_rows_o)
        return (
            np.vstack([A_target, block_t]),
            np.vstack([A_other, block_o]),
            np.concatenate([d_in, np.zeros(len(smoo_rows_t))]),
        )

    # Smooth ocean velocity, then CTD velocity (MATLAB calls lainsmoo twice)
    A_ocean, A_ctd, d = _smoo_one(A_ocean, A_ctd, d, smoofac)
    A_ctd, A_ocean, d = _smoo_one(A_ctd, A_ocean, d, smoofac)
    return A_ocean, A_ctd, d


def _add_zero_mean(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    weight: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Constrain mean(u_ocean) = 0 when no external velocity reference (lainocean).

    Appends one row: sum(A_ocean columns) * weight / n_zbins = 0.
    """
    n_zbins = A_ocean.shape[1]
    scale = float(np.mean(np.abs(A_ocean).sum(axis=0)))
    row_o = np.ones(n_zbins) * weight * scale / n_zbins
    row_c = np.zeros(A_ctd.shape[1])
    return (
        np.vstack([A_ocean, row_o[np.newaxis, :]]),
        np.vstack([A_ctd, row_c[np.newaxis, :]]),
        np.concatenate([d, [0.0]]),
    )


def _add_bottom_track(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    bvel: np.ndarray,
    bvels: np.ndarray,
    *,
    botfac: float = 1.0,
    velerr: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Constrain CTD velocity to bottom-track velocity where available (lainbott).

    Each valid ensemble e gets one row: A_ctd[new_row, e] = weight_e
    with d[new_row] = bvel[e] * weight_e.
    Weight = botfac * velerr / bvels[e] scaled by sqrt(sum(|A_ctd columns|)).

    Parameters
    ----------
    bvel  : (n_se,) bottom-track velocity component (NaN = no measurement)
    bvels : (n_se,) bottom-track velocity std (m/s)
    """
    n_se = A_ctd.shape[1]
    valid = np.isfinite(bvel) & np.isfinite(bvels) & (bvels > 0)
    if not valid.any():
        return A_ocean, A_ctd, d

    col_scale = np.sqrt(np.abs(A_ctd).sum(axis=0))  # (n_se,)

    rows_o, rows_c, rhs = [], [], []
    for e in np.where(valid)[0]:
        weight_e = botfac * (velerr / bvels[e]) * col_scale[e]
        row_c = np.zeros(n_se)
        row_c[e] = weight_e
        rows_c.append(row_c)
        rows_o.append(np.zeros(A_ocean.shape[1]))
        rhs.append(bvel[e] * weight_e)

    return (
        np.vstack([A_ocean, rows_o]),
        np.vstack([A_ctd, rows_c]),
        np.concatenate([d, rhs]),
    )


def _add_sadcp(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    *,
    sadcp_z: np.ndarray,
    sadcp_vel: np.ndarray,
    sadcp_err: np.ndarray,
    dz: float,
    sadcpfac: float = 1.0,
    velerr: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Constrain u_ocean at SADCP depth bins (lainsadcp from getinv.m).

    Each valid SADCP measurement at depth z_j gets one row in A_ocean at
    column round(z_j / dz) - 1 with weight sadcpfac * velerr / sadcp_err[j].

    Parameters
    ----------
    sadcp_z   : (n_sadcp,) positive depth m
    sadcp_vel : (n_sadcp,) velocity component m/s
    sadcp_err : (n_sadcp,) velocity std m/s
    """
    n_zbins = A_ocean.shape[1]
    valid = np.isfinite(sadcp_z) & np.isfinite(sadcp_vel) & np.isfinite(sadcp_err)
    if not valid.any():
        return A_ocean, A_ctd, d

    col_scale = np.sqrt(np.abs(A_ocean).sum(axis=0))  # (n_zbins,)

    rows_o, rows_c, rhs = [], [], []
    for k in np.where(valid)[0]:
        j = int(np.round(sadcp_z[k] / dz)) - 1
        j = min(max(j, 0), n_zbins - 1)
        w = sadcpfac * (velerr / max(sadcp_err[k], 1e-6)) * col_scale[j]
        row_o = np.zeros(n_zbins)
        row_o[j] = w
        rows_o.append(row_o)
        rows_c.append(np.zeros(A_ctd.shape[1]))
        rhs.append(sadcp_vel[k] * w)

    return (
        np.vstack([A_ocean, rows_o]),
        np.vstack([A_ctd, rows_c]),
        np.concatenate([d, rhs]),
    )


def _add_barotropic(
    A_ocean_u: np.ndarray,
    A_ctd_u: np.ndarray,
    d_u: np.ndarray,
    A_ocean_v: np.ndarray,
    A_ctd_v: np.ndarray,
    d_v: np.ndarray,
    *,
    u_ship: float,
    v_ship: float,
    dt: np.ndarray,
    barofac: float = 1.0,
    nav_error: float = 30.0,
    velerr: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Constrain time-mean CTD velocity = GPS-derived ship velocity (lainbaro).

    Appends one row: sum(A_ctd * dt / T) = u_ship.

    Weight matches MATLAB getinv.m:
      barvelerr = 2 * nav_error / T          (velocity uncertainty from GPS positioning)
      fac_nav   = velerr / barvelerr         ("normalized barotropic constraint weight")
      fac       = sqrt(sum(|A_ctd|))         (column-sum scale, MATLAB lainbaro fac)
      row_weight = barofac * fac_nav * fac
    """
    T = float(dt.sum())
    # MATLAB: barvelerr = 2 * nav_error / dt_profile  (factor of 2 = error at both endpoints)
    barvelerr = 2.0 * nav_error / T
    fac_nav = velerr / barvelerr          # "normalized barotropic constraint weight"
    # MATLAB lainbaro: fac = sqrt(sum(abs(Ac)))  — sqrt of total absolute sum (no inner sqrt)
    fac = float(np.sqrt(np.abs(A_ctd_u).sum()))

    w = barofac * fac_nav * fac

    # Barotropic row: dt[e]/T per column, weighted
    row_c = dt / T * w
    row_o = np.zeros(A_ocean_u.shape[1])

    A_ocean_u2 = np.vstack([A_ocean_u, row_o[np.newaxis, :]])
    A_ctd_u2 = np.vstack([A_ctd_u, row_c[np.newaxis, :]])
    d_u2 = np.concatenate([d_u, [-u_ship * w]])

    A_ocean_v2 = np.vstack([A_ocean_v, row_o[np.newaxis, :]])
    A_ctd_v2 = np.vstack([A_ctd_v, row_c[np.newaxis, :]])
    d_v2 = np.concatenate([d_v, [-v_ship * w]])

    return A_ocean_u2, A_ctd_u2, d_u2, A_ocean_v2, A_ctd_v2, d_v2


def _solve_lsq(
    A: np.ndarray,
    d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the least-squares system d = A*m (replicates lesqfit() + lainsolv()).

    Returns
    -------
    m  : (n_params,) solution vector
    me : (n_params,) 1-sigma parameter error estimates

    Error formula matches MATLAB lesqfit:
        me = sqrt(diag(inv(A'A)) * ||d - Am||² / (n - p))
    """
    m, _, _, _ = scipy.linalg.lstsq(A, d, check_finite=False)

    # Error estimate via normal equations
    dm = A @ m
    n, p = A.shape
    dof = max(n - p, 1)
    sigma2 = float(np.sum((d - dm) ** 2) / dof)
    try:
        AtA_inv = np.linalg.inv(A.T @ A)
        me = np.sqrt(np.abs(np.diag(AtA_inv)) * sigma2)
    except np.linalg.LinAlgError:
        me = np.full(p, np.nan)

    return m, me


@dataclass
class InverseParams:
    """Tuning parameters for compute_inverse() (getinv.m ps struct)."""
    dz: float = 10.0          # depth bin size m
    botfac: float = 1.0       # bottom-track constraint weight (0 = disable)
    sadcpfac: float = 1.0     # SADCP constraint weight (0 = disable); wired in Task 7
    barofac: float = 1.0      # GPS barotropic constraint weight (0 = disable)
    smoofac: float = 0.0      # curvature smoothing weight (0 = disabled, no rows added)
    velerr: float = 0.05      # nominal velocity error m/s
    weightmin: float = 0.05   # minimum observation weight threshold
    nav_error: float = 30.0   # navigation error m (for barvelerr computation)
    down_up: bool = True       # also solve down-cast and up-cast separately


@dataclass
class InverseResult:
    """Output of compute_inverse() (getinv.m dr struct)."""
    z: np.ndarray       # (n_zbins,) depth m (positive, increasing downward)
    u: np.ndarray       # (n_zbins,) eastward velocity m/s
    v: np.ndarray       # (n_zbins,) northward velocity m/s
    uerr: np.ndarray    # (n_zbins,) velocity error estimate m/s
    nvel: np.ndarray    # (n_zbins,) number of observations per depth bin
    u_do: np.ndarray    # (n_zbins,) downcast-only eastward velocity
    v_do: np.ndarray    # (n_zbins,) downcast-only northward velocity
    u_up: np.ndarray    # (n_zbins,) upcast-only eastward velocity
    v_up: np.ndarray    # (n_zbins,) upcast-only northward velocity
    u_ctd: np.ndarray   # (n_se,) CTD eastward velocity m/s
    v_ctd: np.ndarray   # (n_se,) CTD northward velocity m/s
    ubar: float         # depth-mean eastward velocity
    vbar: float         # depth-mean northward velocity
    zctd: np.ndarray    # (n_se,) CTD depth time series m
    wctd: np.ndarray    # (n_se,) CTD vertical velocity m/s


def compute_inverse(
    se: SuperEnsemble,
    *,
    params: InverseParams | None = None,
    u_ship: float | None = None,
    v_ship: float | None = None,
    sadcp_z: np.ndarray | None = None,
    sadcp_u: np.ndarray | None = None,
    sadcp_v: np.ndarray | None = None,
    sadcp_err: np.ndarray | None = None,
) -> InverseResult:
    """Solve the LADCP inverse velocity problem (replicates getinv.m).

    Parameters
    ----------
    se      : SuperEnsemble from prepare_superensembles().
    params  : Tuning parameters; defaults to InverseParams().
    u_ship  : Mean eastward ship velocity m/s (from GPS). None = no GPS available.
              Pass 0.0 for a stationary ship — that is still a valid GPS constraint.
    v_ship  : Mean northward ship velocity m/s (from GPS). None = no GPS available.
    sadcp_z : SADCP depth m (positive) or None.
    sadcp_u : SADCP eastward velocity m/s or None.
    sadcp_v : SADCP northward velocity m/s or None.
    sadcp_err : SADCP velocity std m/s or None.
    """
    if params is None:
        params = InverseParams()

    n_se = se.izm.shape[1]

    # --- Flatten observations and build A matrices ---
    d_u, d_v, izv, jprof, wm = _flatten_obs(se, params.velerr, params.weightmin)
    if len(d_u) < 10:
        raise ValueError("Too few valid observations for inversion")

    A_ocean_sp = _build_obs_matrix(izv, params.dz)
    A_ctd_sp = _build_ctd_matrix(jprof, n_se)
    n_zbins = A_ocean_sp.shape[1]

    A_o_u, A_c_u, dw_u, idx_down, idx_up = _apply_weights(
        A_ocean_sp, A_ctd_sp, d_u, wm
    )
    A_o_v, A_c_v, dw_v, _, _ = _apply_weights(A_ocean_sp, A_ctd_sp, d_v, wm)

    # --- Depth vector for output ---
    z = np.arange(1, n_zbins + 1) * params.dz  # positive, 1-indexed depth bins

    # --- Smoothness constraints (applied to both U and V identically) ---
    A_o_u, A_c_u, dw_u = _add_smoothness(A_o_u, A_c_u, dw_u, params.smoofac)
    A_o_v, A_c_v, dw_v = _add_smoothness(A_o_v, A_c_v, dw_v, params.smoofac)

    # --- Bottom-track constraint ---
    has_btrack = params.botfac > 0 and np.any(np.isfinite(se.bvel[0]))
    if has_btrack:
        A_o_u, A_c_u, dw_u = _add_bottom_track(
            A_o_u, A_c_u, dw_u,
            se.bvel[0], se.bvels[0], botfac=params.botfac, velerr=params.velerr,
        )
        A_o_v, A_c_v, dw_v = _add_bottom_track(
            A_o_v, A_c_v, dw_v,
            se.bvel[1], se.bvels[1], botfac=params.botfac, velerr=params.velerr,
        )

    # --- SADCP (ship ADCP) constraint ---
    has_sadcp = (params.sadcpfac > 0 and sadcp_z is not None
                 and np.any(np.isfinite(sadcp_z)))
    if has_sadcp:
        A_o_u, A_c_u, dw_u = _add_sadcp(
            A_o_u, A_c_u, dw_u,
            sadcp_z=sadcp_z, sadcp_vel=sadcp_u, sadcp_err=sadcp_err,
            dz=params.dz, sadcpfac=params.sadcpfac, velerr=params.velerr,
        )
        A_o_v, A_c_v, dw_v = _add_sadcp(
            A_o_v, A_c_v, dw_v,
            sadcp_z=sadcp_z, sadcp_vel=sadcp_v, sadcp_err=sadcp_err,
            dz=params.dz, sadcpfac=params.sadcpfac, velerr=params.velerr,
        )

    # --- Barotropic (GPS) constraint ---
    # Gate on whether GPS was provided (u_ship is not None), not on ship speed.
    # A stationary ship (u_ship=0.0) is still a valid GPS-derived constraint.
    has_baro = params.barofac > 0 and u_ship is not None
    if has_baro:
        A_o_u, A_c_u, dw_u, A_o_v, A_c_v, dw_v = _add_barotropic(
            A_o_u, A_c_u, dw_u, A_o_v, A_c_v, dw_v,
            u_ship=float(u_ship), v_ship=float(v_ship),  # type: ignore[arg-type]
            dt=se.dt, barofac=params.barofac,
            nav_error=params.nav_error, velerr=params.velerr,
        )

    # --- Zero-mean fallback when no external constraint ---
    if not has_btrack and not has_baro and not has_sadcp:
        A_o_u, A_c_u, dw_u = _add_zero_mean(A_o_u, A_c_u, dw_u)
        A_o_v, A_c_v, dw_v = _add_zero_mean(A_o_v, A_c_v, dw_v)

    # --- Solve full-cast system ---
    A_full_u = np.hstack([A_o_u, A_c_u])
    A_full_v = np.hstack([A_o_v, A_c_v])
    m_u, me_u = _solve_lsq(A_full_u, dw_u)
    m_v, me_v = _solve_lsq(A_full_v, dw_v)

    u_ocean = m_u[:n_zbins]
    v_ocean = m_v[:n_zbins]
    u_ctd_neg = m_u[n_zbins:]
    v_ctd_neg = m_v[n_zbins:]
    # Sign convention: solved u_ctd_neg = -u_CTD  (see MATLAB dr.uctd = -real(uctd))
    u_ctd = -u_ctd_neg[:n_se] if len(u_ctd_neg) >= n_se else np.full(n_se, np.nan)
    v_ctd = -v_ctd_neg[:n_se] if len(v_ctd_neg) >= n_se else np.full(n_se, np.nan)

    uerr = np.sqrt(me_u[:n_zbins] ** 2 + me_v[:n_zbins] ** 2)
    nvel = np.asarray(A_ocean_sp.sum(axis=0)).ravel()

    # --- Down/up cast separately (ps.down_up=1) ---
    _BAROCLINIC_FAC = 10.0  # MATLAB: baroclinfac = 10 (large = forces zero mean)

    def _solve_subset(idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # idx rows are from the pre-augmentation observation block; safe because
        # constraint rows are appended after _apply_weights returns idx_down/idx_up.
        A_os = A_o_u[idx, :n_zbins]
        A_cs = A_c_u[idx, :]
        ds_u = dw_u[idx]
        A_os2 = A_o_v[idx, :n_zbins]
        A_cs2 = A_c_v[idx, :]
        ds_v = dw_v[idx]

        # Remove zero-constrained CTD columns
        active = np.where(np.abs(A_cs).sum(axis=0) > 0)[0]
        A_cs = A_cs[:, active]
        A_cs2 = A_cs2[:, active]

        # Add zero-mean baroclinic constraint
        A_os = np.vstack([A_os, (_BAROCLINIC_FAC * np.ones(n_zbins))[np.newaxis, :]])
        A_cs = np.vstack([A_cs, np.zeros((1, A_cs.shape[1]))])
        ds_u = np.concatenate([ds_u, [0.0]])
        A_os2 = np.vstack([A_os2, (_BAROCLINIC_FAC * np.ones(n_zbins))[np.newaxis, :]])
        A_cs2 = np.vstack([A_cs2, np.zeros((1, A_cs2.shape[1]))])
        ds_v = np.concatenate([ds_v, [0.0]])

        if A_os.shape[0] < A_os.shape[1] + 2:
            return np.full(n_zbins, np.nan), np.full(n_zbins, np.nan)

        mu, _ = _solve_lsq(np.hstack([A_os, A_cs]), ds_u)
        mv, _ = _solve_lsq(np.hstack([A_os2, A_cs2]), ds_v)
        u_sub = mu[:n_zbins]
        v_sub = mv[:n_zbins]
        # Clip unreasonably large values
        u_sub[np.abs(u_sub) > 5.0] = np.nan
        v_sub[np.abs(v_sub) > 5.0] = np.nan
        return u_sub, v_sub

    if params.down_up and len(idx_down) > 5 and len(idx_up) > 5:
        u_do, v_do = _solve_subset(idx_down)
        u_up, v_up = _solve_subset(idx_up)
    else:
        u_do = u_ocean.copy()
        v_do = v_ocean.copy()
        u_up = u_ocean.copy()
        v_up = v_ocean.copy()

    return InverseResult(
        z=z,
        u=u_ocean, v=v_ocean, uerr=uerr, nvel=nvel,
        u_do=u_do, v_do=v_do, u_up=u_up, v_up=v_up,
        u_ctd=u_ctd, v_ctd=v_ctd,
        ubar=float(np.nanmean(u_ocean)),
        vbar=float(np.nanmean(v_ocean)),
        zctd=se.z,
        wctd=-np.nanmean(se.rw, axis=0),
    )
