from __future__ import annotations

import dataclasses
import math

import numpy as np

from ladcp.solution.inverse import EnsembleData, _ref_medianan


def edit_sidelobes(
    ens: EnsembleData,
    *,
    zbottom: float | None = None,
    theta_deg: float = 20.0,
    cell_size_m: float,
) -> EnsembleData:
    """Zero-weight ADCP bins contaminated by surface and bottom acoustic side-lobes.

    Based on LDEO_IX edit_data.m lines 142–186 (Eric Firing convention).
    Unlike edit_data.m, the surface edit is applied unconditionally — the MATLAB
    reference only applies it for uplooker configurations, but surface sidelobe
    contamination is real regardless. The 1.5× cell_size margin keeps the edit
    conservative. Returns a new EnsembleData; the input is not modified.
    """
    f = 1.0 - math.cos(math.radians(theta_deg))
    margin = 1.5 * cell_size_m

    if zbottom is None:
        derived = float(np.nanmax(-ens.z + ens.hbot))
        zbottom = derived if math.isfinite(derived) else None

    # Surface mask: bins shallower than zlim_surface are contaminated.
    # zlim_surface shape: (n_ens,); ens.izm shape: (n_bins, n_ens) — broadcasts correctly.
    zlim_surface = f * ens.z - margin
    bad_surface = ens.izm > zlim_surface

    if zbottom is not None:
        hab = ens.z + zbottom                       # (n_ens,) height above floor
        zlim_bot = -zbottom + f * hab + margin      # (n_ens,)
        bad_bottom = ens.izm < zlim_bot
    else:
        bad_bottom = np.zeros_like(ens.izm, dtype=bool)

    new_weight = ens.weight.copy()
    new_weight[bad_surface | bad_bottom] = np.nan
    return dataclasses.replace(ens, weight=new_weight)


def edit_large_velocities(
    ens: EnsembleData,
    *,
    maxspeed: float = 2.5,
) -> EnsembleData:
    """Mask bins where horizontal speed exceeds maxspeed.

    Mirrors MATLAB loadrdi.m lines 235–242: sqrt(u²+v²) > p.vlim=2.5 m/s sets
    that bin's weight to NaN. Applied after Earth-frame rotation, before
    super-ensemble averaging.
    """
    speed_sq = np.where(np.isfinite(ens.u) & np.isfinite(ens.v),
                        ens.u ** 2 + ens.v ** 2,
                        0.0)
    bad = speed_sq > maxspeed ** 2
    new_weight = ens.weight.copy()
    new_weight[bad] = np.nan
    return dataclasses.replace(ens, weight=new_weight)


def _rms(x: np.ndarray) -> float:
    """MATLAB outlier.m rms(): root-mean-square without mean removal, NaN-aware."""
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(x * x)))


def tilt_from_pitch_roll(
    pitch_deg: np.ndarray, roll_deg: np.ndarray
) -> np.ndarray:
    """Combined tilt angle in degrees; loadrdi.m lines 414-415."""
    pit_r = np.radians(pitch_deg)
    rol_r = np.radians(roll_deg)
    return np.degrees(np.arcsin(np.clip(
        np.sqrt(np.sin(pit_r) ** 2 + np.sin(rol_r) ** 2), 0.0, 1.0)))


def edit_outliers(
    ens: EnsembleData,
    *,
    n: tuple[float, ...] = (4.0, 3.0),
    block_minutes: float = 5.0,
    do_bvel: bool = True,
    nblock: int | None = None,
) -> EnsembleData:
    """Reject velocity outliers per 5-minute block; port of loadrdi.m::outlier().

    For the DL and UL bin blocks independently: subtract the per-ensemble
    median (medianan) from u/v/w, then within each block of ~block_minutes
    NaN every cell whose anomaly exceeds n[sweep] * rms(block anomalies).
    Two sweeps by default (4 sigma then 3 sigma), the second computed on the
    already-cleaned fields. Unlike the weight-only editors, this NaNs the
    velocity components themselves as well as the weight — matching
    outlier.m lines 114-117 (d.ru/rv/rw += dummy).

    Bottom track (DL section of outlier.m): per block, NaN whole bvel rows
    (and hbot) where the median-removed bvel u/v, the w anomaly (bvel w minus
    the water-column median w), or the median-removed hbot exceeds the same
    threshold. outlier.m gates this on a 4-column bvel; our bvel is (nens, 3)
    so it is applied unconditionally.

    Returns a new EnsembleData; the input is not modified.
    """
    if nblock is None:
        dt_min = float(np.nanmean(np.diff(ens.time_jul))) * 24.0 * 60.0
        nblock = max(1, int(math.ceil(block_minutes / dt_min)))
    n_ens = ens.u.shape[1]
    starts = range(0, n_ens, nblock)

    u = ens.u.copy()
    v = ens.v.copy()
    w = ens.w.copy()
    weight = ens.weight.copy()
    bvel = ens.bvel.copy()
    hbot = ens.hbot.copy()

    for rows, do_bvel_block in ((ens.izd, do_bvel), (ens.izu, False)):
        if len(rows) == 0:
            continue
        ru = u[rows, :].copy()
        rv = v[rows, :].copy()
        rw = w[rows, :].copy()
        dummy = np.zeros_like(rw)
        bdummy = np.zeros_like(bvel) if do_bvel_block else None
        bv = bvel.copy() if do_bvel_block else None

        for ni in n:
            rwm = _ref_medianan(rw)
            rw = rw - rwm[np.newaxis, :]
            ru = ru - _ref_medianan(ru)[np.newaxis, :]
            rv = rv - _ref_medianan(rv)[np.newaxis, :]
            if do_bvel_block:
                bv[:, 2] = bv[:, 2] - rwm
            for s in starts:
                sel = slice(s, min(s + nblock, n_ens))
                for anom in (rw[:, sel], ru[:, sel], rv[:, sel]):
                    bad = np.abs(anom) > ni * _rms(anom)
                    dummy[:, sel][bad] = np.nan
                if do_bvel_block:
                    bu_a = bv[sel, 0] - np.nanmedian(bv[sel, 0])
                    bv[sel, 0] = bu_a
                    bdummy[sel][np.abs(bu_a) > ni * _rms(bu_a)] = np.nan
                    bv_a = bv[sel, 1] - np.nanmedian(bv[sel, 1])
                    bv[sel, 1] = bv_a
                    bdummy[sel][np.abs(bv_a) > ni * _rms(bv_a)] = np.nan
                    bw_a = bv[sel, 2]
                    bdummy[sel][np.abs(bw_a) > ni * _rms(bw_a)] = np.nan
                    hb_a = hbot[sel] - np.nanmedian(hbot[sel])
                    bdummy[sel][np.abs(hb_a) > ni * _rms(hb_a)] = np.nan
            # propagate this sweep's rejections into the next sweep's stats
            rw = rw + dummy
            ru = ru + dummy
            rv = rv + dummy

        weight[rows, :] = weight[rows, :] + dummy
        u[rows, :] = u[rows, :] + dummy
        v[rows, :] = v[rows, :] + dummy
        w[rows, :] = w[rows, :] + dummy
        if do_bvel_block:
            bvel = bvel + bdummy
            hbot = hbot + bdummy[:, 0]

    return dataclasses.replace(
        ens, u=u, v=v, w=w, weight=weight, bvel=bvel, hbot=hbot
    )


def build_ldeo_weights(
    cm: np.ndarray,
    ts: np.ndarray,
    pitch_dl: np.ndarray,
    roll_dl: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    izd: np.ndarray,
    izu: np.ndarray,
    *,
    tiltmax: tuple[float, float] = (22.0, 4.0),
    weighbin1: float = 1.0,
) -> np.ndarray:
    """Build the per-cell weight field; port of loadrdi.m lines 408-533.

    Steps, in loadrdi.m order:
      1. weight = cm (median-over-beams correlation, combined UL+DL rows),
         normalized by medianan(maxnan(weight)) -- the median over ensembles
         of the per-ensemble max.
      2. NaN whole ensembles whose tilt exceeds tiltmax[0] degrees or whose
         tilt derivative exceeds tiltmax[1] (tilt from DL pitch/roll only,
         matching l.pit(1,:)/l.rol(1,:)).
      3. Echo-amplitude penalty per row: where ts exceeds the row median,
         weight *= (1 - (ts_anom/max(ts_anom))**1.5)  (crosstalk / bottom
         echo suppression).
      4. Non-pinging ensemble removal per instrument block: ensembles whose
         median |bin-to-bin gradient| of w AND v are both < 0.005 m/s.
         NOTE: loadrdi.m line 478 tests `dru` twice and `drv` never (a real
         bug -- only the rw- and rv-gradients are effective); replicated
         here for parity.
      5. Multiply bin-1 weight by weighbin1 (default 1: no-op).

    Args:
        cm: (nbin, nens) median-over-beams correlation, combined array.
        ts: (nbin, nens) median-over-beams echo amplitude, combined array.
        pitch_dl, roll_dl: (nens,) DL attitude in degrees.
        v, w: (nbin, nens) combined earth-frame velocities.
        izd, izu: combined-array row indices per instrument.

    Returns (nbin, nens) weight array (new; inputs unmodified).
    """
    weight = np.asarray(cm, dtype=np.float64).copy()
    col_max = np.nanmax(weight, axis=0)
    weight = weight / np.nanmedian(col_max)

    tilt = tilt_from_pitch_roll(pitch_dl, roll_dl)
    weight[:, tilt > tiltmax[0]] = np.nan

    def _edge_diff(x: np.ndarray) -> np.ndarray:
        # mean of |backward diff| and |forward diff|, zero-padded at the ends
        # (loadrdi.m: mean(abs(diff([0,x;x,0]'))')).
        back = np.abs(np.diff(np.concatenate(([0.0], x))))
        fwd = np.abs(np.diff(np.concatenate((x, [0.0]))))
        return 0.5 * (back + fwd)

    tiltd = np.sqrt(_edge_diff(roll_dl) ** 2 + _edge_diff(pitch_dl) ** 2)
    weight[:, tiltd > tiltmax[1]] = np.nan

    ts = np.asarray(ts, dtype=np.float64)
    for i in range(weight.shape[0]):
        ts_anom = ts[i] - np.nanmedian(ts[i])
        pos = ts_anom > 0
        if pos.any():
            weight[i, pos] *= 1.0 - (ts_anom[pos] / np.nanmax(ts_anom)) ** 1.5

    for rows in (izd, izu):
        if len(rows) < 2:
            continue
        gw = np.nanmedian(np.abs(np.diff(w[rows, :], axis=0)), axis=0)
        gv = np.nanmedian(np.abs(np.diff(v[rows, :], axis=0)), axis=0)
        nonping = (np.abs(gw) < 0.005) & (np.abs(gv) < 0.005)
        weight[np.ix_(rows, np.flatnonzero(nonping))] = np.nan

    if weighbin1 != 1.0:
        if len(izd) > 0:
            weight[izd[0], :] *= weighbin1
        if len(izu) > 0:
            weight[izu[0], :] *= weighbin1

    return weight


def edit_mask_bins(
    ens: EnsembleData,
    *,
    dn_bins: list[int] | tuple[int, ...] = (),
    up_bins: list[int] | tuple[int, ...] = (),
) -> EnsembleData:
    """NaN the weight of whole instrument bins; port of edit_data.m bin masking.

    dn_bins/up_bins are 0-based per-instrument bin indices (edit_data.m's
    p.edit_mask_dn_bins/up_bins are 1-based). edit_data.m auto-adds bin 1
    (our bin 0) for any instrument whose blanking distance is zero — callers
    should pass dn_bins=[0] / up_bins=[0] when rdi.blnk_m == 0. Only the
    weight is NaN'd (edit_data.m line 120/129); velocities are untouched.
    Out-of-range indices are ignored (edit_data.m checks <= nbin).
    """
    weight = ens.weight.copy()
    for bins, rows in ((dn_bins, ens.izd), (up_bins, ens.izu)):
        for b in bins:
            if 0 <= b < len(rows):
                weight[rows[b], :] = np.nan
    return dataclasses.replace(ens, weight=weight)


def edit_w_outliers(
    ens: EnsembleData,
    *,
    wlim: float = 0.2,
    wrange: int = 5,
) -> EnsembleData:
    """Mask bins with anomalous vertical velocity relative to near-instrument reference.

    Mirrors MATLAB loadrdi.m lines 185–201. For each ensemble the reference w is
    the median of the `wrange` bins nearest each transducer (DL and UL separately).
    Any bin where |w_measured - w_ref| > wlim is masked.

    Physical meaning: since the ADCP measures velocity relative to itself and the
    instrument descends at ~1 m/s, near-transducer bins give w_ref ≈ +1 m/s.
    Bins deviating by more than wlim are contaminated (reflections, vibration,
    beam interference) and should be excluded before shear-based inversion.
    """
    w_ref = np.full(ens.w.shape, np.nan)

    # DL reference: first wrange DL bins (nearest to DL transducer = shallowest DL)
    if len(ens.izd) > 0:
        dl_ref_idx = ens.izd[:wrange]
        w_ref_dl = np.nanmedian(ens.w[dl_ref_idx, :], axis=0)  # (n_ens,)
        w_ref[ens.izd, :] = w_ref_dl[np.newaxis, :]

    # UL reference: first wrange UL bins (nearest to UL transducer = deepest UL)
    if len(ens.izu) > 0:
        ul_ref_idx = ens.izu[:wrange]
        w_ref_ul = np.nanmedian(ens.w[ul_ref_idx, :], axis=0)  # (n_ens,)
        w_ref[ens.izu, :] = w_ref_ul[np.newaxis, :]

    # Mask bins deviating more than wlim from their reference
    bad = np.abs(ens.w - w_ref) > wlim
    bad &= np.isfinite(ens.w)

    new_weight = ens.weight.copy()
    new_weight[bad] = np.nan
    return dataclasses.replace(ens, weight=new_weight)
