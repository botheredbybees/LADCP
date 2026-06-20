"""Shear-based horizontal velocity solution.

Reference: docs/legacy/getshear2.m (Visbeck, LDEO 1997).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ladcp._typing import NDArray


def _central_diff_shear(
    u: NDArray,
    v: NDArray,
    w: NDArray,
    izm: NDArray,
    weight_mask: NDArray,
) -> tuple[NDArray, NDArray, NDArray]:
    """Compute stride-2 central-difference shear, weighted and NaN-padded.

    Replicates MATLAB diff2(): x[2:,:] - x[:-2,:] divided by izm[2:] - izm[:-2].
    First and last rows are NaN because no two-sided neighbour exists.
    weight_mask is 1.0 where valid, NaN where excluded.
    """
    du = u[2:, :] - u[:-2, :]     # (nbin-2, nens)
    dv = v[2:, :] - v[:-2, :]
    dw = w[2:, :] - w[:-2, :]
    dz = izm[2:, :] - izm[:-2, :]  # depth increment between skip-one neighbours

    shear_u = np.full_like(u, np.nan)
    shear_v = np.full_like(v, np.nan)
    shear_w = np.full_like(w, np.nan)

    shear_u[1:-1, :] = du / dz
    shear_v[1:-1, :] = dv / dz
    shear_w[1:-1, :] = dw / dz

    shear_u *= weight_mask
    shear_v *= weight_mask
    shear_w *= weight_mask

    return shear_u, shear_v, shear_w


def _bin_average_shear(
    shear_u: NDArray,
    shear_v: NDArray,
    shear_w: NDArray,
    izm: NDArray,
    z_bins: NDArray,
    dz: float,
    stdf: float = 2.0,
) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray, NDArray]:
    """Average shear estimates into depth bins with 2σ outlier editing.

    For each output bin centred at z_bins[k]:
      1. Collect all (bin, ensemble) pairs where depth is within dz of z_bins[k]+dz/2.
      2. Discard NaN estimates.
      3. If n ≤ 2: leave mean/std as NaN.
      4. Reject estimates more than stdf*std from the median.
      5. If ≥ 2 remain: store mean and std.

    Window condition replicates MATLAB: abs(izm - (z + dz/2)) <= dz.
    """
    nz = len(z_bins)
    iz_flat = izm.ravel()
    su_flat = shear_u.ravel()
    sv_flat = shear_v.ravel()
    sw_flat = shear_w.ravel()

    usm = np.full(nz, np.nan)
    vsm = np.full(nz, np.nan)
    wsm = np.full(nz, np.nan)
    use = np.full(nz, np.nan)
    vse = np.full(nz, np.nan)
    wse = np.full(nz, np.nan)
    nn = np.zeros(nz, dtype=np.intp)

    for k, center in enumerate(z_bins):
        in_window = np.abs(iz_flat - (center + dz / 2)) <= dz
        finite = in_window & np.isfinite(su_flat + sv_flat)
        su = su_flat[finite]
        sv = sv_flat[finite]
        sw = sw_flat[finite & np.isfinite(sw_flat)]
        n = len(su)
        nn[k] = n
        if n <= 2:
            continue

        for arr, mean_out, std_out in [
            (su, usm, use),
            (sv, vsm, vse),
        ]:
            med = np.median(arr)
            std = np.std(arr, ddof=1)
            # <= not < so std=0 (constant window) keeps all values rather than rejecting all
            keep = np.abs(arr - med) <= stdf * std
            if keep.sum() > 1:
                mean_out[k] = np.mean(arr[keep])
                std_out[k] = np.std(arr[keep], ddof=1)

        # w may have fewer finite values — filter independently
        if len(sw) > 2:
            med_w = np.median(sw)
            std_w = np.std(sw, ddof=1)
            # <= not < so std=0 (constant window) keeps all values rather than rejecting all
            keep_w = np.abs(sw - med_w) <= stdf * std_w
            if keep_w.sum() > 1:
                wsm[k] = np.mean(sw[keep_w])
                wse[k] = np.std(sw[keep_w], ddof=1)

    return usm, vsm, wsm, use, vse, wse, nn


def _integrate_shear(
    usm: NDArray,
    vsm: NDArray,
    wsm: NDArray,
    dz: float,
) -> tuple[NDArray, NDArray, NDArray]:
    """Integrate shear from bottom up to produce zero-mean relative velocities.

    NaN bins are replaced with 0 (no shear contribution) before integrating.
    Replicates MATLAB: flipud(cumsum(flipud(usm))) * dz, then subtract mean.
    """
    u = np.where(np.isnan(usm), 0.0, usm)
    v = np.where(np.isnan(vsm), 0.0, vsm)
    w = np.where(np.isnan(wsm), 0.0, wsm)

    ur = np.flipud(np.cumsum(np.flipud(u))) * dz
    vr = np.flipud(np.cumsum(np.flipud(v))) * dz
    wr = np.flipud(np.cumsum(np.flipud(w))) * dz

    ur -= np.mean(ur)
    vr -= np.mean(vr)
    wr -= np.mean(wr)

    return ur, vr, wr


@dataclass
class ShearProfile:
    """Depth-binned shear profile and integrated relative velocity.

    All arrays have shape (nz,) where nz = ceil(z_max / dz).
    z[k] is the centre of the k-th depth bin in metres (positive down).
    """

    z: NDArray           # (nz,) bin centre depths, m
    u_shear: NDArray     # (nz,) mean du/dz per bin, s⁻¹
    v_shear: NDArray     # (nz,) mean dv/dz per bin, s⁻¹
    w_shear: NDArray     # (nz,) mean dw/dz per bin, s⁻¹
    u_shear_err: NDArray # (nz,) 1-σ shear uncertainty, s⁻¹
    v_shear_err: NDArray # (nz,) 1-σ shear uncertainty, s⁻¹
    w_shear_err: NDArray # (nz,) 1-σ shear uncertainty, s⁻¹
    n: NDArray           # (nz,) int — estimates in window before outlier editing
    u_rel: NDArray       # (nz,) relative eastward velocity, m/s (zero-mean)
    v_rel: NDArray       # (nz,) relative northward velocity, m/s (zero-mean)
    w_rel: NDArray       # (nz,) relative vertical velocity, m/s (zero-mean)


def compute_shear(
    u: NDArray,
    v: NDArray,
    w: NDArray,
    izm: NDArray,
    weight: NDArray,
    dz: float = 10.0,
    *,
    stdf: float = 2.0,
    weight_min: float = 0.1,
) -> ShearProfile:
    """Compute depth-binned shear profile and integrate to relative velocities.

    Parameters
    ----------
    u, v, w : ndarray, shape (nbin, nens)
        Earth-frame velocity components (East, North, Up) in m/s.
    izm : ndarray, shape (nbin, nens)
        Depth of each ADCP bin in metres (positive down).
    weight : ndarray, shape (nbin, nens)
        Quality weight in [0, 1]. Bins with weight ≤ weight_min are excluded.
    dz : float
        Output depth-bin size in metres. Default 10 m.
    stdf : float
        Outlier rejection threshold in units of std. Default 2 (2σ editing).
    weight_min : float
        Minimum weight for a bin to be included. Default 0.1.

    Returns
    -------
    ShearProfile
        Depth-binned shear and integrated relative velocity profile.
    """
    weight_mask = np.where(weight > weight_min, 1.0, np.nan)

    shear_u, shear_v, shear_w = _central_diff_shear(u, v, w, izm, weight_mask)

    z_max = float(np.nanmax(izm))
    z_bins = np.arange(dz / 2, z_max, dz)

    usm, vsm, wsm, use, vse, wse, nn = _bin_average_shear(
        shear_u, shear_v, shear_w, izm, z_bins, dz, stdf
    )

    ur, vr, wr = _integrate_shear(usm, vsm, wsm, dz)

    return ShearProfile(
        z=z_bins,
        u_shear=usm, v_shear=vsm, w_shear=wsm,
        u_shear_err=use, v_shear_err=vse, w_shear_err=wse,
        n=nn,
        u_rel=ur, v_rel=vr, w_rel=wr,
    )
