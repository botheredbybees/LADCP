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
    n: NDArray           # (nz,) int — estimates used per bin after editing
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
    raise NotImplementedError
