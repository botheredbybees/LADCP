"""Shear-based horizontal velocity solution.

Reference: docs/legacy/getshear2.m (Visbeck, LDEO 1997).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ladcp._typing import NDArray


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
