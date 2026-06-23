from __future__ import annotations

import dataclasses
import math

import numpy as np

from ladcp.solution.inverse import EnsembleData


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
