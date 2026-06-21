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
