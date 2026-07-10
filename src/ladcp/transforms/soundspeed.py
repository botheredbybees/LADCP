"""Sound-speed corrections; ports of LDEO_IX sounds.m, press.m, and the
correction blocks in getdpthi.m (lines 182-207 velocities, 428-441 izm).

The ADCP converts Doppler shifts to velocities using an assumed sound speed
at the transducer (d.sv, from the fixed/variable leader). When the true
sound speed ss (from CTD pressure/temperature) differs, every velocity and
every along-beam distance is off by the ratio ss/sv, per ensemble and per
instrument. LDEO applies this in getdpthi.m; the correction is gated there
on d.soundc == 0 (i.e. the instrument did NOT already track sound speed
itself) -- callers are responsible for that check.
"""
from __future__ import annotations

import dataclasses

import numpy as np
from numpy.typing import NDArray

from ladcp.solution.inverse import EnsembleData


def sound_speed(
    p_dbar: NDArray[np.float64],
    temp_c: NDArray[np.float64],
    salinity: float | NDArray[np.float64] = 34.5,
) -> NDArray[np.float64]:
    """Sound speed in seawater (m/s); Chen & Millero (1977), UNESCO 44.

    Port of sounds.m (itself a translation of the UNESCO FORTRAN SVEL).
    Inputs: pressure in dbar, temperature in degC (IPTS-68), salinity PSS-78.
    Check value: 1731.995 m/s at S=40, T=40, P=10000 dbar. LDEO's getdpthi.m
    calls this with constant salinity 34.5 when no CTD salinity is available.
    """
    P = np.asarray(p_dbar, dtype=np.float64) / 10.0  # dbar -> bars
    T = np.asarray(temp_c, dtype=np.float64)
    S = np.asarray(salinity, dtype=np.float64)
    SR = np.sqrt(np.abs(S))
    # S**2 term
    D = 1.727e-3 - 7.8936e-6 * P
    # S**(3/2) term
    B1 = 7.3637e-5 + 1.7945e-7 * T
    B0 = -1.922e-2 - 4.42e-5 * T
    B = B0 + B1 * P
    # S**1 term
    A3 = (-3.389e-13 * T + 6.649e-12) * T + 1.100e-10
    A2 = ((7.988e-12 * T - 1.6002e-10) * T + 9.1041e-9) * T - 3.9064e-7
    A1 = (((-2.0122e-10 * T + 1.0507e-8) * T - 6.4885e-8) * T - 1.2580e-5) * T \
        + 9.4742e-5
    A0 = (((-3.21e-8 * T + 2.006e-6) * T + 7.164e-5) * T - 1.262e-2) * T + 1.389
    A = ((A3 * P + A2) * P + A1) * P + A0
    # S**0 term
    C3 = (-2.3643e-12 * T + 3.8504e-10) * T - 9.7729e-9
    C2 = (((1.0405e-12 * T - 2.5335e-10) * T + 2.5974e-8) * T - 1.7107e-6) * T \
        + 3.1260e-5
    C1 = (((-6.1185e-10 * T + 1.3621e-7) * T - 8.1788e-6) * T + 6.8982e-4) * T \
        + 0.153563
    C0 = ((((3.1464e-9 * T - 1.47800e-6) * T + 3.3420e-4) * T - 5.80852e-2) * T
          + 5.03711) * T + 1402.388
    C = ((C3 * P + C2) * P + C1) * P + C0
    return C + (A + B * SR + D * S) * S


def depth_to_pressure(depth_m: NDArray[np.float64]) -> NDArray[np.float64]:
    """Depth (m) to pressure (dbar); GEOSECS formula, port of press.m."""
    z = np.asarray(depth_m, dtype=np.float64)
    return 2.398599584e05 - np.sqrt(5.753279964e10 - 4.833657881e05 * z)


def apply_sound_speed_correction(
    ens: EnsembleData,
    *,
    ss: NDArray[np.float64],
    sv_dl: NDArray[np.float64],
    sv_ul: NDArray[np.float64] | None = None,
) -> EnsembleData:
    """Scale velocities, bottom track, and izm bin offsets by ss/sv.

    Port of getdpthi.m lines 182-207 (velocities: DL rows by ss/sv_dl, UL
    rows by ss/sv_ul; bvel and hbot by the DL ratio) and lines 428-441 (the
    bin-offset part of izm scales by the same ratio; the CTD-derived
    instrument depth z does NOT scale, so izm_new = z + (izm - z) * sc).

    Args:
        ens: combined-ensemble data (izm already assigned).
        ss: (nens,) true sound speed at the instrument.
        sv_dl / sv_ul: (nens,) sound speed each instrument assumed
            (RDIData.sound_vel_ms, UL aligned to DL columns).

    Returns a new EnsembleData; the input is not modified.
    """
    sc_dl = np.asarray(ss, dtype=np.float64) / np.asarray(sv_dl, dtype=np.float64)

    u = ens.u.copy()
    v = ens.v.copy()
    w = ens.w.copy()
    izm = ens.izm.copy()

    for rows, sc in (
        (ens.izd, sc_dl),
        (ens.izu, None if sv_ul is None
         else np.asarray(ss, dtype=np.float64) / np.asarray(sv_ul, dtype=np.float64)),
    ):
        if len(rows) == 0 or sc is None:
            continue
        u[rows, :] *= sc[np.newaxis, :]
        v[rows, :] *= sc[np.newaxis, :]
        w[rows, :] *= sc[np.newaxis, :]
        izm[rows, :] = ens.z[np.newaxis, :] + (izm[rows, :] - ens.z[np.newaxis, :]) \
            * sc[np.newaxis, :]

    bvel = ens.bvel * sc_dl[:, np.newaxis]
    hbot = ens.hbot * sc_dl

    return dataclasses.replace(ens, u=u, v=v, w=w, izm=izm, bvel=bvel, hbot=hbot)
