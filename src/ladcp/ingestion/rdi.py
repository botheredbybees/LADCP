"""Read Teledyne RDI PD0 binary files. Reference: docs/legacy/loadrdi.m."""

from pathlib import Path

import numpy as np

from ladcp.ingestion._pd0 import parse_pd0
from ladcp.ingestion._types import RDIData


def best_ul_shift(
    w_dl: np.ndarray,
    w_ul: np.ndarray,
    ul_idx: np.ndarray,
    *,
    max_lag: int = 20,
    min_corr: float = 0.9,
) -> tuple[int, float]:
    """Refine the UL->DL ensemble pairing by w cross-correlation.

    Port of loadrdi.m's merge-lag step (lines 932-956): after the
    nearest-time match (ul_idx), both instruments' per-ensemble median
    vertical velocities are cross-correlated over integer ensemble shifts;
    the best shift is applied to the pairing. Physically both instruments
    ride the same rosette, so their w series are near-identical — the
    recorded clocks, however, can differ by a substantial fraction of a
    ping interval, making nearest-RECORDED-time pick the wrong neighbor
    (P16N 003: UL clock ~0.6 s off, ping interval 1.33-1.58 s, true pairing
    is nearest-time MINUS 1 for ~80% of ensembles).

    The shift is a SEQUENCE shift, exactly like loadrdi.m's ``iu = iu(iiu)``:
    the corrected pairing for DL ensemble k is ``ul_idx[k + shift]``, NOT
    ``ul_idx[k] + shift``. With staggered pinging the nearest-time index
    sequence jumps by 0 or 2 between neighbors ~20% of the time, and only
    the sequence shift reproduces Octave's merge there (verified on P16N
    003 by heading fingerprint: sequence shift -1 matches 100.0% of merged
    ensembles; the value shift only 78%).

    Args:
        w_dl: (nbin_dl, nens_dl) DL earth-frame vertical velocity.
        w_ul: (nbin_ul, nens_ul) UL earth-frame vertical velocity.
        ul_idx: (nens_dl,) nearest-time UL ensemble index per DL ensemble.
        max_lag: shift search half-window (loadrdi.m maxlag = 20).
        min_corr: below this, keep the nearest-time pairing
            (loadrdi.m: "best lag not obvious").

    Returns:
        (shift, corr): the corrected pairing is
        ``ul_idx[np.clip(np.arange(len(ul_idx)) + shift, 0, len(ul_idx)-1)]``;
        corr is the w correlation at that shift.
    """
    wb_d = np.nanmedian(w_dl, axis=0)
    wb_u_full = np.nanmedian(w_ul, axis=0)
    n_dl = len(ul_idx)
    k = np.arange(n_dl)

    best = (0, -np.inf)
    for s in range(-max_lag, max_lag + 1):
        wb_u = wb_u_full[ul_idx[np.clip(k + s, 0, n_dl - 1)]]
        ok = np.isfinite(wb_d) & np.isfinite(wb_u)
        if ok.sum() < 10:
            continue
        a = wb_d[ok] - wb_d[ok].mean()
        b = wb_u[ok] - wb_u[ok].mean()
        denom = np.sqrt(np.dot(a, a) * np.dot(b, b))
        if denom == 0.0:
            continue
        c = float(np.dot(a, b) / denom)
        if c > best[1]:
            best = (s, c)

    shift, corr = best
    if corr < min_corr:
        return 0, corr
    return shift, corr


def load_rdi(path: Path) -> RDIData:
    """Load one RDI PD0 binary (.000) file.

    Returns an RDIData matching the MATLAB ``d`` struct from loadrdi.m.
    Velocities are in m/s; bad values are NaN. Time is Julian days.
    Reference: docs/legacy/loadrdi.m::rdread, rdflead, rdvlead, rdbtrack.
    """
    data = Path(path).read_bytes()
    ensembles = parse_pd0(data)
    if not ensembles:
        raise ValueError(f"No valid PD0 ensembles found in {path}")

    fl = ensembles[0]["fixed_leader"]
    nbin = fl["nbin"]
    nens = len(ensembles)

    u = np.full((nbin, nens), np.nan)
    v = np.full((nbin, nens), np.nan)
    w = np.full((nbin, nens), np.nan)
    e = np.full((nbin, nens), np.nan)
    heading = np.full(nens, np.nan)
    pitch = np.full(nens, np.nan)
    roll = np.full(nens, np.nan)
    time_julian = np.full(nens, np.nan)
    temp_c = np.full(nens, np.nan)
    sound_vel = np.full(nens, np.nan)
    echo = np.zeros((nbin, nens, 4), dtype=np.uint8)
    corr = np.zeros((nbin, nens, 4), dtype=np.uint8)
    pg = np.zeros((nbin, nens, 4), dtype=np.uint8)
    bt_range = np.full((4, nens), np.nan)
    bt_vel = np.full((4, nens), np.nan)

    for i, ens in enumerate(ensembles):
        vl = ens["variable_leader"]
        vel = ens["velocity"]
        n = min(nbin, vel.shape[0])

        u[:n, i] = vel[:n, 0]
        v[:n, i] = vel[:n, 1]
        w[:n, i] = vel[:n, 2]
        e[:n, i] = vel[:n, 3]

        heading[i] = vl["heading_deg"]
        pitch[i] = vl["pitch_deg"]
        roll[i] = vl["roll_deg"]
        time_julian[i] = vl["time_julian"]
        temp_c[i] = vl["temp_c"]
        sound_vel[i] = vl["sound_vel_ms"]

        if ens["echo"] is not None:
            echo[:n, i, :] = ens["echo"][:n, :]
        if ens["correlation"] is not None:
            corr[:n, i, :] = ens["correlation"][:n, :]
        if ens["percent_good"] is not None:
            pg[:n, i, :] = ens["percent_good"][:n, :]
        if ens["bottom_track"] is not None:
            bt_range[:, i] = ens["bottom_track"]["range_m"]
            bt_vel[:, i] = ens["bottom_track"]["vel_ms"]

    return RDIData(
        u=u,
        v=v,
        w=w,
        e=e,
        heading=heading,
        pitch=pitch,
        roll=roll,
        time_julian=time_julian,
        temp_c=temp_c,
        sound_vel_ms=sound_vel,
        echo=echo,
        corr=corr,
        pg=pg,
        btrack_range_m=bt_range,
        btrack_vel_ms=bt_vel,
        nbin=nbin,
        nens=nens,
        blen_m=fl["blen_m"],
        blnk_m=fl["blnk_m"],
        dist_m=fl["dist_m"],
        npng=fl["npng"],
        coord_transform=fl["coord_transform"],
        serial=fl["serial"],
        sysconfig=fl["sysconfig"],
        beams_up=fl["beams_up"],
        beam_angle_deg=fl["beam_angle_deg"],
        hdg_align_deg=fl["hdg_align_deg"],
        hdg_bias_deg=fl["hdg_bias_deg"],
    )
