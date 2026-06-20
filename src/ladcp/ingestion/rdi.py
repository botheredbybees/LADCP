"""Read Teledyne RDI PD0 binary files. Reference: docs/legacy/loadrdi.m."""

from pathlib import Path

import numpy as np

from ladcp.ingestion._pd0 import parse_pd0
from ladcp.ingestion._types import RDIData


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
        serial=fl["serial"],
    )
