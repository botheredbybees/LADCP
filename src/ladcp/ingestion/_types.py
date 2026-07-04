"""RDIData dataclass — result of load_rdi().

Reference: MATLAB 'd' struct in loadrdi.m.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass
class RDIData:
    """Parsed RDI PD0 data for one instrument file (downlooker or uplooker).

    Array shapes use axis convention (bins, ensembles) to match MATLAB loadrdi.m.

    Coordinate assumption: u/v/w/e are labeled as Earth-frame velocities.
    This is valid only when the PD0 file was recorded in Earth coordinates
    (EX command = 11xxx). Beam or instrument coordinate files will have their
    raw beam velocities stored here without transformation.
    """

    u: NDArray[np.float64]  # eastward velocity  (nbin, nens) m/s
    v: NDArray[np.float64]  # northward velocity (nbin, nens) m/s
    w: NDArray[np.float64]  # vertical velocity  (nbin, nens) m/s
    e: NDArray[np.float64]  # error velocity     (nbin, nens) m/s
    heading: NDArray[np.float64]  # degrees            (nens,)
    pitch: NDArray[np.float64]  # degrees            (nens,)
    roll: NDArray[np.float64]  # degrees            (nens,)
    time_julian: NDArray[np.float64]  # Julian days        (nens,)
    temp_c: NDArray[np.float64]  # Celsius            (nens,)
    sound_vel_ms: NDArray[np.float64]  # m/s                (nens,)
    echo: NDArray[np.uint8]  # amplitude          (nbin, nens, 4)
    corr: NDArray[np.uint8]  # correlation        (nbin, nens, 4)
    pg: NDArray[np.uint8]  # percent good       (nbin, nens, 4)
    btrack_range_m: NDArray[np.float64]  # (4, nens)
    btrack_vel_ms: NDArray[np.float64]  # (4, nens)
    nbin: int
    nens: int
    blen_m: float
    blnk_m: float
    dist_m: float
    npng: int
    coord_transform: int  # EX byte: bits 4-3 = frame, 2=tilt, 1=binmap, 0=3beam
    serial: list[int]
    sysconfig: int = 0  # fixed-leader system-configuration word (raw uint16)
    beams_up: bool = False  # sysconfig LSB bit 7: True = up-facing (uplooker)
    beam_angle_deg: float = 20.0  # sysconfig MSB bits 0-1: 15/20/30 deg
    hdg_align_deg: float = 0.0  # EA heading-alignment word, degrees
    hdg_bias_deg: float = 0.0  # EB heading-bias word, degrees
