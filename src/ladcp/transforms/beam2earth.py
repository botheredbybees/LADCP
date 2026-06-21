"""Janus beam → Earth coordinate transforms.

Reference: docs/legacy/ADCPtools/janus2earth.m, janus2xyz.m (Apaloczy et al.).
TRDI convention: positive along-beam velocity = toward transducer face.
Heading increases clockwise from the y-axis (North when heading=0).
"""

import numpy as np

from ladcp._typing import NDArray


def beam2xyz(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
    b4: np.ndarray,
    theta_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert 4-beam Janus along-beam velocities to instrument frame (Vx, Vy, Vz).

    Vectorized numpy replacement for janus2xyz.m (which has a nested nz×nt loop).
    All inputs broadcast together; NaN in any beam propagates to all outputs for
    that bin-ensemble.

    Parameters
    ----------
    b1, b2, b3, b4 : array_like, shape (nbin, nens)
        Along-beam velocities m/s. Positive = toward transducer face (TRDI convention).
        Beam layout: 1=+x, 2=−x, 3=+y, 4=−y (beam angle theta from vertical).
    theta_deg : float
        Beam angle from vertical in degrees. 20.0° for RDI Workhorse 300/600 kHz.

    Returns
    -------
    Vx, Vy, Vz : ndarray, shape (nbin, nens)
        Instrument-frame velocity components m/s.
        x: beam-1 direction; y: beam-3 direction; z: upward.
    """
    theta = np.radians(theta_deg)
    uvfac = 1.0 / (2.0 * np.sin(theta))
    wfac = 1.0 / (4.0 * np.cos(theta))
    Vx = uvfac * (-b1 + b2)
    Vy = uvfac * (-b3 + b4)
    Vz = wfac * (-b1 - b2 - b3 - b4)
    return Vx, Vy, Vz


def beam2earth(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
    b4: np.ndarray,
    heading: np.ndarray,
    pitch: np.ndarray,
    roll: np.ndarray,
    theta_deg: float,
    gimbaled: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert 4-beam Janus velocities to Earth frame (u=East, v=North, w=Up).

    Implements Appendix A of Dewey & Stringer (2007) as coded in janus2earth.m.
    TRDI convention: positive beam velocity = toward transducer face.
    Heading: clockwise from y-axis (y=North when heading=0).

    Parameters
    ----------
    b1, b2, b3, b4 : ndarray, shape (nbin, nens)
        Along-beam velocities m/s.
    heading, pitch, roll : ndarray, shape (nens,)
        Instrument orientation angles in degrees.
    theta_deg : float
        Beam angle from vertical in degrees (20.0° for Workhorse 300/600 kHz).
    gimbaled : bool
        If True, apply gimbaled heading correction (D&S 2007 eq. A2); if False,
        apply fixed-mount pitch correction (eq. A1). Default True (LADCP standard).

    Returns
    -------
    u, v, w : ndarray, shape (nbin, nens)
        Earth-frame velocity components m/s: East, North, Up.
    """
    h = np.radians(heading)
    p = np.radians(pitch)
    r = np.radians(roll)

    Sph1 = np.sin(h)
    Cph1 = np.cos(h)
    Sph2 = np.sin(p)
    Cph2 = np.cos(p)
    Sph3 = np.sin(r)
    Cph3 = np.cos(r)

    if gimbaled:
        Sph2Sph3 = Sph2 * Sph3
        h = h + np.arcsin(Sph2Sph3 / np.sqrt(Cph2**2 + Sph2Sph3**2))
        Sph1 = np.sin(h)
        Cph1 = np.cos(h)

    cx1 = Cph1 * Cph3 + Sph1 * Sph2 * Sph3
    cx2 = Sph1 * Cph3 - Cph1 * Sph2 * Sph3
    cx3 = Cph2 * Sph3
    cy1 = Sph1 * Cph2
    cy2 = Cph1 * Cph2
    cy3 = Sph2
    cz1 = Cph1 * Sph3 - Sph1 * Sph2 * Cph3
    cz2 = Sph1 * Sph3 + Cph1 * Sph2 * Cph3
    cz3 = Cph2 * Cph3

    Vx, Vy, Vz = beam2xyz(b1, b2, b3, b4, theta_deg)

    u = +Vx * cx1 + Vy * cy1 + Vz * cz1
    v = -Vx * cx2 + Vy * cy2 - Vz * cz2
    w = -Vx * cx3 + Vy * cy3 + Vz * cz3

    return u, v, w


def uvrot(
    u: np.ndarray,
    v: np.ndarray,
    drot_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate velocity vectors counterclockwise by drot_deg degrees.

    Equivalent to MATLAB's uvrot(u, v, ang): applies R(drot) = [cos -sin; sin cos]
    to each (u, v) pair.  Used to convert magnetic-North-referenced velocities to
    true-North: pass drot_deg = +magnetic_declination_east.

    Parameters
    ----------
    u, v : ndarray
        East and North velocity components (any shape, must broadcast).
    drot_deg : float
        Rotation angle in degrees, positive = counterclockwise (East declination positive).

    Returns
    -------
    u_rot, v_rot : ndarray
        Rotated East and North components, same shape as inputs.
    """
    cr = np.cos(np.radians(drot_deg))
    sr = np.sin(np.radians(drot_deg))
    return u * cr - v * sr, u * sr + v * cr


def janus5beam2earth(
    heading: NDArray,
    pitch: NDArray,
    roll: NDArray,
    theta: float,
    b1: NDArray,
    b2: NDArray,
    b3: NDArray,
    b4: NDArray,
    *,
    gimbaled: bool = False,
    binmap: str = "none",
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Convert along-beam velocities to Earth coordinates.

    Returns ``(u, v, w, w5)`` — eastward, northward, upward,
    vertical-beam-only vertical velocity.

    ``gimbaled`` and ``binmap`` match the ``Gimbaled`` / ``Binmap`` kwargs
    in the ADCPtools MATLAB reference (docs/legacy/ADCPtools/README.md).
    """
    raise NotImplementedError
