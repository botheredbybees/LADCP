"""Janus beam → Earth coordinate transforms.

Primary reference: docs/legacy/loadrdi.m::b2earth (the LDEO_IX transform that
produced the validation reference outputs; Martini/Pluddeman algorithm).
Secondary: docs/legacy/ADCPtools/janus2earth.m, janus2xyz.m (Apaloczy et al.)
— same Earth rotation matrix, but its beam matrix is the UP-looking one and it
has no up/down switch; use ``beams_up`` here to select the correct matrix.
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
    beams_up: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert 4-beam Janus along-beam velocities to instrument frame (Vx, Vy, Vz).

    Faithful to loadrdi.m::b2earth Step 4 (lines 1736-1748), which is the
    transform that produced the LDEO_IX reference outputs:

        upward looking convex:    VX = (-b1+b2)  VY = (-b3+b4)  VZ = -(sum)
        downward looking convex:  VX = (+b1-b2)  VY = (-b3+b4)  VZ = +(sum)

    The two matrices differ by the physical 180-degree flip of the instrument
    about its y (beam-3) axis: x and z change sign, y does not.  Using the
    wrong matrix produces a MIRRORED horizontal velocity (a reflection about
    the instrument y-axis), which cannot be corrected by any heading rotation.

    NaN in any beam propagates to all outputs for that bin-ensemble.

    Parameters
    ----------
    b1, b2, b3, b4 : array_like, shape (nbin, nens)
        Along-beam velocities m/s. Positive = toward transducer face (TRDI convention).
    theta_deg : float
        Beam angle from vertical in degrees. 20.0 for RDI Workhorse 300/600 kHz.
    beams_up : bool
        False (default) for a down-looking instrument (DL); True for an
        up-looking instrument (UL).  Matches ``a.beams_up`` (the fixed-leader
        sysconfig up/down bit) in loadrdi.m.

    Returns
    -------
    Vx, Vy, Vz : ndarray, shape (nbin, nens)
        Instrument-frame velocity components m/s (z positive up in the
        Earth-facing sense used by the Step-5 rotation).
    """
    theta = np.radians(theta_deg)
    uvfac = 1.0 / (2.0 * np.sin(theta))
    wfac = 1.0 / (4.0 * np.cos(theta))
    if beams_up:
        Vx = uvfac * (-b1 + b2)
        Vy = uvfac * (-b3 + b4)
        Vz = wfac * (-b1 - b2 - b3 - b4)
    else:
        Vx = uvfac * (+b1 - b2)
        Vy = uvfac * (-b3 + b4)
        Vz = wfac * (+b1 + b2 + b3 + b4)
    return Vx, Vy, Vz


def reconstruct_3beam(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
    b4: np.ndarray,
) -> tuple[list[np.ndarray], int]:
    """Fill cells where exactly ONE beam is NaN by assuming zero error velocity.

    Port of loadrdi.m::b2earth lines 1713-1726 (p.allow_3beam_solutions,
    default on). The Janus error velocity is proportional to b1+b2-b3-b4;
    setting it to zero solves for the single missing beam:

        b1 = -b2 + b3 + b4        b3 = b1 + b2 - b4
        b2 = -b1 + b3 + b4        b4 = b1 + b2 - b3

    (identical for both beam orientations — the VE row of the matrix does
    not change sign with beams_up). Cells with 0 or >=2 missing beams are
    left untouched. Returns ([b1, b2, b3, b4] copies, n_3beam_cells).
    """
    beams = [np.asarray(b, dtype=np.float64).copy() for b in (b1, b2, b3, b4)]
    nan_mask = [np.isnan(b) for b in beams]
    n_missing = sum(m.astype(np.int8) for m in nan_mask)
    single = n_missing == 1
    n_3beam = int(single.sum())
    if n_3beam == 0:
        return beams, 0
    signs = {0: (-1, 1, 1), 1: (-1, 1, 1), 2: (1, 1, -1), 3: (1, 1, -1)}
    for i in range(4):
        fill = single & nan_mask[i]
        if not fill.any():
            continue
        others = [beams[j] for j in range(4) if j != i]
        s = signs[i]
        beams[i][fill] = (
            s[0] * others[0][fill] + s[1] * others[1][fill] + s[2] * others[2][fill]
        )
    return beams, n_3beam


def janus_error_velocity(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
    b4: np.ndarray,
    theta_deg: float,
) -> np.ndarray:
    """Janus error velocity VE = (b1 + b2 - b3 - b4) / (4 cos(theta)).

    loadrdi.m::b2earth Step 4 (lines 1741/1747; VES = VZS = 1/(4 C30)) --
    identical for both beam orientations. The redundancy residual of the
    4-beam solution: large |VE| flags inhomogeneous flow or noise-driven
    cells (in weak scattering, bins beyond the effective range return
    plausible-looking velocity from pure noise with |VE| ~ 0.7 m/s).
    NaN where any beam is missing -- 3-beam-reconstructed cells have
    VE = 0 by construction, so both are exempt from the elim edit,
    matching MATLAB.
    """
    return (
        np.asarray(b1, dtype=np.float64) + np.asarray(b2, dtype=np.float64)
        - np.asarray(b3, dtype=np.float64) - np.asarray(b4, dtype=np.float64)
    ) / (4.0 * np.cos(np.radians(theta_deg)))


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
    beams_up: bool = False,
    allow_3beam: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert 4-beam Janus velocities to Earth frame (u=East, v=North, w=Up).

    Rotation matrix per Appendix A of Dewey & Stringer (2007), identical to the
    Martini/Pluddeman matrix in loadrdi.m::b2earth Step 5.
    TRDI convention: positive beam velocity = toward transducer face.
    Heading: clockwise from y-axis (y=North when heading=0).

    Parameters
    ----------
    b1, b2, b3, b4 : ndarray, shape (nbin, nens)
        Along-beam velocities m/s.
    heading, pitch, roll : ndarray, shape (nens,)
        Instrument orientation angles in degrees, from the instrument's OWN
        sensors, unmodified.  An up-looking instrument needs no sign fiddling:
        its physical inversion is fully handled by ``beams_up``.
    theta_deg : float
        Beam angle from vertical in degrees (20.0° for Workhorse 300/600 kHz).
    gimbaled : bool
        If True, apply gimbaled heading correction (D&S 2007 eq. A2); if False,
        apply the fixed-sensor pitch correction (D&S eq. A1; Lohrmann et al.
        1990) — ``PP = asin(sin(p)cos(r)/sqrt(1-(sin(p)sin(r))^2))`` — which is
        what loadrdi.m::b2earth always applies ("fixed sensor case").  Use
        ``gimbaled=False`` to replicate LDEO_IX.
    beams_up : bool
        False (default) for a down-looking instrument; True for up-looking.
        Selects the beam-to-instrument matrix (see beam2xyz).  Matches the
        fixed-leader sysconfig up/down bit (``RDIData.beams_up``).
    allow_3beam : bool
        If True, cells with exactly one NaN beam are first reconstructed by
        assuming zero error velocity (see reconstruct_3beam; loadrdi.m's
        p.allow_3beam_solutions, default ON in LDEO_IX).  Default False here
        to preserve existing caller behavior — the P16N pipeline passes True.

    Returns
    -------
    u, v, w : ndarray, shape (nbin, nens)
        Earth-frame velocity components m/s: East, North, Up.
    """
    if allow_3beam:
        (b1, b2, b3, b4), _ = reconstruct_3beam(b1, b2, b3, b4)

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
    else:
        # Fixed-sensor pitch correction (loadrdi.m b2earth Step 1:
        # KA = sqrt(1 - (sin(p)sin(r))^2); PP = asin(sin(p)cos(r)/KA)).
        p = np.arcsin(Sph2 * Cph3 / np.sqrt(1.0 - (Sph2 * Sph3) ** 2))
        Sph2 = np.sin(p)
        Cph2 = np.cos(p)

    cx1 = Cph1 * Cph3 + Sph1 * Sph2 * Sph3
    cx2 = Sph1 * Cph3 - Cph1 * Sph2 * Sph3
    cx3 = Cph2 * Sph3
    cy1 = Sph1 * Cph2
    cy2 = Cph1 * Cph2
    cy3 = Sph2
    cz1 = Cph1 * Sph3 - Sph1 * Sph2 * Cph3
    cz2 = Sph1 * Sph3 + Cph1 * Sph2 * Cph3
    cz3 = Cph2 * Cph3

    Vx, Vy, Vz = beam2xyz(b1, b2, b3, b4, theta_deg, beams_up=beams_up)

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
        Rotation angle in degrees, positive = counterclockwise (East
        declination positive).

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
