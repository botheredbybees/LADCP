"""Beam-to-Earth coordinate transforms. Reference: docs/legacy/ADCPtools/."""

from ladcp._typing import NDArray


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
