"""Shear-based horizontal velocity solution. Reference: docs/legacy/getshear2.m."""

from ladcp._typing import NDArray


def shear_solution(
    u_shear: NDArray,
    v_shear: NDArray,
    depth: NDArray,
) -> tuple[NDArray, NDArray]:
    """Integrate velocity shear profiles to absolute velocities.

    Returns ``(u, v)`` — eastward and northward velocity profiles.
    Reference: docs/legacy/getshear2.m and docs/legacy/getinv.m.
    """
    raise NotImplementedError
