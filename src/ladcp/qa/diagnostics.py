"""Diagnostic plots and QC summaries. Reference: docs/legacy/plotraw.m, plotinv.m."""

from pathlib import Path


def tilt_heading_plot(data: dict, output_path: Path) -> None:
    """Plot tilt and heading time series for a cast.

    Reproduces figures 01–02 in test_data/plots/.
    Reference: docs/legacy/plotraw.m.
    """
    raise NotImplementedError


def residual_plot(data: dict, output_path: Path) -> None:
    """Plot velocity residuals after inversion.

    Reproduces figures 09–11 in test_data/plots/.
    Reference: docs/legacy/plotinv.m.
    """
    raise NotImplementedError
