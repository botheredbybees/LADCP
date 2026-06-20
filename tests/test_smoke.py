"""Smoke tests: package imports and CLI entry point."""

import pytest
from click.testing import CliRunner


def test_version():
    import ladcp

    assert ladcp.__version__ == "0.1.0"


def test_cli_help():
    from ladcp.cli import app

    result = CliRunner().invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "LADCP processing toolkit" in result.output


def test_cli_process_help():
    from ladcp.cli import app

    result = CliRunner().invoke(app, ["process", "--help"])
    assert result.exit_code == 0
    assert "cast_file" in result.output.lower()


def test_cli_check_help():
    from ladcp.cli import app

    result = CliRunner().invoke(app, ["check", "--help"])
    assert result.exit_code == 0
    assert "cast_file" in result.output.lower()


def test_stubs_importable():
    """All stub modules import without error."""
    from ladcp.ingestion import rdi  # noqa: F401
    from ladcp.qa import diagnostics  # noqa: F401
    from ladcp.solution import shear  # noqa: F401
    from ladcp.transforms import beam2earth  # noqa: F401


def test_stubs_raise_not_implemented():
    """Stubs raise NotImplementedError, not silently pass."""
    from pathlib import Path

    import numpy as np

    from ladcp.qa.diagnostics import tilt_heading_plot
    from ladcp.solution.shear import shear_solution
    from ladcp.transforms.beam2earth import janus5beam2earth

    dummy = np.zeros((10, 8))
    dummy_1d = np.zeros(10)

    with pytest.raises(NotImplementedError):
        janus5beam2earth(dummy_1d, dummy_1d, dummy_1d, 20.0, dummy, dummy, dummy, dummy)

    with pytest.raises(NotImplementedError):
        shear_solution(dummy, dummy, dummy_1d)

    with pytest.raises(NotImplementedError):
        tilt_heading_plot({}, Path("out.pdf"))
