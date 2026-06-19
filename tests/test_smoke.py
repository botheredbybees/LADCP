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
