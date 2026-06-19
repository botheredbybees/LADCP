"""Shared fixtures and integration test configuration."""

import os
from pathlib import Path

import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "integration: requires real cast data in TEST_DATA_DIR"
    )


@pytest.fixture
def test_data_dir() -> Path:
    """Path to the test data directory. Skips test if not populated."""
    path = Path(os.environ.get("TEST_DATA_DIR", "test_data"))
    if not path.exists():
        pytest.skip("TEST_DATA_DIR not populated — see test_data/sources.md")
    return path
