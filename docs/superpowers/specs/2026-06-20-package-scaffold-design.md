# Design: Python Package Scaffolding for LADCP Toolkit

**Date:** 2026-06-20
**Status:** Approved

## Summary

Scaffold a working `ladcp` Python package using `uv`, `src/` layout, Python 3.11, `click` CLI, `ruff` + `pytest`, a two-stage Dockerfile, and a `Makefile`. The scaffold installs and runs (`ladcp --help`) before any algorithm is implemented. A real-data integration test layer is included but skips unless test data is present.

---

## File tree

```
ladcp/
├── src/
│   └── ladcp/
│       ├── __init__.py          # version + top-level exports
│       ├── cli.py               # click app, entry point
│       ├── _typing.py           # shared type aliases (NDArray, etc.)
│       ├── ingestion/
│       │   ├── __init__.py
│       │   └── rdi.py           # Teledyne RDI raw binary reader stub
│       ├── transforms/
│       │   ├── __init__.py
│       │   └── beam2earth.py    # janus5beam2earth stub (ref: docs/legacy/ADCPtools/)
│       ├── solution/
│       │   ├── __init__.py
│       │   └── shear.py         # shear-based solution stub (inverse solver added later)
│       └── qa/
│           ├── __init__.py
│           └── diagnostics.py   # diagnostic plot stubs
├── tests/
│   ├── conftest.py              # fixtures + integration skip logic
│   ├── test_smoke.py            # import + CLI --help (always runs)
│   └── integration/
│       └── test_cast.py         # real GO-SHIP cast tests (skipped if data absent)
├── test_data/                   # gitignored except sources.md
│   └── sources.md               # data acquisition instructions (already exists)
├── Dockerfile
├── .dockerignore
├── Makefile
└── pyproject.toml
```

No `deployment/` sub-package — Docker and the CLI entry point *are* the deployment layer.

---

## `pyproject.toml`

```toml
[project]
name = "ladcp"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "click>=8.1",
    "numpy>=1.26",
    "xarray>=2024.1",
    "netCDF4>=1.7",
    "scipy>=1.12",
    "matplotlib>=3.8",
]

[project.scripts]
ladcp = "ladcp.cli:app"

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/ladcp"]

[tool.uv]
dev-dependencies = ["pytest>=8", "pytest-cov", "ruff>=0.4"]

[tool.ruff.lint]
select = ["E", "F", "I", "NPY", "UP"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- `hatchling` is the standard uv-compatible build backend.
- `NPY` ruff rules catch NumPy-specific antipatterns.
- `UP` enforces modern Python syntax.

---

## CLI

Two stub subcommands matching the acquisition commands operators already know (`Lcheck` → `check`, full processing → `process`):

```python
# src/ladcp/cli.py
import click

@click.group()
@click.version_option()
def app():
    """LADCP processing toolkit."""

@app.command()
@click.argument("cast_file", type=click.Path(exists=True))
def process(cast_file):
    """Process a single LADCP cast file."""
    raise NotImplementedError("ingestion layer not yet implemented")

@app.command()
@click.argument("cast_file", type=click.Path(exists=True))
def check(cast_file):
    """Integrate vertical velocity to estimate zmax and zend."""
    raise NotImplementedError("ingestion layer not yet implemented")
```

`NotImplementedError` (not a silent no-op) makes it obvious a stub is a stub.

---

## Stub module conventions

Each stub module has a one-line docstring pointing to its MATLAB reference. Function signatures match the reference implementation so validation diffs are obvious.

```python
# src/ladcp/transforms/beam2earth.py
"""Beam-to-Earth coordinate transforms. Reference: docs/legacy/ADCPtools/."""

import numpy as np
from numpy.typing import NDArray

def janus5beam2earth(
    heading: NDArray, pitch: NDArray, roll: NDArray,
    theta: float,
    b1: NDArray, b2: NDArray, b3: NDArray, b4: NDArray,
    *,
    gimbaled: bool = False,
    binmap: str = "none",
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    # returns (u, v, w, w5) — eastward, northward, upward, vertical-beam-only
    raise NotImplementedError
```

The `gimbaled` and `binmap` kwargs match the ADCPtools `Gimbaled` and `Binmap` options exactly.

---

## Dockerfile

Two-stage build: builder installs deps via `uv`, runtime image is slim with no build tools.

```dockerfile
FROM python:3.11-slim AS builder
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
WORKDIR /app
COPY pyproject.toml .
RUN uv sync --no-install-project

FROM python:3.11-slim
COPY --from=builder /app/.venv /app/.venv
COPY src/ /app/src/
ENV PATH="/app/.venv/bin:$PATH"
ENTRYPOINT ["ladcp"]
```

---

## Makefile

```makefile
test:
	uv run pytest

test-integration:
	TEST_DATA_DIR=test_data uv run pytest -m integration

lint:
	uv run ruff check src tests
	uv run ruff format --check src tests

fmt:
	uv run ruff format src tests

docker-build:
	docker build -t ladcp .
```

---

## Real-data test layer

Test data is not committed. `test_data/` is gitignored (except `sources.md`). Integration tests are skipped unless `TEST_DATA_DIR` is set and populated (see `test_data/sources.md` for download instructions).

```python
# tests/conftest.py
import os
import pytest
from pathlib import Path

def pytest_configure(config):
    config.addinivalue_line("markers", "integration: requires real cast data")

@pytest.fixture
def test_data_dir():
    path = Path(os.environ.get("TEST_DATA_DIR", "test_data"))
    if not path.exists():
        pytest.skip("TEST_DATA_DIR not populated — see test_data/sources.md")
    return path
```

Recommended first dataset: one GO-SHIP cruise from NCEI GOSHIP-LADCP collection (raw + processed + CTD + SADCP). This gives a closed system for end-to-end validation against the LDEO MATLAB reference outputs.

---

## `.gitignore` additions

```
test_data/*
!test_data/sources.md
```
