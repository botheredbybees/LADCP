# Package Scaffold Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a working `ladcp` Python package that installs, passes linting, runs `ladcp --help`, and has a real-data integration test that skips cleanly when test data is absent.

**Architecture:** `src/` layout with `hatchling` build backend and `uv` for environment management. Five sub-packages (`ingestion`, `transforms`, `solution`, `qa`) with `NotImplementedError` stubs and docstrings pointing to MATLAB/Perl reference files. A `click` CLI, two-stage Dockerfile, and `Makefile`. Integration tests gated by `TEST_DATA_DIR` env var.

**Tech Stack:** Python 3.11, uv, hatchling, click, numpy, xarray, netCDF4, scipy, matplotlib, pytest, ruff, Docker

---

### Task 1: Update `.gitignore` for test data

**Files:**
- Modify: `.gitignore`

- [ ] **Step 1: Add test_data exclusions to `.gitignore`**

Open `.gitignore` and add these lines at the top, after the existing `docs/legacy/` entry:

```
# test data (large binary files — see test_data/sources.md for acquisition)
test_data/*
!test_data/sources.md
```

- [ ] **Step 2: Verify test_data is now ignored**

```powershell
git status --short
```

Expected: `test_data/` no longer appears as `??` (untracked). Only `test_data/sources.md` should show as untracked if not yet committed.

- [ ] **Step 3: Commit**

```powershell
git add .gitignore
git commit -m "chore: ignore test_data binaries, keep sources.md"
```

---

### Task 2: Create `pyproject.toml` and initialise `uv`

**Files:**
- Create: `pyproject.toml`

- [ ] **Step 1: Create `pyproject.toml`**

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

- [ ] **Step 2: Run `uv sync` to create lockfile and venv**

```powershell
uv sync
```

Expected: resolves deps, creates `.venv/` and `uv.lock`. No errors.

- [ ] **Step 3: Commit**

```powershell
git add pyproject.toml uv.lock
git commit -m "chore: add pyproject.toml and uv lockfile"
```

---

### Task 3: Package skeleton and CLI — TDD

**Files:**
- Create: `tests/test_smoke.py`
- Create: `src/ladcp/__init__.py`
- Create: `src/ladcp/_typing.py`
- Create: `src/ladcp/cli.py`
- Create: `src/ladcp/ingestion/__init__.py`
- Create: `src/ladcp/transforms/__init__.py`
- Create: `src/ladcp/solution/__init__.py`
- Create: `src/ladcp/qa/__init__.py`

- [ ] **Step 1: Write the failing smoke tests**

Create `tests/test_smoke.py`:

```python
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
```

- [ ] **Step 2: Run to verify tests fail**

```powershell
uv run pytest tests/test_smoke.py -v
```

Expected: 4 errors — `ModuleNotFoundError: No module named 'ladcp'`

- [ ] **Step 3: Create `src/ladcp/__init__.py`**

```python
"""LADCP processing toolkit."""

__version__ = "0.1.0"
```

- [ ] **Step 4: Create `src/ladcp/_typing.py`**

```python
"""Shared type aliases used across sub-packages."""

from numpy.typing import NDArray

__all__ = ["NDArray"]
```

- [ ] **Step 5: Create `src/ladcp/cli.py`**

```python
"""CLI entry point."""

import click


@click.group()
@click.version_option()
def app() -> None:
    """LADCP processing toolkit."""


@app.command()
@click.argument("cast_file", type=click.Path(exists=True))
def process(cast_file: str) -> None:
    """Process a single LADCP cast file."""
    raise NotImplementedError("ingestion layer not yet implemented")


@app.command()
@click.argument("cast_file", type=click.Path(exists=True))
def check(cast_file: str) -> None:
    """Integrate vertical velocity to estimate zmax and zend."""
    raise NotImplementedError("ingestion layer not yet implemented")
```

- [ ] **Step 6: Create empty `__init__.py` for each sub-package**

Create these four files, each with the single line `""""""` (empty docstring):

- `src/ladcp/ingestion/__init__.py`
- `src/ladcp/transforms/__init__.py`
- `src/ladcp/solution/__init__.py`
- `src/ladcp/qa/__init__.py`

Each file:
```python
""""""
```

- [ ] **Step 7: Re-sync so uv picks up the new source tree**

```powershell
uv sync
```

Expected: installs `ladcp` in editable mode, no errors.

- [ ] **Step 8: Run tests to verify they pass**

```powershell
uv run pytest tests/test_smoke.py -v
```

Expected:
```
PASSED tests/test_smoke.py::test_version
PASSED tests/test_smoke.py::test_cli_help
PASSED tests/test_smoke.py::test_cli_process_help
PASSED tests/test_smoke.py::test_cli_check_help
```

- [ ] **Step 9: Verify CLI works directly**

```powershell
uv run ladcp --help
```

Expected output contains `LADCP processing toolkit` and lists `check` and `process` subcommands.

- [ ] **Step 10: Commit**

```powershell
git add src/ tests/test_smoke.py
git commit -m "feat: add package skeleton, CLI entry point, and smoke tests"
```

---

### Task 4: Stub modules — TDD

**Files:**
- Modify: `tests/test_smoke.py`
- Create: `src/ladcp/ingestion/rdi.py`
- Create: `src/ladcp/transforms/beam2earth.py`
- Create: `src/ladcp/solution/shear.py`
- Create: `src/ladcp/qa/diagnostics.py`

- [ ] **Step 1: Add stub tests to `tests/test_smoke.py`**

Append to the existing file:

```python

def test_stubs_importable():
    """All stub modules import without error."""
    from ladcp.ingestion import rdi  # noqa: F401
    from ladcp.qa import diagnostics  # noqa: F401
    from ladcp.solution import shear  # noqa: F401
    from ladcp.transforms import beam2earth  # noqa: F401


def test_stubs_raise_not_implemented():
    """Stubs raise NotImplementedError, not silently pass."""
    from pathlib import Path
    from ladcp.ingestion.rdi import load_rdi
    from ladcp.transforms.beam2earth import janus5beam2earth
    from ladcp.solution.shear import shear_solution
    from ladcp.qa.diagnostics import tilt_heading_plot

    import numpy as np

    dummy = np.zeros((10, 8))
    dummy_1d = np.zeros(10)

    with pytest.raises(NotImplementedError):
        load_rdi(Path("nonexistent.000"))

    with pytest.raises(NotImplementedError):
        janus5beam2earth(dummy_1d, dummy_1d, dummy_1d, 20.0,
                         dummy, dummy, dummy, dummy)

    with pytest.raises(NotImplementedError):
        shear_solution(dummy, dummy, dummy_1d)

    with pytest.raises(NotImplementedError):
        tilt_heading_plot({}, Path("out.pdf"))
```

- [ ] **Step 2: Run to verify new tests fail**

```powershell
uv run pytest tests/test_smoke.py::test_stubs_importable tests/test_smoke.py::test_stubs_raise_not_implemented -v
```

Expected: `ModuleNotFoundError` for each stub module.

- [ ] **Step 3: Create `src/ladcp/ingestion/rdi.py`**

```python
"""Read Teledyne RDI PD0 binary files. Reference: docs/legacy/loadrdi.m."""

from pathlib import Path


def load_rdi(path: Path) -> dict:
    """Load one RDI PD0 binary (.000) file.

    Returns a dict matching the MATLAB ``d`` struct from loadrdi.m:
    keys will include ``vel`` (velocity, ensembles × bins × beams),
    ``heading``, ``pitch``, ``roll``, ``time``, ``pressure``.
    """
    raise NotImplementedError
```

- [ ] **Step 4: Create `src/ladcp/transforms/beam2earth.py`**

```python
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
```

- [ ] **Step 5: Create `src/ladcp/solution/shear.py`**

```python
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
```

- [ ] **Step 6: Create `src/ladcp/qa/diagnostics.py`**

```python
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
```

- [ ] **Step 7: Run all tests**

```powershell
uv run pytest tests/test_smoke.py -v
```

Expected: all 6 tests pass.

- [ ] **Step 8: Run linter**

```powershell
uv run ruff check src tests
uv run ruff format --check src tests
```

Expected: no errors. If format errors appear, run `uv run ruff format src tests` then re-check.

- [ ] **Step 9: Commit**

```powershell
git add src/ladcp/ tests/test_smoke.py
git commit -m "feat: add stub modules with NotImplementedError and reference docstrings"
```

---

### Task 5: Integration test layer

**Files:**
- Create: `tests/conftest.py`
- Create: `tests/integration/__init__.py`
- Create: `tests/integration/test_cast.py`

- [ ] **Step 1: Create `tests/conftest.py`**

```python
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
```

- [ ] **Step 2: Create `tests/integration/__init__.py`**

Empty file — required for pytest to discover the sub-directory:

```python
""""""
```

- [ ] **Step 3: Create `tests/integration/test_cast.py`**

```python
"""Integration tests against I7N GO-SHIP cast 002 reference data."""

import pytest
import xarray as xr


@pytest.mark.integration
def test_reference_netcdf_loads(test_data_dir):
    """Reference LDEO_IX output for cast 002 is readable and has data."""
    nc_path = test_data_dir / "data" / "002.nc"
    assert nc_path.exists(), f"Reference NetCDF not found at {nc_path}"

    ds = xr.open_dataset(nc_path)
    assert len(ds.data_vars) > 0, "NetCDF has no data variables"
    ds.close()


@pytest.mark.integration
def test_reference_netcdf_has_velocity(test_data_dir):
    """Reference NetCDF contains at least one velocity variable."""
    ds = xr.open_dataset(test_data_dir / "data" / "002.nc")
    vel_vars = [v for v in ds.data_vars if any(
        tok in v.lower() for tok in ("u", "v", "vel", "east", "north")
    )]
    assert vel_vars, (
        f"No velocity variable found. Variables present: {list(ds.data_vars)}"
    )
    ds.close()


@pytest.mark.integration
def test_reference_cast_depth(test_data_dir):
    """Cast 002 reaches approximately 4892 m (within 50 m tolerance)."""
    ds = xr.open_dataset(test_data_dir / "data" / "002.nc")
    # depth coordinate may be named 'depth', 'z', 'pressure', or similar
    depth_coords = [c for c in ds.coords if any(
        tok in c.lower() for tok in ("depth", "z", "pressure", "p")
    )]
    assert depth_coords, f"No depth coordinate found. Coords: {list(ds.coords)}"
    import numpy as np
    max_depth = float(np.nanmax(ds[depth_coords[0]].values))
    assert abs(max_depth - 4892) < 50, (
        f"Expected max depth ~4892 m, got {max_depth:.1f} m"
    )
    ds.close()
```

- [ ] **Step 4: Run unit tests to confirm they still pass**

```powershell
uv run pytest tests/test_smoke.py -v
```

Expected: all 6 pass.

- [ ] **Step 5: Run integration tests without data — verify skip**

```powershell
uv run pytest tests/integration/ -v -m integration
```

Expected:
```
SKIPPED tests/integration/test_cast.py::test_reference_netcdf_loads
SKIPPED tests/integration/test_cast.py::test_reference_netcdf_has_velocity
SKIPPED tests/integration/test_cast.py::test_reference_cast_depth
```
(All skip with message "TEST_DATA_DIR not populated")

- [ ] **Step 6: Run integration tests with data — verify they pass**

```powershell
$env:TEST_DATA_DIR = "test_data"; uv run pytest tests/integration/ -v -m integration
```

Expected: all 3 integration tests pass. If a depth test fails, inspect the coordinate names:

```powershell
$env:TEST_DATA_DIR = "test_data"; uv run python -c "import xarray as xr; ds = xr.open_dataset('test_data/data/002.nc'); print(list(ds.coords)); print(list(ds.data_vars))"
```

Adjust the coordinate name detection strings in `test_cast.py` to match what's actually in the file.

- [ ] **Step 7: Commit**

```powershell
git add tests/conftest.py tests/integration/
git commit -m "test: add integration test layer with TEST_DATA_DIR skip guard"
```

---

### Task 6: Dockerfile, `.dockerignore`, and `Makefile`

**Files:**
- Create: `Dockerfile`
- Create: `.dockerignore`
- Create: `Makefile`

- [ ] **Step 1: Create `Dockerfile`**

```dockerfile
FROM python:3.11-slim AS builder
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
WORKDIR /app
COPY pyproject.toml .
COPY src/ src/
RUN uv venv && uv pip install .

FROM python:3.11-slim
COPY --from=builder /app/.venv /app/.venv
ENV PATH="/app/.venv/bin:$PATH"
ENTRYPOINT ["ladcp"]
```

- [ ] **Step 2: Create `.dockerignore`**

```
.venv/
__pycache__/
*.pyc
.pytest_cache/
.ruff_cache/
test_data/
docs/
*.md
Makefile
uv.lock
```

- [ ] **Step 3: Create `Makefile`**

> **Important:** Makefile recipes must use a hard tab character (not spaces) for indentation. If your editor auto-converts tabs, set it to insert literal tabs for `.mk` / `Makefile` files.

```makefile
.PHONY: test test-integration lint fmt docker-build

test:
	uv run pytest

test-integration:
	TEST_DATA_DIR=test_data uv run pytest -m integration -v

lint:
	uv run ruff check src tests
	uv run ruff format --check src tests

fmt:
	uv run ruff format src tests

docker-build:
	docker build -t ladcp .
```

- [ ] **Step 4: Run `make test` to confirm the full suite passes**

```powershell
uv run pytest
```

(Use this directly on Windows since `make` may not be installed. `make test` works in Docker/Linux.)

Expected: 6 tests pass, 0 failures.

- [ ] **Step 5: Run `make lint` equivalent**

```powershell
uv run ruff check src tests; uv run ruff format --check src tests
```

Expected: no errors.

- [ ] **Step 6: Build Docker image**

```powershell
docker build -t ladcp .
```

Expected: image builds successfully. Final layer is `python:3.11-slim`.

- [ ] **Step 7: Verify Docker CLI works**

```powershell
docker run --rm ladcp --help
```

Expected: prints `LADCP processing toolkit` and lists `check` and `process` subcommands.

- [ ] **Step 8: Commit**

```powershell
git add Dockerfile .dockerignore Makefile
git commit -m "chore: add Dockerfile, .dockerignore, and Makefile"
```
