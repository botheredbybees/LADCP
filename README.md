# LADCP Processing Toolkit

A modern Python toolkit for **LADCP (Lowered Acoustic Doppler Current Profiler)** data processing, replacing a fragmented stack of MATLAB scripts used in oceanographic research. The scientific reference is the LDEO MATLAB workflow (Thurnherr, Lamont-Doherty Earth Observatory).

## Status: Early Development

| Layer | Status |
|---|---|
| Ingestion — parse Teledyne RDI PD0 binary files | **Implemented** |
| Coordinate transforms (beam → Earth) | Planned |
| Shear-based and inverse velocity solutions | Planned |
| QA / diagnostics | Planned |
| End-to-end CLI (`ladcp process`) | Planned |

The `load_rdi()` function parses real PD0 files and returns validated numpy arrays. Everything downstream is a stub that raises `NotImplementedError`.

## Installation

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```bash
uv sync           # install all dependencies from uv.lock
uv pip install -e .   # install the ladcp package in editable mode
```

Or with Docker (see `Makefile` targets):

```bash
make build
make test
```

## Quick Start

```python
from pathlib import Path
from ladcp.ingestion import load_rdi

# Load a Teledyne RDI Workhorse PD0 file
d = load_rdi(Path("003DL000.000"))

print(f"Ensembles: {d.nens}, Bins: {d.nbin}")
print(f"Time span: {d.time_julian[0]:.4f} – {d.time_julian[-1]:.4f} (Julian day, midnight-based)")
print(f"Heading range: {d.heading.min():.1f}°–{d.heading.max():.1f}°")

import numpy as np
finite_frac = np.isfinite(d.u).mean()
print(f"Valid velocity cells: {finite_frac:.0%}")
```

## Running Tests

```bash
uv run pytest                          # unit tests only (no data files needed)
TEST_DATA_DIR=test_data uv run pytest  # + integration tests against real PD0 files
```

The integration tests require raw PD0 files. See `test_data/sources.md` for data provenance.

## Documentation

- [PROGRAMMERS_NOTES.md](PROGRAMMERS_NOTES.md) — architecture, design decisions, how to extend
- [OCEANOGRAPHERS_NOTES.md](OCEANOGRAPHERS_NOTES.md) — scientific context, MATLAB equivalents, validation strategy
- [USER_NOTES.md](USER_NOTES.md) — installation, input files, current limitations

## Project Layout

```
src/ladcp/
  ingestion/         # Implemented: PD0 parser, RDIData, load_rdi()
  transforms/        # Stub: beam → Earth coordinate transform
  solution/          # Stub: shear-based velocity inversion
  qa/                # Stub: diagnostics and plots
  cli.py             # Entry point: `ladcp process` / `ladcp check` (stubs)
docs/legacy/         # Read-only MATLAB reference (LDEO_IX, LADCP_w, ADCPtools)
test_data/           # Validation data (processed outputs + raw PD0 for one cruise)
```

## License

Not yet specified.
