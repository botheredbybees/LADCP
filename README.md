# LADCP Processing Toolkit

A modern Python toolkit for **LADCP (Lowered Acoustic Doppler Current Profiler)** data processing, replacing a fragmented stack of MATLAB scripts used in oceanographic research. The scientific reference is the LDEO MATLAB workflow (Thurnherr, Lamont-Doherty Earth Observatory).

## Status: Active Development

| Layer | Status |
|---|---|
| Ingestion — parse Teledyne RDI PD0 binary; load SBE CNV (ASCII & binary) | **Implemented** |
| Coordinate transforms — beam → Earth (gimbaled), magnetic declination rotation | **Implemented** |
| QA / editing — sidelobe masking, large-velocity and vertical-velocity outlier removal | **Implemented** |
| Shear-based solution | **Implemented** |
| Inverse velocity solution — constrained least-squares with GPS, SADCP, bottom-track | **Implemented** |
| End-to-end CLI (`ladcp process`) | Planned |

The inverse solver produces velocity profiles for full-water-column casts. Validation against
LDEO MATLAB reference output (P16N 2015, cast 003) gives u RMSE ≈ 0.07 m/s; the 0–1000 m
range correlates at r ≈ +0.90. A systematic anti-correlation at 1000–2000 m is under investigation
and requires MATLAB intermediate arrays to resolve (see `docs/superpowers/plans/2026-06-22-rmse-closure.md`).

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
from ladcp.ingestion.rdi import load_rdi
from ladcp.ingestion.ctd import load_ctd, assign_bin_depths
from ladcp.transforms.beam2earth import beam2earth, uvrot
from ladcp.qa.editing import edit_sidelobes, edit_large_velocities, edit_w_outliers
from ladcp.solution.inverse import EnsembleData, prepare_superensembles, compute_inverse

rdi = load_rdi(Path("003DL000.000"))
ctd = load_ctd(Path("003_01.cnv"))

u, v, w = beam2earth(rdi.u, rdi.v, rdi.w, rdi.e,
                     rdi.heading, rdi.pitch, rdi.roll,
                     theta_deg=20.0, gimbaled=True)
u, v = uvrot(u, v, angle_deg=-12.3)  # magnetic declination correction

z, izm = assign_bin_depths(rdi, ctd, looker="down")

ens = EnsembleData(u=u, v=v, w=w, ...)
ens = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=rdi.blen_m)
ens = edit_large_velocities(ens)

se = prepare_superensembles(ens)
result = compute_inverse(se, u_ship=0.002, v_ship=-0.001)

print(f"Profile depth range: {result.z.min():.0f}–{result.z.max():.0f} m")
```

## Running Tests

```bash
uv run pytest                          # unit tests only (no data files needed)
TEST_DATA_DIR=test_data uv run pytest  # + integration tests against real PD0 files
```

Integration tests for P16N cast 003 require raw PD0, CTD, and reference NetCDF files.
See `test_data/sources.md` for data provenance.

## Documentation

- [PROGRAMMERS_NOTES.md](PROGRAMMERS_NOTES.md) — architecture, design decisions, how to extend
- [OCEANOGRAPHERS_NOTES.md](OCEANOGRAPHERS_NOTES.md) — scientific context, MATLAB equivalents, validation status
- [USER_NOTES.md](USER_NOTES.md) — installation, input files, current limitations

## Project Layout

```
src/ladcp/
  ingestion/         # RDI PD0 parser, SBE CNV loader, bin-depth assignment
  transforms/        # beam2earth (gimbaled Janus), uvrot (magnetic declination)
  solution/          # Super-ensemble formation + constrained inverse solver
  qa/                # Sidelobe / velocity / vertical-velocity editing
  cli.py             # Entry point stubs: `ladcp process`, `ladcp check`
docs/legacy/         # Read-only MATLAB reference (LDEO_IX, LADCP_w, ADCPtools)
docs/superpowers/    # Design specs and implementation plans
test_data/           # P16N 2015 cast 003 raw files + LDEO reference output
```

## License

MIT
