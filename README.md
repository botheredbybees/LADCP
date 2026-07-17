# LADCP Processing Toolkit

A modern Python toolkit for **LADCP (Lowered Acoustic Doppler Current Profiler)** data processing, replacing a fragmented stack of MATLAB scripts used in oceanographic research. The scientific reference is the LDEO MATLAB workflow (Thurnherr, Lamont-Doherty Earth Observatory).

## Status: Active Development

| Layer | Status |
|---|---|
| Ingestion — RDI PD0 binary, SBE CNV (ASCII & binary), SBE hex time-series, SADCP NetCDF, UH/CLIVAR CTD time-series | **Implemented** |
| Coordinate transforms — beam → Earth (gimbaled), magnetic declination rotation, sound-speed correction | **Implemented** |
| QA / editing — sidelobe masking, large/error-velocity and vertical-velocity outlier removal, PPI editing | **Implemented** |
| Shear-based solution | **Implemented** |
| Inverse velocity solution — constrained least-squares with GPS, SADCP, bottom-track | **Implemented** |
| NetCDF output writer | **Implemented** |
| End-to-end CLI (`ladcp process`) | Planned |

The inverse solver produces velocity profiles for full-water-column casts. Validation against the
LDEO MATLAB reference on **P16N 2015 cast 003 meets both targets** (u RMSE 0.045 m/s, v RMSE
0.033 m/s, both under the 0.05 m/s GO-SHIP tolerance — hard test assertions, not `xfail`).

Multi-cast validation across three unseen cruises (263 tests total, 255 passed / 8 skipped):

| Cruise | Casts | Pass both (u, v < 0.05) | u RMSE median | Notes |
|---|---|---|---|---|
| P16N 2015 | 1 (tuning cast) | — | 0.045 | Reference cast the pipeline was built against |
| I7N 2018 | 124/124 | 53 (43%) | 0.043 | 10 casts numerically "explode" (RMSE ~10⁶–10¹⁰) — an open ill-conditioning lead |
| A16N 2013 | 95/95 | 15 (16%) | 0.34 | Deep casts (>4 km) fail with large mid-column swings — open investigation, several leads ruled out |

See `docs/HANDOVER.md` for the current session-to-session status and
`docs/validation/BULK_VALIDATION_REPORT.md` for full per-cast numbers.

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
uv run python -m pytest                          # unit tests only (no data files needed)
TEST_DATA_DIR=test_data uv run python -m pytest  # + integration tests against real PD0 files
```

**Note:** `uv run pytest` (without `python -m`) fails on some machines with a broken
entry-point shim — use the `python -m pytest` form above.

Integration tests require raw PD0, CTD, and reference NetCDF files for P16N, I7N, and
A16N. See `test_data/sources.md` for data provenance.

## Documentation

- [PROGRAMMERS_NOTES.md](PROGRAMMERS_NOTES.md) — architecture, design decisions, how to extend
- [OCEANOGRAPHERS_NOTES.md](OCEANOGRAPHERS_NOTES.md) — scientific context, MATLAB equivalents, validation status
- [USER_NOTES.md](USER_NOTES.md) — installation, input files, current limitations
- [docs/HANDOVER.md](docs/HANDOVER.md) — current session-to-session status and next steps
- [docs/validation/](docs/validation/) — bulk multi-cruise validation brief, report, and the RMSE-closure plan

## Project Layout

```
src/ladcp/
  ingestion/         # RDI PD0, SBE CNV/hex, SADCP NetCDF, UH/CLIVAR CTD time-series loaders
  transforms/        # beam2earth (gimbaled Janus), uvrot (magnetic declination), sound-speed
  solution/          # Super-ensemble formation, shear solution, constrained inverse solver
  qa/                # Sidelobe / velocity / vertical-velocity / PPI editing, diagnostics
  output/            # NetCDF output writer
  cli.py             # Entry point stubs: `ladcp process`, `ladcp check` (not yet wired)
docs/legacy/         # Read-only MATLAB reference (LDEO_IX, LADCP_w, ADCPtools)
docs/superpowers/    # Design specs and implementation plans
docs/validation/     # Bulk validation brief/report, RMSE-closure plan
docs/history/        # Superseded point-in-time investigation plans
docs/HANDOVER.md     # Current session-to-session status
test_data/           # P16N, I7N, A16N (+ S4P) cruise raw files and LDEO reference output
```

## License

MIT
