# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

**In active development.** The `src/ladcp/` package is partially implemented with 142+ passing tests.

Completed layers: ingestion (PD0 binary, CTD SBE ASCII), coordinate transforms (beam→earth with gimbaled heading and bin-mapping), QA/editing (sidelobe masking, large-velocity rejection), shear solver, and inverse solver (GPS, SADCP, bottom-track, uplooker constraints).

Open gaps: SADCP NetCDF loader, SBE hex decoder, NetCDF output writer (`ladcp2cdf` equivalent), CLI wiring (stubs exist but raise `NotImplementedError`).

Stack: Python 3.11, `uv`, `ruff`, `pytest`, `numpy`/`xarray`/`scipy`/`netCDF4`. Docker image scaffolded.

## What this project is

A modern Python toolkit for **LADCP (Lowered Acoustic Doppler Current Profiler)** data processing — replacing a fragmented stack of MATLAB scripts, shell wrappers, and Perl-based acquisition utilities used in oceanographic research. The target users are physical oceanographers doing ship-based, full-water-column current profile work.

The scientific reference implementation to validate against is the **LDEO MATLAB workflow** (Thurnherr, Lamont-Doherty Earth Observatory), documented in `docs/legacy/`.

## Planned architecture

Five layers, in dependency order:

| Layer | Responsibility |
|---|---|
| **Ingestion** | Parse raw LADCP binary files (Teledyne RDI Workhorse format), CTD `.cnv` exports, GPS/nav data, optional shipboard ADCP (SADCP) |
| **Transforms** | Beam → instrument → ship → Earth coordinate transforms; heading, tilt, and rotation corrections; bin-mapping; sound-speed corrections |
| **Solution engine** | Shear-based and inverse/velocity-based solutions (matching the LDEO and JAMSTEC method families); comparison modes |
| **QA / diagnostics** | Tilt/heading plots, residual checks, bottom-track diagnostics, cast summary reports, machine-readable provenance |
| **Deployment** | Python API + CLI (`ladcp process`, `ladcp check`, …) + Docker image for repeatable execution |

## Key domain concepts

- **DL / UL**: Downlooker and Uplooker — the two co-mounted ADCPs on the CTD rosette.
- **Janus geometry**: Four slanted beams + optional vertical beam 5. The `janus5beam2earth()` function in `docs/legacy/ADCPtools/` shows the transform signature to replicate.
- **Bottom track**: Acoustic return from the sea floor used as an absolute velocity reference boundary condition.
- **SADCP**: Shipboard ADCP — another boundary condition for the inversion.
- **CTD time-series vs. CTD profile**: Raw un-binned CTD data (time-series) is used for depth/pressure; depth-binned CTD profiles are used for sound-speed correction. Both are separate inputs.
- **Timing**: Clock drift between DL and UL instruments and between the ADCP and CTD is a first-class processing concern. The legacy code tracks `params.timoff` and `params.timoff_uplooker` offsets explicitly.

## Legacy reference material (`docs/legacy/`)

These files are **read-only reference** — do not modify them.

There are two distinct legacy software systems, each handling a different output:

### LDEO_IX (MATLAB) — horizontal velocity (u, v)
The primary software used for GO-SHIP processing. Key files:
- `loadrdi.m` — reads RDI PD0 binary (`.000`) files into MATLAB structs `d`, `p`, `de`. The authoritative reference for the Python ingestion layer.
- `edit_data.m` — quality control, bin masking, large-velocity checks.
- `getshear2.m` / `prepinv.m` / `getinv.m` — shear calculation and inverse solver.
- `getbtrack.m` — bottom-track processing.
- `fixcompass.m` / `checktilt.m` — heading and tilt corrections.
- `loadctdprof.m` / `loadsadcp.m` — ancillary data loaders.
- `ladcp2cdf.m` — NetCDF output writer (defines the output schema to match).

### LADCP_w (Perl + ANTSlib) — vertical velocity (w) and VKE
- `docs/legacy/ANTSlib/` — Thurnherr's **ANTS** (Antilean Numerical Tool Suite), a Perl scientific computing framework. This is the computational backbone of the vertical velocity software.
- `docs/legacy/plot_mean_residuals.pl` — example of an ANTS pipeline script.
- The `ProcessingParams` ancillary file (see `test_data/ancillary/`) shows the Perl-based configuration syntax.

### Other legacy material
- `docs/legacy/ADCPtools/` — MATLAB coordinate transform library (apaloczy). The `janus5beam2earth` function signature and its `Gimbaled` / `Binmap` options are the authoritative reference for the Python transform layer.
- `docs/legacy/ladcp/` — IFM-GEOMAR/LDEO MATLAB LADCP Processing v10.16.2 (Krahmann et al.).
- `docs/legacy/LADCP_processing.md` — Real-world cruise report (P02E) describing acquisition commands.

## Test data (`test_data/`)

### 2018 S4P — primary validation dataset (`test_data/2018_S4P/`)

GO-SHIP cruise S4P (NBP1802), Southern Ocean, processed by A.M. Thurnherr (LDEO).
Full details: **`test_data/2018_S4P/DATA_SUMMARY.md`** — read that file first.

Key contents:
- `001.nc`, `002.nc`, `003.nc` — LDEO_IX processed outputs, **primary validation targets** for the inverse solver. Each file embeds GPS, CTD, SADCP, BT, and per-instrument DL/UL profiles alongside the final `u`/`v`, so they can drive solver tests without raw PD0 data.
- `processed_uv/` — processed NC outputs for 55+ casts (the full cruise).
- `CTD/320620180309_ctd.nc` — CCHDO calibrated CTD profiles (118 stations), useful for sound-speed correction.
- `CTD/00101.hex` etc. — raw SBE 24 Hz hex time-series with GPS; paired `.XMLCON` calibration files. Need SBE hex decoder (not yet implemented) to produce the `001.1Hz` ASCII files that LDEO_IX ingests.
- `SADCP/os75nb/contour/os75nb_short.nc` — OS75 SADCP NetCDF, 17 722 time steps × 60 depth cells, covers full cruise. **Longitude stored offset by −360° — normalize with `lon % 360`.**
- `set_cast_params.m` — LDEO_IX parameter file; documents raw file naming convention and final processing version (v8: DL+UL IMPed, GPS + SADCP + BT constraints).

Raw LADCP PD0 binary files (`001DL.101` etc.) are **not present**; they were in an archive that was not downloaded.

### I7N cruise data (`test_data/cruise_data/`, `test_data/data/`)

One cast from the I7N GO-SHIP cruise (2018 Indian Ocean), processed by A.M. Thurnherr (LDEO):

- `test_data/data/002.nc` — processed horizontal velocity NetCDF (LDEO_IX output).
- `test_data/cruise_data/` — compressed archives of the full cruise (raw PD0 files, SADCP `.mat` files, processed outputs).

### P16N cast003 — integration test data (`test_data/cruise_data/` or similar)

The primary dataset for end-to-end pipeline tests (ingestion → transforms → shear/inverse). Raw PD0 files and CTD data are available for this cast and are used by tests in `tests/integration/`.

## Validation-first principle

The proposal explicitly treats validation as a core deliverable, not an afterthought. When implementing:

1. Reproduce a known-good LDEO MATLAB output first before adding new features.
2. Tolerance bands for velocity profiles must be documented in tests.
3. Solutions should support running shear-based and inverse methods side-by-side for comparison.

## Development guidelines

- First phase: support one instrument family (Teledyne RDI Workhorse 300 kHz) and one cruise data convention before generalizing.
- Defer GUI development until the scientific core is stable.
- Keep coordinate transform math explicit and inspectable — avoid hiding rotation assumptions inside opaque functions.
- The Python package name is `ladcp` (implied by project name and domain conventions).
