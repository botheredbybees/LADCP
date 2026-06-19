# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

This repository is in a **pre-implementation phase**. No Python source code exists yet — the root contains only `docs/` and `.gitignore`. The next concrete step is scaffolding the Python package described in the proposal.

The planned stack (inferred from `.gitignore`): Python with `uv` or `pixi` for environment management, `ruff` for linting/formatting, and `pytest` for tests. Docker containerization is planned for deployment.

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

One cast from the I7N GO-SHIP cruise (2018 Indian Ocean), processed by A.M. Thurnherr (LDEO):

- `test_data/data/002.nc` — processed horizontal velocity NetCDF (LDEO_IX output). **This is the primary validation target.**
- `test_data/data/002.mat` / `002.txt` — same data in MATLAB and ASCII formats.
- `test_data/plots/` — 10 diagnostic PDFs produced by LDEO_IX (figures 01–14 with gaps).
- `test_data/ancillary/set_cast_params.m` — LDEO_IX processing parameters for this cruise. Shows the raw file naming convention: `{stn}/{stn}DL000.000` (downlooker) and `{stn}/{stn}UL000.000` (uplooker); CTD at `{stn}.1Hz` (30 header lines, 12 fields per line).
- `test_data/ancillary/ProcessingParams` — LADCP_w (Perl) configuration. CTD at `{stn}.6Hz` for vertical velocity processing.
- `test_data/cruise_data/` — compressed archives of the full cruise (raw PD0 files, SADCP `.mat` files, processed outputs).

The raw PD0 binary files and CTD ASCII files for cast 002 are in the cruise archives.

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
