# Oceanographer's Notes

Context for physical oceanographers evaluating or contributing to this toolkit.

## What This Is

A Python re-implementation of the LDEO LADCP processing workflow, targeting the same outputs as the MATLAB software used for GO-SHIP and similar hydrographic cruises. The goal is to reproduce known-good processed results numerically before adding any new capabilities.

The target validation dataset is I7N 2018, cast 002, processed by A.M. Thurnherr (LDEO). The reference output `test_data/data/002.nc` is the benchmark this software must match.

## What Is Implemented (v0.1.0)

**Ingestion only.** The parser reads Teledyne RDI Workhorse PD0 binary files (`.000` format) and returns a Python object containing:

- Velocity matrices (u, v, w, error) — (nbin × nens) arrays, m/s, NaN where instrument flagged bad
- Heading, pitch, roll, temperature, sound velocity — (nens,) arrays
- Ensemble timestamps — (nens,) Julian days, midnight-based (matching MATLAB `julian.m`)
- Echo intensity, correlation, percent-good — (nbin × nens × 4beam) arrays
- Bottom-track ranges and velocities — (4beam × nens)

Nothing downstream of ingestion exists yet: no coordinate transforms, no shear calculation, no velocity inversion, no CTD integration, no QA plots.

## MATLAB Equivalents

| MATLAB (LDEO_IX) | Python (this toolkit) | Status |
|---|---|---|
| `loadrdi.m` — read PD0 binary into structs `d`, `p`, `de` | `load_rdi(path)` → `RDIData` | **Done** |
| `janus5beam2earth()` (ADCPtools) — beam→Earth transform | `janus5beam2earth()` | Stub |
| `getshear2.m` — compute shear profiles | `compute_shear()` → `ShearProfile` | **Done** |
| `getinv.m` / `prepinv.m` — velocity inversion | (not yet designed) | Planned |
| `getbtrack.m` — bottom-track processing | (not yet designed) | Planned |
| `fixcompass.m` / `checktilt.m` — heading/tilt corrections | (not yet designed) | Planned |
| `loadctdprof.m` / `loadsadcp.m` — ancillary data | (not yet designed) | Planned |
| `ladcp2cdf.m` — NetCDF output | (not yet designed) | Planned |

The Perl-based LADCP_w software (vertical velocity, VKE) is a separate system not yet targeted.

## Key Scientific Conventions

### DL and UL

The downlooker (DL) and uplooker (UL) are two Workhorse ADCPs co-mounted on the CTD rosette, looking down and up respectively. Each has its own PD0 file (e.g. `003DL000.000` and `003UL000.000`). The `load_rdi()` function reads one file at a time — call it twice for a DL/UL pair.

### Timing

Clock drift between the DL and UL, and between either ADCP and the CTD, is a first-class processing concern. The legacy code tracks `params.timoff` and `params.timoff_uplooker` explicitly. Time in `RDIData` is stored as Julian days using the MATLAB midnight convention — not the astronomical noon-based Julian Day Number. This matters when computing DL/UL time overlaps.

### Coordinate frames

`RDIData.u` and `.v` are East and North *only if* the instrument was configured for Earth-frame output (RDI EX command `EX=11xxx`). GO-SHIP LADCP casts are normally configured this way, but the parser does not verify this — it reads whatever the instrument recorded. If you load a file from an instrument configured for beam-frame output, `u` contains beam 1 velocities, not East velocities.

The coordinate transform step (not yet implemented) will apply the Janus geometry, heading, pitch, and roll to convert from instrument frame to Earth frame for instruments not already in Earth-frame mode, and optionally apply tilt corrections (gimbaled) and bin-mapping.

### Velocity and bin geometry

The RDI Workhorse 300 kHz uses 4 Janus beams at 20° from vertical. Standard GO-SHIP configuration: 8 m bins, 25 bins, first bin centre ~8 m from the transducer face. `RDIData` carries `blen_m`, `nbin`, `dist_m` parsed directly from the fixed leader — no assumptions needed.

### Bottom track

Bottom-track data is in `RDIData.btrack_range_m` (4-beam slant range, m) and `btrack_vel_ms` (4-beam velocity, m/s). Not all ensembles have valid BT data — the instrument only returns BT when the acoustic return is strong enough. `NaN` in these arrays means no return in that ensemble. The BT velocities serve as an absolute velocity boundary condition in the inversion.

## Validation Status

The parser passes **sanity checks** against 2015 P16N cast 003:
- Ensemble count within expected range (~8970 for a ~3.6-hour cast)
- Bin geometry matches RDI 300 kHz Workhorse specification
- Headings in [0°, 360°]
- >50% of velocity cells are finite
- Time is monotonically increasing
- DL and UL time spans overlap (same cast, same wire)
- Some bottom-track data present

It has **not yet been numerically validated** against the reference output `test_data/data/002.nc`. That validation — checking that processed u/v profiles match to within documented tolerance bands — is the next milestone and requires completing the transform and solution layers first.

## The Two Reference Software Systems

This project targets two distinct MATLAB/Perl pipelines:

**LDEO_IX (MATLAB)** — processes horizontal velocity (u, v). This is the primary GO-SHIP processing tool. Source in `docs/legacy/LDEO_IX/` (if present) and referenced throughout `docs/legacy/loadrdi.m`, `getshear2.m`, etc.

**LADCP_w (Perl + ANTSlib)** — processes vertical velocity (w) and VKE (vertical kinetic energy). Uses Thurnherr's ANTS framework. Configuration via `ProcessingParams` files (see `test_data/ancillary/`). This system is documented in `docs/legacy/ANTSlib/` but is not yet targeted for Python re-implementation.

The two systems use different CTD time-series resolutions: LDEO_IX uses a 1 Hz CTD, LADCP_w uses 6 Hz. Both reference `test_data/ancillary/` for the processing parameter conventions.

## Roadmap

1. Complete ingestion validation (reproduce `002.nc` profiles to within tolerance)
2. Implement Janus beam→Earth coordinate transform (`janus5beam2earth`)
3. Implement shear-based solution (`getshear2` equivalent)
4. Implement inverse/velocity-based solution (`getinv` equivalent)
5. CTD and SADCP ancillary data loaders
6. QA diagnostics and plots (reproducing the 10 diagnostic PDFs in `test_data/plots/`)
7. End-to-end CLI (`ladcp process <cast>`)

## Contributing Scientific Validation

The most useful contribution right now is numerical validation: run the LDEO_IX MATLAB software on the I7N 2018 cast 002 data and document the intermediate arrays (shear profiles, bottom-track time series, final u/v profiles) so the Python implementation can be validated step by step. Raw PD0 files for cast 002 are not currently in the repository — contact the maintainer about data access.
