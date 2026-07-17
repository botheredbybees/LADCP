# Oceanographer's Notes

Context for physical oceanographers evaluating or contributing to this toolkit.

## What This Is

A Python re-implementation of the LDEO LADCP processing workflow, targeting the same outputs as the MATLAB software used for GO-SHIP and similar hydrographic cruises. The goal is to reproduce known-good processed results numerically before adding any new capabilities.

The primary tuning dataset was **2015 P16N, cast 003**, processed by A.M. Thurnherr (LDEO); the reference output `test_data/2015_P16N/003.nc` is the cast the pipeline's design decisions were validated against. To guard against overfitting to that one cast, the pipeline has since been run against two entire cruises it never saw during development — **I7N 2018** (124 casts) and **A16N 2013** (95 casts) — plus spot checks against **S4P 2018**. See "Validation Status" below for cross-cruise numbers.

## What Is Implemented

The full processing pipeline is implemented end-to-end:

**Ingestion:**
- `load_rdi()` — reads Teledyne RDI Workhorse PD0 binary files (`.000`)
- `load_ctd()` — reads SeaBird CNV files (ASCII and binary formats) and UH/CLIVAR-archive CTD time-series (A16N-style, GPS-dday time base), parsing pressure, temperature, salinity, and optionally GPS lat/lon
- `assign_bin_depths()` — converts ADCP bin ranges to absolute water depths using CTD pressure (Saunders p→z)
- `estimate_ctd_adcp_lag()` — measures the CTD-ADCP clock offset (besttlag equivalent)
- `compute_ship_velocity()` — derives ship velocity from GPS position time series
- SBE hex time-series decoder and SADCP NetCDF loader (shipboard ADCP profiles for the inverse constraint)

**Coordinate transforms:**
- `beam2earth()` — converts beam-coordinate ADCP velocities to Earth frame using heading, pitch, roll; supports the gimbaled heading correction (Janus geometry, RDI Workhorse convention) and 3-beam (single-missing-beam) reconstruction
- `uvrot()` — rotates East/North velocities for magnetic declination
- Sound-speed correction (`sounds.m`/`press.m`/`getdpthi.m` ports)

**Quality editing:**
- `edit_sidelobes()` — masks ADCP bins contaminated by acoustic sidelobes from the surface and bottom (LDEO `edit_data.m` convention)
- `edit_large_velocities()` — removes bins with horizontal speed > 2.5 m/s
- `edit_error_velocity()` — removes bins with anomalous error (4th-beam) velocity ("elim" editing)
- `edit_ppi()` — previous-ping-interference editing (`edit_data.m`)
- `edit_w_outliers()` — removes bins where vertical velocity deviates anomalously from near-instrument reference bins
- `build_ldeo_weights()` — full LDEO weight construction (correlation, echo, tilt, non-pinging removal)

**Inverse velocity solution:**
- `prepare_superensembles()` — depth-window averaging into super-ensembles with reference-bin subtraction; replicates `prepinv.m`
- `compute_inverse()` — constrained least-squares inversion producing `u(z)`, `v(z)` profiles; replicates `getinv.m`; supports GPS barotropic constraint, shipboard ADCP (SADCP) constraint, acoustic bottom-track constraint, and curvature smoothing (`smoofac`, off by default — matches LDEO's own archived setting on every validated cruise so far)
- `compute_shear()` — the "old-fashioned" shear-integration comparison method (`getshear2.m`)

**Output:**
- `write_ladcp_nc()` — NetCDF output writer (`ladcp2cdf.m` equivalent)

## MATLAB Equivalents

| MATLAB (LDEO_IX) | Python (this toolkit) | Status |
|---|---|---|
| `loadrdi.m` — read PD0 binary, UL/DL merge, weight construction | `load_rdi()`, `best_ul_shift()`, `build_ldeo_weights()` | Done |
| `janus5beam2earth()` (ADCPtools) | `beam2earth()` | Done (gimbaled, 3-beam) |
| `edit_data.m` — sidelobe / velocity / PPI editing | `edit_sidelobes()`, `edit_large_velocities()`, `edit_error_velocity()`, `edit_ppi()`, `edit_w_outliers()` | Done |
| `sounds.m` / `press.m` / `getdpthi.m` scaling | `soundspeed.py` | Done |
| `getshear2.m` — shear profiles | `compute_shear()` → `ShearProfile` | Done |
| `prepinv.m` — super-ensemble formation | `prepare_superensembles()` | Done |
| `getinv.m` — inverse solver | `compute_inverse()` | Done (deep-cast A16N divergence under investigation — see below) |
| `loadctdprof.m` / SBE hex decoder | `load_ctd()`, `assign_bin_depths()`, `sbe_hex.py` | Done |
| `loadsadcp.m` | SADCP NetCDF loader (`sadcp.py`) | Done |
| `getbtrack.m` — bottom-track processing | Integrated into `compute_inverse()` | Done |
| `fixcompass.m` / `checktilt.m` | `uvrot()` (declination only) | Partial |
| `ladcp2cdf.m` — NetCDF output | `write_ladcp_nc()` | Done |
| CLI (`ladcp process`, `ladcp check`) | `cli.py` | Stub — raises `NotImplementedError` |

The Perl-based LADCP_w software (vertical velocity, VKE) is a separate system not yet targeted.

## Validation Status

### P16N 2015 cast 003 — both targets MET (2026-07-11)

The 1000–2000 m anti-correlation described in earlier versions of this document is
resolved (root causes: a `medianan` vs `nanmedian` reference-selection bug, a
`beam2earth` up/down beam-matrix convention bug, a depth-registration offset, and
several ported edit/weight-construction steps — full trail in
`PROGRAMMERS_NOTES.md`). Current numbers, hard test assertions (not `xfail`):

**u RMSE = 0.0415 m/s, v RMSE = 0.0447 m/s** (target < 0.05 m/s). Updated
2026-07-17 after the rank-deficient-`lstsq` fix (`docs/HANDOVER.md`) shifted these
slightly from the original 0.0450/0.0333 — both still comfortably under target.

### Multi-cruise bulk validation (2026-07-11 to 2026-07-12)

To check the pipeline generalizes rather than being overfit to cast 003, it was run
against two entire cruises it never saw during development, using each cruise's own
archived SADCP/barotropic constraints (read back from the reference NetCDF attrs):

| Cruise | Casts | Pass both (u,v<0.05) | Pass u only | u RMSE median | u RMSE p90 |
|---|---|---|---|---|---|
| I7N 2018 | 124/124 | 53 (43%) | 20 | 0.043 | 0.117 |
| A16N 2013 | 95/95 | 15 (16%) | 4 | 0.34 | 1.13 |

Full per-cast tables: `docs/validation/BULK_VALIDATION_REPORT.md`.

**I7N**: strong evidence the core pipeline generalizes — u median RMSE close to the
tuning cast's own result on a 124-cast unseen cruise. Two open findings, one now
mostly resolved: (1) 10 casts "explode" numerically (u RMSE ~10⁶–10¹⁰, a bimodal
ill-conditioning failure mode, not gradual disagreement — depths 3440–4560 m,
mid-pack not deepest). **Root-caused and fixed 2026-07-17**: 9 of the 10 shared an
identical rank-deficient `lstsq` bug (an unconstrained depth bin's near-zero
singular value wasn't truncated by scipy's default cutoff); switching to
`numpy.linalg.lstsq(rcond=None)` fixes them (verified on 060, 042; the other 6
share the identical signature but weren't individually re-run). One cast (018) is
a different, still-open issue — well-conditioned but genuinely extreme output.
The pass-rate numbers in the table above predate this fix; see `docs/HANDOVER.md`.
(2) the 41 "marginal" casts are dominated by v-misses, consistent with the
un-ported `lanarrow` outlier-trim step (LDEO step 11).

**A16N**: shallow/mid casts (≤4000 m) perform well (19/36 pass u); all 59 casts deeper
than 4000 m fail (u RMSE ≥ 0.05), with error structure described in
`test_data/2013_A16N/DOWNLOAD_NOTES.md` — alternating-sign 0.3–0.8 m/s swings through
a weakly-constrained mid-column. As of 2026-07-17 the following have each been ruled
out with a direct check (details in `DOWNLOAD_NOTES.md`): the `ps.shear`/`smallfac`
regularization constraints (dead code paths in the exact LDEO software snapshot used
for this cruise — confirmed via the archived per-cast `.mat` parameter structs), the
minimum-norm solve method, mid-column observation-count starvation (Python actually
retains *more* data per bin than the reference, not less), and ship/winch heave
contamination (checked three independent ways: CTD descent-rate residual, ADCP tilt,
and W-anomaly editing rejection rate — none show an elevated signature on the deep
failing casts). The root mechanism is not yet found; the most concrete open lead is a
reproducible asymmetry where merging down-cast and up-cast observations helps shallow
casts a lot but appears to hurt deep casts relative to a (more crudely constrained)
down-only or up-only solve.

## Key Scientific Conventions

### DL and UL

The downlooker (DL) and uplooker (UL) are two Workhorse ADCPs co-mounted on the CTD rosette, looking down and up respectively. Each has its own PD0 file (e.g. `003DL000.000` and `003UL000.000`). Load each with `load_rdi()` and combine before passing to the inverse solver.

The UL is mounted face-up (inverted). Its pitch axis reads the **opposite sign** from the DL — pass `-rdi_ul.pitch` to `beam2earth()`.

### Timing

Clock drift between the DL and UL, and between either ADCP and the CTD, is a first-class processing concern. The current implementation aligns UL to DL by nearest-timestamp lookup. The MATLAB parameters `timoff` and `timoff_uplooker` for explicit clock-drift correction are not yet implemented.

### Coordinate frames

`RDIData.u/v` are in beam coordinates from the raw PD0 file. The `beam2earth()` call converts to Earth frame using heading, pitch, and roll. Magnetic declination correction requires a separate `uvrot()` call with the local East declination.

### Super-ensembles

The inverse solver works on "super-ensembles" — depth-window averages that group raw pings spanning one bin-length of CTD depth change (≈8 m on P16N cast 003; auto-computed per cast from the median bin spacing, matching MATLAB's `avdz` default). Within each window, the velocity is referenced to the two DL bins closest to the transducer face, removing the mean instrument velocity. What remains (the super-ensemble relative velocity `ru`) is approximately `u_ocean(z) − mean_u_instrument(window)`. The full-cast inverse then jointly solves for `u_ocean(z)` and `u_instrument(t)` across all windows.

### Boundary conditions

The inverse solver accepts three types of external velocity constraints:

- **Bottom-track**: ADCP acoustic return from the sea floor gives the instrument's absolute velocity near the bottom. The DL bottom-track is the strongest constraint and is enabled by default.
- **GPS barotropic**: GPS position fixes give the time-mean ship velocity over the cast. This constrains the depth-mean of `u_instrument` and is the primary absolute reference for the deep water column.
- **SADCP**: Shipboard ADCP near-surface measurements (0–300 m typically) constrain `u_ocean` in the surface layer and are used as an additional reference.

## Roadmap

1. ~~Ingestion~~ ✓
2. ~~Coordinate transforms~~ ✓
3. ~~QA editing~~ ✓
4. ~~Shear-based solution~~ ✓
5. ~~Inverse solution~~ ✓ — RMSE targets met on P16N cast 003 (u 0.0415, v 0.0447)
6. ~~NetCDF output (`ladcp2cdf.m` equivalent)~~ ✓
7. ~~Multi-cruise bulk validation (I7N, A16N)~~ ✓ — see "Multi-cruise bulk validation" above
8. A16N deep-cast (>4 km) divergence — **open investigation**, several leads ruled out (see above and `test_data/2013_A16N/DOWNLOAD_NOTES.md`)
9. ~~I7N "exploded" casts (10/124, numerical ill-conditioning)~~ ✓ 9/10 fixed 2026-07-17 (rank-deficient `lstsq`) — cast 018 remains, different cause, **open**
10. `lanarrow` outlier-trim port (LDEO step 11) — lead candidate for closing the I7N v-RMSE gap
11. QA diagnostics and plots — Planned
12. End-to-end CLI (`ladcp process <cast>`) — Planned
