# Oceanographer's Notes

Context for physical oceanographers evaluating or contributing to this toolkit.

## What This Is

A Python re-implementation of the LDEO LADCP processing workflow, targeting the same outputs as the MATLAB software used for GO-SHIP and similar hydrographic cruises. The goal is to reproduce known-good processed results numerically before adding any new capabilities.

The validation dataset is **2015 P16N, cast 003**, processed by A.M. Thurnherr (LDEO). The reference output `test_data/2015_P16N/003.nc` is the benchmark this software currently works against.

## What Is Implemented

The full processing pipeline is implemented end-to-end:

**Ingestion:**
- `load_rdi()` — reads Teledyne RDI Workhorse PD0 binary files (`.000`)
- `load_ctd()` — reads SeaBird CNV files (ASCII and binary formats), parsing pressure, temperature, salinity, and optionally GPS lat/lon
- `assign_bin_depths()` — converts ADCP bin ranges to absolute water depths using CTD pressure
- `compute_ship_velocity()` — derives ship velocity from GPS position time series

**Coordinate transforms:**
- `beam2earth()` — converts beam-coordinate ADCP velocities to Earth frame using heading, pitch, roll; supports the gimbaled heading correction (Janus geometry, RDI Workhorse convention)
- `uvrot()` — rotates East/North velocities for magnetic declination

**Quality editing:**
- `edit_sidelobes()` — masks ADCP bins contaminated by acoustic sidelobes from the surface and bottom (LDEO `edit_data.m` convention)
- `edit_large_velocities()` — removes bins with horizontal speed > 2.5 m/s
- `edit_w_outliers()` — removes bins where vertical velocity deviates anomalously from near-instrument reference bins

**Inverse velocity solution:**
- `prepare_superensembles()` — depth-window averaging into super-ensembles with reference-bin subtraction; replicates `prepinv.m`
- `compute_inverse()` — constrained least-squares inversion producing `u(z)`, `v(z)` profiles; replicates `getinv.m`; supports GPS barotropic constraint, shipboard ADCP (SADCP) constraint, acoustic bottom-track constraint, and smoothness regularisation

## MATLAB Equivalents

| MATLAB (LDEO_IX) | Python (this toolkit) | Status |
|---|---|---|
| `loadrdi.m` — read PD0 binary | `load_rdi()` → `RDIData` | Done |
| `janus5beam2earth()` (ADCPtools) | `beam2earth()` | Done (gimbaled) |
| `edit_data.m` — sidelobe / velocity editing | `edit_sidelobes()`, `edit_large_velocities()`, `edit_w_outliers()` | Done |
| `getshear2.m` — shear profiles | `compute_shear()` → `ShearProfile` | Done |
| `prepinv.m` — super-ensemble formation | `prepare_superensembles()` | Done |
| `getinv.m` — inverse solver | `compute_inverse()` | Done (RMSE work in progress) |
| `loadctdprof.m` | `load_ctd()`, `assign_bin_depths()` | Done |
| `loadsadcp.m` | SADCP fixture loader (tests only) | Done |
| `getbtrack.m` — bottom-track processing | Integrated into `compute_inverse()` | Done |
| `fixcompass.m` / `checktilt.m` | `uvrot()` (declination only) | Partial |
| `ladcp2cdf.m` — NetCDF output | Not yet implemented | Planned |

The Perl-based LADCP_w software (vertical velocity, VKE) is a separate system not yet targeted.

## Validation Status

### Current numbers (P16N 2015, cast 003)

The inverse solver is running and producing profiles for the full 4500 m cast.

**Overall:** u RMSE = 0.072 m/s (target: < 0.05 m/s for GO-SHIP quality)

| Depth range | u correlation | Notes |
|---|---|---|
| 0–500 m | +0.91 | Good |
| 500–1000 m | +0.89 | Good |
| 1000–1500 m | **−0.39** | Anti-correlated — under investigation |
| 1500–2000 m | **−0.41** | Anti-correlated — under investigation |
| 2000–3000 m | +0.20 | Moderate |
| 3000–4500 m | −0.28 | Noisy |

The anti-correlation at 1000–2000 m is the dominant driver of excess RMSE. This depth range corresponds to a portion of the cast where the CTD instrument appears to be drifting rapidly eastward (u_instrument ≈ +0.10–0.15 m/s) relative to the ocean. The inverse solver is not correctly separating instrument motion from ocean velocity at those depths.

**What has been ruled out:** GPS constraint weighting (removing GPS does not change the anti-correlation), bottom-track constraint, adaptive velocity error weighting.

**What is suspected:** A difference in super-ensemble formation between MATLAB (827 super-ensembles, inferred dz ≈ 8 m) and Python (524 super-ensembles, dz = 16 m hardcoded). This requires MATLAB's intermediate arrays (`di.ru`, `di.izm`) for direct comparison.

See `docs/superpowers/plans/2026-06-22-rmse-closure.md` for the full investigation writeup and what data is needed to proceed.

## Key Scientific Conventions

### DL and UL

The downlooker (DL) and uplooker (UL) are two Workhorse ADCPs co-mounted on the CTD rosette, looking down and up respectively. Each has its own PD0 file (e.g. `003DL000.000` and `003UL000.000`). Load each with `load_rdi()` and combine before passing to the inverse solver.

The UL is mounted face-up (inverted). Its pitch axis reads the **opposite sign** from the DL — pass `-rdi_ul.pitch` to `beam2earth()`.

### Timing

Clock drift between the DL and UL, and between either ADCP and the CTD, is a first-class processing concern. The current implementation aligns UL to DL by nearest-timestamp lookup. The MATLAB parameters `timoff` and `timoff_uplooker` for explicit clock-drift correction are not yet implemented.

### Coordinate frames

`RDIData.u/v` are in beam coordinates from the raw PD0 file. The `beam2earth()` call converts to Earth frame using heading, pitch, and roll. Magnetic declination correction requires a separate `uvrot()` call with the local East declination.

### Super-ensembles

The inverse solver works on "super-ensembles" — depth-window averages that group raw pings spanning ≈16 m of CTD depth change. Within each window, the velocity is referenced to the two DL bins closest to the transducer face, removing the mean instrument velocity. What remains (the super-ensemble relative velocity `ru`) is approximately `u_ocean(z) − mean_u_instrument(window)`. The full-cast inverse then jointly solves for `u_ocean(z)` and `u_instrument(t)` across all windows.

### Boundary conditions

The inverse solver accepts three types of external velocity constraints:

- **Bottom-track**: ADCP acoustic return from the sea floor gives the instrument's absolute velocity near the bottom. The DL bottom-track is the strongest constraint and is enabled by default.
- **GPS barotropic**: GPS position fixes give the time-mean ship velocity over the cast. This constrains the depth-mean of `u_instrument` and is the primary absolute reference for the deep water column.
- **SADCP**: Shipboard ADCP near-surface measurements (0–300 m typically) constrain `u_ocean` in the surface layer and are used as an additional reference.

## Data Requirements for Continued Validation

To resolve the 1000–2000 m anti-correlation, the following MATLAB intermediate arrays are needed for P16N cast 003:

```matlab
% Run in MATLAB after processing cast 003 with LDEO_IX
save('di_cast003.mat', 'di');   % prepinv.m output
save('dr_cast003.mat', 'dr');   % getinv.m output
```

Key variables: `di.ru`, `di.rv` (super-ensemble relative velocities, all bins × all SEs), `di.izm` (bin depths), `dr.uctd`, `dr.zctd` (instrument velocity time series).

## Roadmap

1. ~~Ingestion~~ ✓
2. ~~Coordinate transforms~~ ✓
3. ~~QA editing~~ ✓
4. ~~Shear-based solution~~ ✓
5. ~~Inverse solution (implemented; RMSE closure in progress)~~ ⚠
6. NetCDF output (`ladcp2cdf.m` equivalent) — Planned
7. QA diagnostics and plots — Planned
8. End-to-end CLI (`ladcp process <cast>`) — Planned
9. RMSE < 0.05 m/s on P16N cast 003 — **Blocked on MATLAB intermediate arrays**
