# User Notes

For shipboard and shore-based data officers working with LADCP data.

## Current State

**This software is not yet operational for end-to-end production processing.**

The core scientific pipeline — ingestion, coordinate transforms, QA editing, and the inverse velocity solver — is implemented and produces velocity profiles for full-water-column casts. However, the CLI (`ladcp process`, `ladcp check`) is not yet wired to the pipeline, and there is no NetCDF output stage. Profiles can be produced by calling the Python API directly (see example below).

Validation against LDEO MATLAB reference output gives u RMSE ≈ 0.07 m/s. The 0–1000 m range matches well (r ≈ +0.90). A known issue at 1000–2000 m is under active investigation. For production cruise processing, continue using LDEO_IX MATLAB software until the RMSE target (< 0.05 m/s) is met.

## Installation

Requires Python 3.11 or later and [uv](https://docs.astral.sh/uv/).

```bash
git clone <repo-url>
cd LADCP
uv sync
uv pip install -e .
```

Verify the installation:

```bash
uv run ladcp --help
```

## Input File Formats

| File type | Description | Example |
|---|---|---|
| `{stn}DL000.000` | Downlooker PD0 binary (beam mode, EX=0x04) | `003DL000.000` |
| `{stn}UL000.000` | Uplooker PD0 binary | `003UL000.000` |
| `{stn}_01.cnv` | CTD time-series (SeaBird CNV format, ASCII or binary) | `003_01.cnv` |
| `sadcp.npz` | Shipboard ADCP (optional, NumPy archive) | `sadcp_003.npz` |

The naming convention `{stn}/{stn}DL000.000` comes from the GO-SHIP acquisition convention documented in `test_data/ancillary/set_cast_params.m`.

## Running the Inverse Solver (Python API)

```python
from pathlib import Path
import numpy as np
from ladcp.ingestion.rdi import load_rdi
from ladcp.ingestion.ctd import load_ctd, assign_bin_depths, compute_ship_velocity
from ladcp.transforms.beam2earth import beam2earth, uvrot
from ladcp.qa.editing import edit_sidelobes, edit_large_velocities, edit_w_outliers
from ladcp.solution.inverse import EnsembleData, prepare_superensembles, compute_inverse

THETA = 20.0          # RDI Workhorse 300 kHz beam angle
DECLINATION = -12.3   # local East magnetic declination in degrees; negate for uvrot

base = Path("path/to/cast/files")
rdi    = load_rdi(base / "003DL000.000")
rdi_ul = load_rdi(base / "003UL000.000")
ctd    = load_ctd(base / "003_01.cnv")

# --- Downlooker ---
u_dl, v_dl, w_dl = beam2earth(rdi.u, rdi.v, rdi.w, rdi.e,
                               rdi.heading, rdi.pitch, rdi.roll,
                               THETA, gimbaled=True)
u_dl, v_dl = uvrot(u_dl, v_dl, DECLINATION)
z_m, izm_dl = assign_bin_depths(rdi, ctd, looker="down")
wt_dl = np.nanmean(rdi.corr.astype(float), axis=2) / 128.0

# --- Uplooker (negate pitch for inverted mounting) ---
u_ul, v_ul, w_ul = beam2earth(rdi_ul.u, rdi_ul.v, rdi_ul.w, rdi_ul.e,
                               rdi_ul.heading, -rdi_ul.pitch, rdi_ul.roll,
                               THETA, gimbaled=True)
u_ul, v_ul = uvrot(u_ul, v_ul, DECLINATION)
_, izm_ul = assign_bin_depths(rdi_ul, ctd, looker="up")
wt_ul = np.nanmean(rdi_ul.corr.astype(float), axis=2) / 128.0

# --- Time-align UL to DL ---
ul_idx = np.argmin(np.abs(rdi_ul.time_julian[:, None] - rdi.time_julian[None, :]), axis=0)
n_ul, n_dl = rdi_ul.nbin, rdi.nbin

# Combined array: UL reversed (shallow→deep) then DL
u_c  = np.vstack([u_ul[::-1, :][:, ul_idx],  u_dl])
v_c  = np.vstack([v_ul[::-1, :][:, ul_idx],  v_dl])
w_c  = np.vstack([w_ul[::-1, :][:, ul_idx],  w_dl])
wt_c = np.vstack([wt_ul[::-1, :][:, ul_idx], wt_dl])
izm_c = np.vstack([-izm_ul[::-1, :][:, ul_idx], -izm_dl])

izu = np.arange(n_ul - 1, -1, -1, dtype=int)
izd = np.arange(n_ul, n_ul + n_dl, dtype=int)

# --- Bottom track ---
bt_u, bt_v, bt_w = beam2earth(rdi.btrack_vel_ms[0], rdi.btrack_vel_ms[1],
                               rdi.btrack_vel_ms[2], rdi.btrack_vel_ms[3],
                               rdi.heading, rdi.pitch, rdi.roll, THETA, gimbaled=True)
bt_u, bt_v = uvrot(bt_u, bt_v, DECLINATION)
bvel  = np.stack([bt_u, bt_v, bt_w], axis=1)
bvels = np.full_like(bvel, 0.02)
hbot  = np.nanmean(rdi.btrack_range_m, axis=0)

# --- GPS ship velocity (if available from CTD ASCII file) ---
if ctd.lat is not None:
    slat = np.interp(rdi.time_julian, ctd.time_julian, ctd.lat)
    slon = np.interp(rdi.time_julian, ctd.time_julian, ctd.lon)
    u_ship, v_ship = compute_ship_velocity(ctd.lat, ctd.lon, ctd.time_julian)
else:
    slat = slon = np.full(rdi.nens, np.nan)
    u_ship = v_ship = None

ens = EnsembleData(u=u_c, v=v_c, w=w_c, weight=wt_c,
                   izm=izm_c, z=-z_m, time_jul=rdi.time_julian,
                   bvel=bvel, bvels=bvels, hbot=hbot,
                   izd=izd, izu=izu, slat=slat, slon=slon)

ens = edit_sidelobes(ens, theta_deg=THETA, cell_size_m=rdi.blen_m)
ens = edit_large_velocities(ens)
ens = edit_w_outliers(ens)

se     = prepare_superensembles(ens)
result = compute_inverse(se, u_ship=u_ship, v_ship=v_ship)

print(f"Profile: {len(result.z)} depth bins, {result.z.min():.0f}–{result.z.max():.0f} m")
print(f"u range: {result.u.min():.3f} to {result.u.max():.3f} m/s")
```

## Reading a PD0 File (Quick Inspection)

```python
from pathlib import Path
from ladcp.ingestion.rdi import load_rdi
import numpy as np

d = load_rdi(Path("003DL000.000"))
print(f"Ensembles: {d.nens}, Bins: {d.nbin}, bin length: {d.blen_m:.1f} m")
print(f"Cast duration: {(d.time_julian[-1] - d.time_julian[0]) * 24:.2f} hours")
print(f"Valid velocity cells: {np.isfinite(d.u).mean():.0%}")
print(f"Ensembles with bottom track: {np.isfinite(d.btrack_range_m).mean():.0%}")
```

## Running Tests

```bash
uv run pytest                              # unit tests only
TEST_DATA_DIR=test_data uv run pytest     # + integration tests (requires raw files)
```

## Docker

```bash
make build   # builds the Docker image
make test    # runs the test suite inside Docker
```

## Getting Help

For questions or bug reports, open an issue in the repository.

For LADCP processing questions or guidance on the LDEO_IX MATLAB workflow, refer to the GO-SHIP LADCP documentation and `docs/legacy/LADCP_processing.md`.
