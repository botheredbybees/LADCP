# User Notes

For shipboard and shore-based data officers working with LADCP data.

## Current State

**This software is not yet operational for end-to-end processing.**

The command-line tool (`ladcp process`, `ladcp check`) and the Python API for transforms, shear solutions, and QA plots do not work yet — they raise `NotImplementedError`. Only the file reader works.

What you can do today:
- Read a Teledyne RDI Workhorse PD0 binary file and inspect its contents in Python
- Verify that your `.000` files are readable and contain the expected number of ensembles
- Check heading, pitch, roll, and timing from the variable leader

For operational LADCP processing, continue using the LDEO_IX MATLAB software. This project will eventually replace it, but is not ready for cruise use.

## Installation

Requires Python 3.11 or later and [uv](https://docs.astral.sh/uv/).

```bash
# Clone the repository
git clone <repo-url>
cd LADCP

# Install dependencies and the ladcp package
uv sync
uv pip install -e .
```

Verify the installation:

```bash
uv run ladcp --help
```

You should see the help text. The subcommands (`process`, `check`) are listed but not yet functional.

## Input File Formats

The eventual processing pipeline will accept:

| File type | Description | Example |
|---|---|---|
| `{stn}DL000.000` | Downlooker PD0 binary | `003DL000.000` |
| `{stn}UL000.000` | Uplooker PD0 binary | `003UL000.000` |
| `{stn}.1Hz` | CTD time-series (30 header lines, 12 fields) | `003.1Hz` |
| `{stn}.cnv` | CTD profile (SeaBird `.cnv` format) | `003.cnv` |
| SADCP `.mat` | Shipboard ADCP (optional) | `sadcp.mat` |

The naming convention `{stn}/{stn}DL000.000` (station-named subdirectory) comes from the GO-SHIP acquisition convention documented in `test_data/ancillary/set_cast_params.m`.

## Reading a PD0 File (Python)

While the CLI is not operational, you can inspect files directly in Python:

```python
from pathlib import Path
from ladcp.ingestion import load_rdi
import numpy as np

# Load the downlooker file
d = load_rdi(Path("003DL000.000"))

# Basic cast summary
print(f"Ensembles: {d.nens}")
print(f"Bins: {d.nbin}, bin length: {d.blen_m:.1f} m")
print(f"First bin centre: {d.dist_m:.1f} m from transducer")
print(f"Cast duration: {(d.time_julian[-1] - d.time_julian[0]) * 24:.2f} hours")

# Heading statistics
print(f"Mean heading: {np.nanmean(d.heading):.1f}°")

# Velocity data coverage
finite_frac = np.isfinite(d.u).mean()
print(f"Valid velocity cells: {finite_frac:.0%}")

# Bottom-track availability
bt_frac = np.isfinite(d.btrack_range_m).mean()
print(f"Ensembles with bottom track: {bt_frac:.0%}")
```

## What the Output Will Look Like (Planned)

When processing is complete, `ladcp process` will produce a NetCDF file matching the LDEO_IX output format defined by `docs/legacy/ladcp2cdf.m`. The primary variables are:

- `u`, `v` — East and North velocity profiles (m/s)
- Depth bins (m)
- Quality flags and ensemble statistics

The reference output for cast 002 of the I7N 2018 cruise is at `test_data/data/002.nc` — this is the target the software must reproduce.

## Docker

A Dockerfile is provided for reproducible execution. Build and run tests with:

```bash
make build   # builds the Docker image
make test    # runs the test suite inside Docker
```

## Running Tests

To verify your installation works with the provided test data:

```bash
TEST_DATA_DIR=test_data uv run pytest tests/integration/
```

This requires raw PD0 files in `test_data/2015_P16N/`. If you have `.000` files from another cruise, you can point `TEST_DATA_DIR` at any directory containing Teledyne RDI Workhorse files and the basic file-integrity tests will run.

## Getting Help

This software is in active development. For questions or bug reports, contact the maintainer or open an issue in the repository.

For LADCP processing questions or guidance on the LDEO_IX MATLAB workflow, refer to the GO-SHIP LADCP documentation and the cruise report in `docs/legacy/LADCP_processing.md`.
