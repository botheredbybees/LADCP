# Design: CTD Loading + Depth Assignment

**Date:** 2026-06-20
**Status:** Approved
**Phase:** Scoped slice — CTD time-series load → instrument depth → bin absolute depths.
**Prerequisite:** Ingestion layer (load_rdi, RDIData) and transforms layer (beam2earth) complete.

---

## Goal

Produce `z_m(nens)` (instrument depth at each ADCP ensemble) and `izm(nbin, nens)` (absolute depth of each ADCP bin at each ensemble) from a CTD time-series file and a loaded `RDIData` object. This is the missing prerequisite for the shear solution phase.

**Not in scope:** shear calculation, super-ensemble formation, sound-speed correction, CTD profile loading, GPS/nav data, SADCP.

---

## New file

All new code lives in one new module: `src/ladcp/ingestion/ctd.py`.

No changes to existing modules. `RDIData` is not mutated; callers receive numpy arrays.

---

## `CTDTimeSeries` dataclass

```python
@dataclass
class CTDTimeSeries:
    time_julian: NDArray[np.float64]    # (nctd,) Julian days
    pressure_dbar: NDArray[np.float64]  # (nctd,) positive down
    temp_c: NDArray[np.float64]         # (nctd,) NaN if absent
    salinity: NDArray[np.float64]       # (nctd,) NaN if absent
```

Salinity is included now because the shear phase will need it for sound speed correction.

---

## `load_ctd(path, **kwargs) → CTDTimeSeries`

### Format detection (in order)

1. Scan ASCII header lines (lines starting with `*` or `#`).
2. If any header line matches `# file_type = binary` → **binary SBE reader**.
3. Else if any header line matches `# name 0 =` → **ASCII SBE reader**.
4. Else → **generic ASCII reader** (uses kwargs; see below).

### Binary SBE reader

- Header ends at line `*END*`.
- `# nquan = N` gives number of float32 columns per scan.
- Column names from `# name N = varname: …` lines — parse `varname` before the colon.
- Data section: `N × float32` little-endian per scan. Use `np.frombuffer(data, dtype='<f4').reshape(-1, N)`.
- `# bad_flag = <value>` → mask those values to NaN before returning.

### ASCII SBE reader

- Same header parse as binary (header ends at `*END*`; `# nquan`, `# name N =` lines).
- Data section: whitespace-delimited float rows; use `np.loadtxt`.

### Generic ASCII reader (kwargs)

| kwarg | default | meaning |
|---|---|---|
| `skip_rows` | `0` | header lines to skip |
| `col_time` | `0` | 0-based column index for time |
| `col_pressure` | `1` | 0-based column for pressure (dbar) |
| `col_temp` | `2` | 0-based column for temperature |
| `col_salinity` | `None` | 0-based column for salinity; None → NaN array |
| `time_base` | `"elapsed_s"` | `"elapsed_s"` or `"julian"` |
| `time_start_julian` | `None` | required when `time_base="elapsed_s"` |

The I7N `.1Hz` files (30 header lines, 12 fields, elapsed-time) are covered by this path.

### Column name → semantic mapping (SBE readers)

| Variable prefix | Semantic role |
|---|---|
| `pr` (`prDM`, `prSM`, `prE`, …) | pressure |
| `t[0-9]` (`t090C`, `t190C`, …) | temperature |
| `sal` (`sal00`, `sal11`, …) | salinity |
| `timeJ` | time (Julian days) |
| `timeS` | time (elapsed seconds, requires `time_start_julian`) |

Unknown columns are ignored. If a required column (pressure, time) is absent, raise `ValueError` naming the missing variable.

### Time conversion

- `timeJ` → Julian days directly.
- `timeS` → `time_julian = timeS / 86400 + time_start_julian`. `time_start_julian` is extracted from the header `# start_time = <date string>` line if present, else the caller must supply it via kwargs.

---

## `assign_bin_depths(rdi, ctd, *, looker="down") → tuple[NDArray, NDArray]`

### Arguments

| Arg | Type | Notes |
|---|---|---|
| `rdi` | `RDIData` | from `load_rdi()` |
| `ctd` | `CTDTimeSeries` | from `load_ctd()` |
| `looker` | `str` | `"down"` (DL) or `"up"` (UL) |

### Steps

**Step 1 — Interpolate CTD pressure to ADCP time**

```python
p_interp = np.interp(rdi.time_julian, ctd.time_julian, ctd.pressure_dbar)
```

Clamps to boundary values outside the CTD time range (numpy default). Ensembles outside CTD range are expected to be small in number (deployment/recovery transients).

**Step 2 — Convert pressure to depth**

UNESCO 1983 approximation (no latitude dependency; error <0.5% at any depth):

```python
# depth_m positive down, pressure in dbar
z_m = (9.72659 * p - 2.2512e-5 * p**2 + 2.279e-10 * p**3 - 1.82e-15 * p**4) / (9.780318 * (1 + 5.2788e-3 * sin(lat)**2) + 1.092e-6 * p)
```

For the no-latitude fallback: `z_m ≈ p * 1.00445` (within 0.5% for 0–6000 m). Use the fallback unless the caller supplies `lat_deg` kwarg.

**Step 3 — Compute bin offset vector**

```python
bin_offsets = rdi.dist_m + np.arange(rdi.nbin) * rdi.blen_m  # (nbin,)
```

For DL (`looker="down"`): bins go deeper than instrument → offsets are positive.
For UL (`looker="up"`): bins go shallower than instrument → offsets are negative.

```python
sign = +1.0 if looker == "down" else -1.0
```

**Step 4 — Broadcast to (nbin, nens)**

```python
izm = z_m[np.newaxis, :] + sign * bin_offsets[:, np.newaxis]  # (nbin, nens)
```

No Python loops.

### Returns

`(z_m, izm)` where:
- `z_m`: `(nens,)` float64, instrument depth in metres positive-down
- `izm`: `(nbin, nens)` float64, absolute depth of each bin in metres positive-down

---

## Tests

### Unit tests — `tests/test_ctd_loader.py`

Run without any data files. Synthetic data only.

| Test | What it checks |
|---|---|
| `test_sbe_binary_header_parse` | Extract nquan, column names, bad_flag, time_start from a synthetic header string |
| `test_column_name_mapping_pressure` | `prDM`, `prSM`, `prE` all map to pressure |
| `test_column_name_mapping_time_julian` | `timeJ` → time Julian |
| `test_column_name_mapping_temperature` | `t090C`, `t190C` → temp |
| `test_bad_flag_masked_to_nan` | Values equal to bad_flag become NaN in output arrays |
| `test_generic_ascii_elapsed_time` | Round-trip through generic reader with `time_base="elapsed_s"` |
| `test_assign_bin_depths_shape` | `izm` shape == `(nbin, nens)` for synthetic inputs |
| `test_assign_bin_depths_down_deeper` | DL bins have increasing depth with bin index |
| `test_assign_bin_depths_up_shallower` | UL bins have decreasing depth with bin index |
| `test_assign_bin_depths_nan_propagation` | NaN in `z_m` propagates to `izm` for that ensemble |
| `test_format_dispatch_binary_flag` | Header with `# file_type = binary` routes to binary reader |
| `test_format_dispatch_sbe_ascii` | Header with `# name 0 =` but no binary flag routes to ASCII SBE reader |

### Integration tests — `tests/integration/test_ctd_p16n_cast003.py`

Gated by `TEST_DATA_DIR`. File: `test_data/2015_P16N/003_01.cnv`.

| Test | What it checks |
|---|---|
| `test_load_cnv_returns_ctd_time_series` | Returns `CTDTimeSeries`; no errors |
| `test_cnv_pressure_range` | `pressure_dbar` spans 0 to ~4367 dbar (header says max=4367.032) |
| `test_cnv_scan_count` | nctd == 307198 (from `# nvalues`) |
| `test_cnv_time_monotone` | `time_julian` is non-decreasing |
| `test_cnv_temperature_plausible` | All finite temps in [−2, 32] °C |
| `test_assign_depths_shape` | `izm.shape == (25, nens)` where nens from DL load |
| `test_assign_depths_instrument_range` | `z_m` spans 0 to ~4300 m |
| `test_assign_depths_bin0_deeper_than_instrument` | `izm[0, :] > z_m` for DL at all ensembles |

---

## Spec self-review

- No TBDs or placeholders.
- `time_start_julian` extraction from `# start_time` header line is specified clearly.
- Generic ASCII and SBE ASCII paths are distinct and non-overlapping.
- Pressure→depth formula specified with explicit fallback.
- Binary endianness explicit: little-endian float32 (`<f4`).
- `looker` parameter disambiguates DL/UL bin direction.
- Salinity included for future sound-speed use.
- All tests listed run without data files (unit) or are gated (integration).
