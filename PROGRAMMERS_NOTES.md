# Programmer's Notes

Technical reference for developers working on this codebase. Read alongside the MATLAB source in `docs/legacy/` — the Python implementation is designed to be traceable to the MATLAB reference line-by-line.

## Architecture

Five layers in dependency order.

```
Ingestion  ──▶  Transforms  ──▶  Solution  ──▶  QA / Diagnostics  ──▶  CLI / API
(done)          (done)           (done)          (editing done)         (stubs)
```

### Module Map

```
src/ladcp/
├── __init__.py
├── cli.py                   Click app: `ladcp process` and `ladcp check` (stubs)
├── ingestion/
│   ├── __init__.py          exports load_rdi
│   ├── _pd0.py              parse_pd0(), low-level binary parser
│   ├── _types.py            RDIData dataclass
│   ├── rdi.py               load_rdi(path) → RDIData
│   └── ctd.py               load_ctd(), assign_bin_depths(), compute_ship_velocity()
├── transforms/
│   └── beam2earth.py        beam2earth() (gimbaled Janus), uvrot()
├── solution/
│   ├── shear.py             getshear2 equivalent
│   └── inverse.py           EnsembleData, SuperEnsemble, prepare_superensembles(),
│                            compute_inverse(), InverseParams, InverseResult
└── qa/
    └── editing.py           edit_sidelobes(), edit_large_velocities(), edit_w_outliers()
```

---

## Key Design Decisions

### Julian day convention: midnight-based (matching `julian.m`)

The MATLAB reference uses a non-standard Julian day convention where JD starts at **midnight**, not astronomical noon. This is implemented in `_pd0._to_julian()` using the Fliegel/Van Flandern algorithm:

```python
j = (146097 * c) // 4 + (1461 * yr) // 4 + (153 * mo + 2) // 5 + day + 1721119
return float(j) + hour_frac / 24.0
```

The difference from the Meeus astronomical noon formula is exactly 0.5 JD — using the wrong algorithm shifts every timestamp by 12 hours, corrupting DL/UL clock-drift corrections downstream.

### Byte offsets in the PD0 format

Field positions inside each block are *relative offsets* — every block starts at a base address found in the offset table at bytes 6+ of the ensemble header.

**Fixed leader (type 0x0000):**
- `nbin` at +7 (1 byte); `blen_cm` at +10; `blnk_cm` at +12
- `dist_cm` (distance to first bin centre) at +30 — **NOT +32**
- `serial` at +40 — **NOT +42**

**Variable leader (type 0x0080):**
- Timestamp (7 bytes) at +2
- After timestamp: 3-byte skip, then `sound_vel` at +12 — **NOT +14**
- `heading` at +16; pitch at +18; roll at +20

### Coordinate frame assumption

`RDIData.u/v/w/e` contain whatever the instrument recorded. For GO-SHIP casts the instrument is normally in beam mode (`EX=0x04`). `beam2earth()` in `transforms/beam2earth.py` converts beam coordinates to Earth frame; it must be called explicitly. The parser does not check the EX byte and does not auto-rotate.

### Velocity scaling and NaN

PD0 velocity fields are int16, unit 0.001 m/s, sentinel -32768 → `np.nan`. Bottom-track ranges are uint16 at 0.01 m/LSB with sentinel 0 → `np.nan`.

### UL pitch sign convention

The uplooker (UL) is mounted face-up (inverted). Its pitch sensor reads the **opposite sign** from the downlooker (DL) for the same physical tilt. Pass `-rdi_ul.pitch` to `beam2earth()` for the UL. Failure to negate UL pitch corrupts the gimbaled heading correction and rotates the Earth-frame velocity by the wrong heading.

### Combined DL+UL array layout

In `EnsembleData`, rows are ordered:
```
[UL bins reversed (shallowest to deepest water), DL bins (shallowest to deepest water)]
```
UL bins are reversed so that row index increases with water depth throughout the array — matching MATLAB's combined `d.ru` layout in `prepinv.m`. The index arrays `izu` and `izd` track which rows belong to each instrument.

### Super-ensemble reference bins

`prepare_superensembles()` uses DL bins 1 and 2 (0-indexed, skipping bin 0) as velocity reference, plus UL bins near the UL top. The reference subtraction removes the mean instrument velocity within each depth window, leaving `ru ≈ u_ocean[z_bin] - mean(u_ctd[window])`. Replicates `prepinv.m` lines ~500–600.

### Inverse solver sign convention

The solved system is `[A_ocean | A_ctd] * [u_ocean; u_ctd_neg] = d`. The solution `u_ctd_neg = -u_instrument`. After solving, `u_ctd = -m[n_zbins:]` restores the physical sign. This matches MATLAB's `dr.uctd = -real(uctd(:,1))'`.

### Magnetic declination rotation

`uvrot(u, v, angle_deg)` rotates velocities counter-clockwise by `angle_deg`. East magnetic declination is a clockwise heading shift, so the correction is `uvrot(u, v, -declination_deg)`. For P16N 2015: `uvrot(u, v, -12.318441)`.

---

## Validation: P16N Cast 003

The integration test in `tests/integration/test_inverse_p16n_cast003.py` runs the full pipeline and compares against the LDEO MATLAB reference output `test_data/2015_P16N/003.nc`.

Current status (2026-06-27):
- **u RMSE ≈ 0.0806 m/s** (target: < 0.05 m/s; tests are `xfail`)
- Correlation r ≈ +0.95 at 0–500 m, r ≈ +0.67 at 1000–2000 m (previously anti-correlated)
- Anti-correlation at 1000–2000 m is resolved; residual RMSE gap remains at 1000–4400 m

### Root Cause: Reference Subtraction Median vs MATLAB medianan (2026-06-27)

**Bug confirmed**: `prepare_superensembles()` used `np.nanmedian()` for the per-ensemble reference velocity; MATLAB uses `medianan(x, na=0)` which picks the `round(n/2)`-th sorted value.

For 4 reference bins (DL bin 1, DL bin 2, UL bin 1, UL bin 2):
- `np.nanmedian` returns the **average of the 2nd and 3rd sorted values**
- MATLAB `medianan(na=0)` returns the **2nd sorted value** (`round(4/2) = 2`)

When the UL instrument has a constant compass offset (~87° relative to DL — confirmed from heading data: downcast DL=115.9°/UL=29.7°, upcast DL=34.1°/UL=305.6°), the UL beam2earth produces Earth-frame velocities that are rotated ~87° from the DL's Earth-frame velocities. This causes DL and UL reference bins to give systematically different velocity values.

**Effect during upcast at 1800m**:
- DL reference bins: u ≈ −0.060 m/s
- UL reference bins: u ≈ +0.019 m/s (measuring a rotated component / different shear layer)
- Python `nanmedian`: (−0.058 + 0.018)/2 = **−0.020** → DL bin deref = −0.040 (wrong)
- MATLAB `medianan`: **−0.058** (2nd sorted = DL value) → DL bin deref ≈ 0 (correct)

The contaminated Python reference caused a ~0.10 m/s systematic bias and anti-correlation at 1000–2000 m depth.

**Fix applied**: Added `_ref_medianan()` in `inverse.py` (before `_window_boundaries()`), replacing the three `np.nanmedian` calls at lines 192–194 in `prepare_superensembles()`. The function picks `round(n_valid/2)`-th sorted value per column, matching MATLAB exactly.

**Post-fix status**:
- Anti-correlation at 1000–2000 m eliminated (r went from −0.40 to +0.67)
- Total u RMSE = 0.0806 (with bin masking), 0-1000m RMSE ~0.025 (excellent)
- Remaining gap at 1000–2000m (RMSE ≈ 0.10): likely from UL bins contributing observations contaminated by the DL-UL compass offset (UL measures a rotated Earth-frame velocity relative to DL); and 2000-4400m possible anti-correlation with UL-only depth coverage

### n_se Discrepancy Root Cause (2026-06-27)

**Bug**: test used hardcoded `dz=16.0`; MATLAB default `avdz = medianan(abs(diff(d.izm(:,1))))` = 8.0 m (bin spacing). This halved the super-ensemble count from ~947 to 524.

**Fix**: changed test to `dz=None` (auto-computed from `median(|diff(izm[:, 0])|)` = 8.0 m).

**MATLAB oversample not yet fully replicated**: MATLAB's `prepinv.m` expands each window symmetrically around its center with `i1l = length(i1)/2 * oversample` (default oversample=1). This creates overlapping windows and increases effective step to N+1 per window. `_window_boundaries` now implements `oversample=1.0` by default but produces n_se=947 vs MATLAB's 827. Remaining ~14% gap is from exact MATLAB rounding differences and the depth variable (MATLAB uses `d.izm(1,:)` = shallowest bin depth; Python uses `ens.z` = CTD depth; both change at ~1.09 m/ens so not a major factor).

### Remaining Gap: DL–UL Compass Offset and rotup2down

The remaining RMSE gap at 1000–4400 m is likely caused by the DL–UL compass misalignment affecting the inversion:

**Observation**: UL heading is consistently ~87° less than DL heading throughout the cast (downcast: 86.2°; upcast: 88.5°). After beam2earth, UL bins produce Earth-frame velocities rotated ~87° from DL's, meaning UL u-observations are approximately measuring v_ocean instead of u_ocean.

**MATLAB's `rotup2down=1`** (default from `default.m`) in `prepinv.m`: computes per-ensemble heading residual `hrotcomp = angle((UL_heading − hoff)/DL_heading)` where `hoff` is the cast-mean DL–UL offset. Rotates DL by −hrotcomp/2 and UL by +hrotcomp/2 to align them. Because `hoff ≈ 87°` removes the constant offset, `hrotcomp` is only the small per-ensemble variation (~±1°), so this correction has negligible effect on a constant-offset instrument pair.

**What Python is missing**: The `_ref_medianan` fix corrects the reference subtraction so DL bins have correct deref. But UL bins (with ~87° heading error) still provide wrong u_ocean observations in the inversion (they see v_ocean). With 24 UL bins vs 5 DL bins (after bin masking), the wrong UL observations dominate the solution at depths covered only by UL.

**Next steps to close the RMSE gap**:
1. Implement `rotup2down=1`: compute mean DL–UL heading offset, then rotate DL by −hrot/2 and UL by +hrot/2 per ensemble before `prepare_superensembles()`. Pass heading arrays (DL and UL) through `EnsembleData` or apply rotation in the test fixture.
2. Alternatively: apply a fixed 87° compass correction to UL headings in beam2earth (re-run UL beam2earth with `rdi_ul.heading + 87°`). This would fix the fundamental UL Earth-frame measurement error.
3. Investigate whether MATLAB's `offsetup2down` (step 12 re-form with first-guess) also contributes.

---

## Testing Approach

Tests live in `tests/`. Three levels of confidence:

**Unit tests** (`tests/test_pd0_parser.py`, `tests/test_inverse.py`, etc.) — no external files needed. Run in CI unconditionally.

**Integration tests** (`tests/integration/`) — real instrument files. Gated by `TEST_DATA_DIR` env var.
- `test_pd0_p16n_cast003.py` — P16N cast 003 file integrity and header checks
- `test_inverse_p16n_cast003.py` — full pipeline vs LDEO reference (145 tests, 2 `xfail`)

Run with:
```bash
TEST_DATA_DIR=test_data uv run pytest
```

---

## Extending the Ingestion Layer

To add support for a new PD0 block type:

1. Add a type-ID dispatch case in `_pd0.parse_pd0()`.
2. Write a `_read_<name>()` helper following the `data[offset:offset+N]` pattern.
3. Add the field to `RDIData` in `_types.py` with shape documentation.
4. Assemble the array in `rdi.load_rdi()`.
5. Add a unit test with a synthetic buffer and an integration assertion.

## Running Linter

```bash
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
```

Rules: E, F, I (imports), NPY (numpy), UP (pyupgrade). Line length 88.
