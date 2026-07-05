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

Current status (2026-07-05, `scripts/diag_rmse_strata.py`, convention fix applied,
rotup2down NOT applied — see below):

| stratum | n | u RMSE | v RMSE | r(u) |
|---|---|---|---|---|
| TOTAL | 520 | 0.0678 | 0.0573 | +0.483 |
| 0–1000 m | 120 | 0.0142 | 0.0195 | +0.970 |
| 1000–2000 m | 121 | 0.1034 | 0.0352 | +0.748 |
| 2000–3000 m | 122 | 0.0463 | 0.0735 | +0.576 |
| 3000–4500 m | 157 | 0.0720 | 0.0737 | +0.014 |

Target: TOTAL RMSE < 0.05 m/s (tests are `xfail`; not yet met). 0–1000 m is already
excellent; the gap is concentrated at depth, especially 1000–2000 m u RMSE.

### UL transform-convention fix (2026-07-05)

The `~87°`/`~90°` DL–UL heading disagreement previously blamed on a compass fault
(see the retired root-cause note below) was in fact a Python transform-convention
bug: `beam2earth` was not applying the correct up/down beam-matrix sign convention
for the uplooker (per `loadrdi.m::b2earth`, which uses a different beam→instrument
matrix depending on `beams_up`). Fixing the convention (commit `f3569c4`) makes each
instrument's Earth-frame velocity correct using its OWN heading — no compass-angle
hardcode needed. `scripts/diag_ul_dl_rotation.py` (E1) confirms the residual UL→DL
rotation after the fix is noise-dominated around 0° (circ mean −8° with prior 2 fold
variability, downcast +5°, upcast −35°, model rho ~0.1–0.2 — i.e. not a coherent
rotation), not a systematic ~87–90° bias. See `VALIDATION_PLAN.md` for the full
reasoning that led to retiring the "hardcode +87°" option.

### rotup2down: implemented, tested, tried, does not help this cast (2026-07-05)

`rotup2down()` in `src/ladcp/solution/inverse.py` is a faithful line-by-line port of
`prepinv.m` lines 294–418 (`rotup2down==1`, harmonize per-ensemble DL/UL heading
fluctuation after removing the cast-mean offset `hoff` via `compoff`). Verified
against the legacy source, including the MATLAB/Python `uvrot` sign-convention
inversion (MATLAB's `uvrot` negates its angle internally; Python's does not, so a
MATLAB call `uvrot(x,y,-hrot/2)` is replicated by Python `uvrot(x,y,+hrot/2)`) and
the bottom-track rotation using the raw (non-NaN-guarded) `hrot`. Three unit tests in
`tests/test_inverse.py` cover the no-op constant-offset case, the fluctuation-split
case, and non-mutation.

Wired into both the integration test fixture and `diag_rmse_strata.py` (`rot=True`
config) and measured:

| config | TOTAL u RMSE | TOTAL v RMSE | r(u) | 0–1000m u RMSE |
|---|---|---|---|---|
| convention fix only (baseline) | 0.0678 | 0.0573 | +0.483 | 0.0142 |
| + rotup2down | 0.0755 | 0.0552 | +0.381 | 0.0327 |

rotup2down **worsens** u RMSE almost everywhere, most severely in the 0–1000 m
stratum (already the best-performing region), and only marginally improves v RMSE
in the two deepest strata. A sign-bug was ruled out empirically: running the same
comparison with the DL/UL rotation directions swapped gives 0.0754 TOTAL u RMSE
(no material difference from the current sign, and still worse than baseline) — if
the correction were undoing a real physical misalignment, flipping the sign would
swing the result the other way instead of landing on the same degradation. The
per-ensemble `hrot` residual itself is small (mean 0.7°, std 7.3°, cast 003) and is
most plausibly per-instrument compass jitter rather than a genuine mechanical
flex between the two frames — "correcting" for it injects that jitter into u/v
instead of removing a real signal.

**Decision: not committed.** The implementation is stashed (`git stash list`, not in
the working tree) rather than discarded, in case a future cast or instrument pair
(e.g. Nuyina's own rosette, with different mounting rigidity) shows a real
per-ensemble heading disagreement worth correcting. Do not re-apply it here without
new evidence the residual is signal, not noise.

### Root Cause: Reference Subtraction Median vs MATLAB medianan (2026-06-27)

**Bug confirmed**: `prepare_superensembles()` used `np.nanmedian()` for the per-ensemble reference velocity; MATLAB uses `medianan(x, na=0)` which picks the `round(n/2)`-th sorted value.

For 4 reference bins (DL bin 1, DL bin 2, UL bin 1, UL bin 2):
- `np.nanmedian` returns the **average of the 2nd and 3rd sorted values**
- MATLAB `medianan(na=0)` returns the **2nd sorted value** (`round(4/2) = 2`)

At the time, the ~87° DL/UL heading disagreement (downcast DL=115.9°/UL=29.7°,
upcast DL=34.1°/UL=305.6°) was attributed to a compass fault; it was later found to
be a `beam2earth` transform-convention bug, fixed in `f3569c4` (see "UL
transform-convention fix" above) — both instruments' raw headings were correct. At
the time of this fix the wrong convention was still in place, so UL beam2earth
still produced Earth-frame velocities rotated ~87° from DL's, and this is what
caused DL and UL reference bins to give systematically different velocity values.

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
- Remaining gap at 1000–2000m (RMSE ≈ 0.10): at the time, attributed to UL bins
  contaminated by the (mis-diagnosed) compass offset. With the transform-convention
  fix now applied, current TOTAL u RMSE is 0.0678 (see "Validation: P16N Cast 003"
  status table above) — improved but still above the 0.05 target; the residual gap
  is not yet root-caused (rotup2down does not close it — see below).

### n_se Discrepancy Root Cause (2026-06-27)

**Bug**: test used hardcoded `dz=16.0`; MATLAB default `avdz = medianan(abs(diff(d.izm(:,1))))` = 8.0 m (bin spacing). This halved the super-ensemble count from ~947 to 524.

**Fix**: changed test to `dz=None` (auto-computed from `median(|diff(izm[:, 0])|)` = 8.0 m).

**MATLAB oversample not yet fully replicated**: MATLAB's `prepinv.m` expands each window symmetrically around its center with `i1l = length(i1)/2 * oversample` (default oversample=1). This creates overlapping windows and increases effective step to N+1 per window. `_window_boundaries` now implements `oversample=1.0` by default but produces n_se=947 vs MATLAB's 827. Remaining ~14% gap is from exact MATLAB rounding differences and the depth variable (MATLAB uses `d.izm(1,:)` = shallowest bin depth; Python uses `ens.z` = CTD depth; both change at ~1.09 m/ens so not a major factor).

### Remaining Gap: DL–UL Compass Offset and rotup2down — superseded (2026-07-05)

**This section is retired.** The "~87° DL–UL compass offset" was not a compass fault;
it was a Python `beam2earth` transform-convention bug (wrong up/down beam-matrix sign
for the uplooker), fixed in commit `f3569c4`. Hardcoding a compass-angle correction
(previously proposed as "Option A") would have been wrong — see `VALIDATION_PLAN.md`
for the evidence. `rotup2down` ("Option B") was subsequently implemented and measured;
see "rotup2down: implemented, tested, tried, does not help this cast" above for the
current status (not committed — it worsens RMSE on this cast). The remaining RMSE gap
at depth is not yet root-caused; see `VALIDATION_PLAN.md` Phase 2 for the solver-only /
transform-only harness planned to isolate it.

---

## Testing Approach

Tests live in `tests/`. Three levels of confidence:

**Unit tests** (`tests/test_pd0_parser.py`, `tests/test_inverse.py`, etc.) — no external files needed. Run in CI unconditionally.

**Integration tests** (`tests/integration/`) — real instrument files. Gated by `TEST_DATA_DIR` env var.
- `test_pd0_p16n_cast003.py` — P16N cast 003 file integrity and header checks
- `test_inverse_p16n_cast003.py` — full pipeline vs LDEO reference (190 tests total repo-wide with `TEST_DATA_DIR` set, 8 skipped, 2 `xfail` on the RMSE checks)

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
