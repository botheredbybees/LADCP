# Handover: P16N Cast 003 RMSE Closure

**Date**: 2026-06-27  
**Status**: Root cause found and partially fixed. One remaining fix needed.

## What Was Done This Session

### Root cause confirmed and fixed: reference median mismatch

`prepare_superensembles()` used `np.nanmedian()` for the per-ensemble reference velocity `ur_t`. MATLAB uses `medianan(x, na=0)` which picks the **`round(n/2)`-th sorted value** — for 4 reference bins this is the 2nd sorted value, not the average of the 2nd and 3rd.

**Why this matters**: The DL compass reads ~115° and the UL compass reads ~30° at the same instant (87° offset throughout the cast). After beam2earth, the two instruments give systematically different Earth-frame velocities. During the upcast at 1800m, DL ref bins ≈ −0.060 m/s and UL ref bins ≈ +0.019 m/s. Python's nanmedian picks −0.020 (midpoint, wrong); MATLAB's medianan picks −0.058 (a DL value, correct). This 0.040 m/s reference error cascaded into ≈ −0.10 m/s systematic bias and anti-correlation at 1000–2000 m.

**Fix applied** (`src/ladcp/solution/inverse.py`):
- Added `_ref_medianan()` function before `_window_boundaries()` (lines ~62–83)
- Replaced the three `np.nanmedian(u_win[izr_valid], axis=0)` calls in `prepare_superensembles()` with `_ref_medianan(u_win[izr_valid])` (and same for v, w)

**Result**: Anti-correlation at 1000–2000 m eliminated (r: −0.40 → +0.67). Total u RMSE with bin masking = 0.0806. Without bin masking not yet tested this session.

## Remaining Work to Reach RMSE < 0.05

### Primary issue: UL compass offset corrupts inversion observations

The UL is consistently 87° offset from DL throughout the cast (confirmed from heading data). After beam2earth with the wrong UL heading, UL u-observations measure ≈ v_ocean instead of u_ocean. With 24 active UL bins vs 5 DL bins (after bin masking), the inversion is heavily influenced by the wrong UL data.

**Two equivalent fixes** (try option A first):

**Option A — Correct UL heading in beam2earth** (cleanest):
In the test fixture and integration test, replace:
```python
u_ul, v_ul, w_ul = beam2earth(..., rdi_ul.heading, ...)
```
with:
```python
# Apply mean DL-UL heading offset before beam2earth
ul_heading_corrected = rdi_ul.heading + 87.0  # adjust to match DL compass frame
u_ul, v_ul, w_ul = beam2earth(..., ul_heading_corrected, ...)
```
The exact offset (87°) can be computed from `np.nanmean(rdi.heading - rdi_ul.heading[ul_idx])` on the time-matched arrays. Try both +87 and using DL heading directly for UL.

**Option B — MATLAB rotup2down=1** (heavier, matches MATLAB exactly):
After beam2earth for both instruments, per ensemble compute `hrot = DL_heading - UL_heading - hoff` (residual after removing mean) and rotate DL by `−hrot/2`, UL by `+hrot/2`. See `docs/legacy/prepinv.m` lines 294–418 for full implementation.

### Bin masking status
The test fixture in `tests/integration/test_inverse_p16n_cast003.py` does NOT apply bin masking — that was tested in diagnostic scripts only. The integration test that drives the `xfail` RMSE check uses no bin masking. Verify which configuration gives better RMSE and apply consistently.

### Scripts to use for diagnosis
- `scripts/diag_downcast_only.py` — compares full vs downcast-only inversion
- `scripts/diag_ru_vs_ref.py` — compares SE ru values to reference profile
- `scripts/diag_colnorms.py` — column norm analysis of A_ocean matrix

## Files Changed This Session

| File | Change |
|------|--------|
| `src/ladcp/solution/inverse.py` | Added `_ref_medianan()`, replaced 3× `np.nanmedian` in reference subtraction |
| `PROGRAMMERS_NOTES.md` | Updated validation status, documented root cause and fix |
| `HANDOVER.md` | This file |
