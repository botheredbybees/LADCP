# Side-lobe Contamination Editing — Design Spec

**Date:** 2026-06-21
**Status:** Approved
**Reference:** `docs/legacy/edit_data.m` lines 142–186

---

## Context

The LADCP inverse solver integration test (`tests/integration/test_inverse_p16n_cast003.py`)
has an xfail RMSE of ~0.37 m/s vs a 0.05 m/s target. One documented pipeline gap is the
absence of any QC pass before super-ensemble formation. The LDEO MATLAB reference runs
`edit_data.m` (step 9 of `process_cast.m`) between data load and `prepinv.m`.

This spec covers the highest-impact, always-on operation from `edit_data.m`:
**side-lobe contamination editing**. The spike filter and PPI filter (both off by default)
and manual bad-ensemble blocking (cast-specific) are out of scope.

---

## What Side-lobe Contamination Is

An ADCP's transducer beams have side-lobes — spurious acoustic energy radiated at large
angles. These hit the sea surface (for an uplooker or shallow ADCP) and the sea floor (for
a downlooker near the bottom) and return echoes that contaminate velocity estimates in bins
near those boundaries.

The contaminated depth range is determined by beam geometry:

- **Surface sidelobe:** for a beam angle θ from vertical, the sidelobe reflection from a
  boundary at range `R` contaminates data within `(1 − cos θ) × R` of that boundary.
  With a margin of 1.5× cell size (matching Eric Firing's convention, reproduced in the
  LDEO code), the surface-contamination limit for an ensemble at depth `z` (negative) is:

  ```
  zlim_surface = (1 − cos θ) × z − 1.5 × cell_size_m
  ```

  Bins shallower than `zlim_surface` (i.e., `izm > zlim_surface`) are contaminated.

- **Bottom sidelobe:** with the ADCP at height-above-bottom `hab = z + zbottom` (positive),
  the bottom-contamination limit is:

  ```
  zlim_bottom = −zbottom + (1 − cos θ) × hab + 1.5 × cell_size_m
  ```

  Bins deeper than `zlim_bottom` (i.e., `izm < zlim_bottom`) are contaminated.

All depths follow the negative-down convention used throughout `EnsembleData`.

---

## Implementation

### New module: `src/ladcp/qa/editing.py`

**Public function:**

```python
def edit_sidelobes(
    ens: EnsembleData,
    *,
    zbottom: float | None = None,
    theta_deg: float = 20.0,
    cell_size_m: float,
) -> EnsembleData:
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `ens` | `EnsembleData` | Input data (negative-down depth convention) |
| `zbottom` | `float \| None` | Sea-floor depth in m (positive). If `None`, derived as `nanmedian(−ens.z + ens.hbot)`. |
| `theta_deg` | `float` | Beam angle from vertical in degrees. Default 20.0° for RDI Workhorse 300/600 kHz. |
| `cell_size_m` | `float` | Bin (cell) size in metres. Pass `rdi.blen_m`. |

**Returns:** A new `EnsembleData` (via `dataclasses.replace`) with contaminated weights
set to `NaN`. Depth and velocity arrays are unchanged.

**Algorithm:**

```
f = 1 − cos(radians(theta_deg))
margin = 1.5 × cell_size_m

# Auto-derive zbottom if not provided
if zbottom is None:
    zbottom = nanmedian(−ens.z + ens.hbot)   # positive floor depth
    if not isfinite(zbottom):
        zbottom = None                        # no BT data → skip bottom edit

# Surface sidelobe mask  (n_bins, n_ens) via broadcasting
zlim_surface = f × ens.z − margin            # (n_ens,)
bad_surface  = ens.izm > zlim_surface        # bins shallower than limit

# Bottom sidelobe mask (skipped if zbottom is None)
if zbottom is not None:
    hab        = ens.z + zbottom             # (n_ens,) height above floor
    zlim_bot   = −zbottom + f × hab + margin
    bad_bottom = ens.izm < zlim_bot
else:
    bad_bottom = zeros_like(ens.izm, dtype=bool)

bad        = bad_surface | bad_bottom
new_weight = ens.weight.copy()
new_weight[bad] = NaN

return dataclasses.replace(ens, weight=new_weight)
```

### Export

`edit_sidelobes` is added to `src/ladcp/qa/__init__.py`.

---

## Tests: `tests/test_editing.py`

| Test | What it checks |
|------|---------------|
| `test_surface_sidelobe_masks_shallow_bins` | Bins above `zlim_surface` become NaN; bins below are untouched |
| `test_bottom_sidelobe_masks_deep_bins` | Bins below `zlim_bottom` become NaN; bins above are untouched |
| `test_no_hbot_skips_bottom_edit` | All-NaN `hbot` with `zbottom=None` → only surface mask applied |
| `test_existing_nan_weights_preserved` | Pre-existing NaN weights remain NaN regardless of sidelobe status |
| `test_explicit_zbottom_overrides_auto` | Caller-supplied `zbottom` is used instead of auto-derive |
| `test_returns_new_ensemble_not_mutated` | Input `ens.weight` is not modified in-place |

---

## Integration wiring

In `tests/integration/test_inverse_p16n_cast003.py`, the `inverse_result` fixture is updated
to call `edit_sidelobes` between `EnsembleData` construction and `prepare_superensembles`:

```python
from ladcp.qa.editing import edit_sidelobes

ens = EnsembleData(...)
ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)
# zbottom auto-derived from ens.hbot
```

The xfail reason strings on `test_inverse_u_rmse` and `test_inverse_v_rmse` are updated
to remove the QC gap item once wired.

---

## Out of scope

- Spike filter (`edit_spike_filter`, off by default in MATLAB)
- Previous-ping interference filter (`edit_PPI`, off by default)
- Manual bad-ensemble blocking (`edit_dn_bad_ensembles`, cast-specific)
- Bottom-track HAB range filtering (from `prepinv.m`) — separate follow-up
