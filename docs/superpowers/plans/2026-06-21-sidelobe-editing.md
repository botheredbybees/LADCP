# Side-lobe Contamination Editing — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `edit_sidelobes()` in `src/ladcp/qa/editing.py`, add 6 unit tests, export it from `src/ladcp/qa/__init__.py`, and wire it into the integration test fixture between `EnsembleData` construction and `prepare_superensembles`.

**Architecture:** A pure-function `edit_sidelobes()` receives an `EnsembleData` and returns a new one (via `dataclasses.replace`) with contaminated weight cells set to `NaN`. Contamination zones are computed from beam geometry: `f = 1 − cos(θ)` times range from boundary, plus a 1.5×cell_size margin. Broadcasting over `(n_bins, n_ens)` vs `(n_ens,)` shapes handles per-ensemble thresholds without loops.

**Tech Stack:** Python 3.11, NumPy, `dataclasses.replace`, `math.cos`/`math.radians`. No new dependencies.

## Global Constraints

- Depth convention throughout: negative = below surface (both `z` and `izm` are ≤ 0)
- `zbottom` is **positive** (sea-floor depth in metres)
- `hab = ens.z + zbottom` — positive when ADCP is above the floor
- All QC functions return new `EnsembleData` via `dataclasses.replace()` — no mutation
- `rdi.blen_m` is the `cell_size_m` to pass from integration tests
- `THETA_DEG = 20.0` for RDI Workhorse 300 kHz (already in integration test)
- `EnsembleData` imported from `ladcp.solution.inverse`

---

### Task 1: Implement `edit_sidelobes()` with unit tests

**Files:**
- Create: `src/ladcp/qa/editing.py`
- Create: `tests/test_editing.py`

**Interfaces:**
- Produces: `edit_sidelobes(ens: EnsembleData, *, zbottom: float | None = None, theta_deg: float = 20.0, cell_size_m: float) -> EnsembleData`

- [ ] **Step 1: Write the 6 failing tests**

Create `tests/test_editing.py` with this complete content:

```python
"""Unit tests for ladcp.qa.editing.edit_sidelobes."""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from ladcp.qa.editing import edit_sidelobes
from ladcp.solution.inverse import EnsembleData


def _make_ens(
    izm: np.ndarray,
    z: np.ndarray,
    weight: np.ndarray | None = None,
    hbot: np.ndarray | None = None,
) -> EnsembleData:
    """Minimal EnsembleData for sidelobe editing tests."""
    n_bins, n_ens = izm.shape
    if weight is None:
        weight = np.ones((n_bins, n_ens))
    if hbot is None:
        hbot = np.full(n_ens, np.nan)
    return EnsembleData(
        u=np.zeros((n_bins, n_ens)),
        v=np.zeros((n_bins, n_ens)),
        w=np.zeros((n_bins, n_ens)),
        weight=weight,
        izm=izm,
        z=z,
        time_jul=np.zeros(n_ens),
        bvel=np.zeros((n_ens, 3)),
        bvels=np.full((n_ens, 3), 0.02),
        hbot=hbot,
        izd=np.arange(n_bins),
        izu=np.array([], dtype=int),
        slat=np.full(n_ens, np.nan),
        slon=np.full(n_ens, np.nan),
    )


def test_surface_sidelobe_masks_shallow_bins():
    # z=-100, theta=20°, cell=8m:
    #   f = 1 - cos(20°) ≈ 0.06031
    #   margin = 1.5 * 8 = 12.0
    #   zlim_surface = 0.06031 * (-100) - 12 ≈ -18.03
    # bin at izm=-5  → -5 > -18.03  → contaminated (shallower than limit)
    # bin at izm=-30 → -30 < -18.03 → clean
    izm = np.array([[-5.0], [-30.0]])   # (2, 1)
    z   = np.array([-100.0])            # (1,)
    ens = _make_ens(izm, z)
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "shallow bin should be masked"
    assert result.weight[1, 0] == 1.0, "deep bin should be untouched"


def test_bottom_sidelobe_masks_deep_bins():
    # z=-100, zbottom=150 (explicit), theta=20°, cell=8m:
    #   hab = -100 + 150 = 50
    #   zlim_bot = -150 + 0.06031*50 + 12 ≈ -134.98
    #   zlim_surface ≈ -18.03  (neither test bin is shallower)
    # bin at izm=-120 → -120 > -134.98 → clean
    # bin at izm=-140 → -140 < -134.98 → bottom-contaminated
    izm = np.array([[-120.0], [-140.0]])
    z   = np.array([-100.0])
    ens = _make_ens(izm, z)
    result = edit_sidelobes(ens, zbottom=150.0, theta_deg=20.0, cell_size_m=8.0)
    assert result.weight[0, 0] == 1.0, "mid-column bin should be untouched"
    assert np.isnan(result.weight[1, 0]), "near-bottom bin should be masked"


def test_no_hbot_skips_bottom_edit():
    # hbot all NaN → auto-derived zbottom = nanmedian(NaN) = NaN → skip bottom mask
    # bin at -5 is surface-contaminated (z=-100, zlim_surface≈-18.03)
    # bin at -140 would be bottom-contaminated with zbottom=150, but is NOT masked here
    izm = np.array([[-5.0], [-140.0]])
    z   = np.array([-100.0])
    ens = _make_ens(izm, z, hbot=np.array([np.nan]))
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "surface bin still masked"
    assert result.weight[1, 0] == 1.0, "bottom mask skipped when no BT data"


def test_existing_nan_weights_preserved():
    # A bin well inside the safe zone has a pre-existing NaN weight.
    # After editing it should still be NaN (copy preserves it).
    # z=-100, zlim_surface≈-18.03: bin at -50 is not surface-contaminated.
    izm    = np.array([[-50.0]])
    z      = np.array([-100.0])
    weight = np.array([[np.nan]])
    ens    = _make_ens(izm, z, weight=weight)
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "pre-existing NaN must not be cleared"


def test_explicit_zbottom_overrides_auto():
    # hbot=[200] → auto zbottom = 50+200 = 250 → zlim_bot ≈ -225.94
    # explicit zbottom=100  → hab=50,          zlim_bot ≈  -84.98
    # bin at izm=-90:
    #   surface: zlim_surface = 0.06031*(-50)-12 ≈ -15.02; -90 < -15.02 → NOT surface-masked
    #   with explicit zbottom=100: -90 < -84.98 → contaminated
    #   with auto  zbottom=250: -90 > -225.94  → clean
    # Passing explicit zbottom=100 should mask the bin.
    izm = np.array([[-90.0]])
    z   = np.array([-50.0])
    ens = _make_ens(izm, z, hbot=np.array([200.0]))
    result = edit_sidelobes(ens, zbottom=100.0, theta_deg=20.0, cell_size_m=8.0)
    assert np.isnan(result.weight[0, 0]), "explicit zbottom must be used, not auto-derived"


def test_returns_new_ensemble_not_mutated():
    # The function must return a new EnsembleData and not alter the input weight array.
    izm = np.array([[-5.0], [-50.0]])   # shallow bin will be surface-masked
    z   = np.array([-100.0])
    ens = _make_ens(izm, z)
    original_weight = ens.weight.copy()
    result = edit_sidelobes(ens, theta_deg=20.0, cell_size_m=8.0)
    np.testing.assert_array_equal(ens.weight, original_weight)  # not mutated
    assert result is not ens                                     # new object
    assert np.isnan(result.weight[0, 0])                        # sanity: mask applied
```

- [ ] **Step 2: Run tests to verify they all fail**

```
pytest tests/test_editing.py -v
```

Expected: 6 failures with `ModuleNotFoundError: No module named 'ladcp.qa.editing'`

- [ ] **Step 3: Implement `edit_sidelobes()`**

Create `src/ladcp/qa/editing.py` with this complete content:

```python
from __future__ import annotations

import dataclasses
import math

import numpy as np

from ladcp.solution.inverse import EnsembleData


def edit_sidelobes(
    ens: EnsembleData,
    *,
    zbottom: float | None = None,
    theta_deg: float = 20.0,
    cell_size_m: float,
) -> EnsembleData:
    """Zero-weight ADCP bins contaminated by surface and bottom acoustic side-lobes.

    Matches LDEO_IX edit_data.m lines 142–186 (Eric Firing convention).
    Returns a new EnsembleData; the input is not modified.
    """
    f = 1.0 - math.cos(math.radians(theta_deg))
    margin = 1.5 * cell_size_m

    if zbottom is None:
        derived = float(np.nanmedian(-ens.z + ens.hbot))
        zbottom = derived if math.isfinite(derived) else None

    # Surface mask: bins shallower than zlim_surface are contaminated.
    # zlim_surface shape: (n_ens,); ens.izm shape: (n_bins, n_ens) — broadcasts correctly.
    zlim_surface = f * ens.z - margin
    bad_surface = ens.izm > zlim_surface

    if zbottom is not None:
        hab = ens.z + zbottom                       # (n_ens,) height above floor
        zlim_bot = -zbottom + f * hab + margin      # (n_ens,)
        bad_bottom = ens.izm < zlim_bot
    else:
        bad_bottom = np.zeros_like(ens.izm, dtype=bool)

    new_weight = ens.weight.copy()
    new_weight[bad_surface | bad_bottom] = np.nan
    return dataclasses.replace(ens, weight=new_weight)
```

- [ ] **Step 4: Run tests to verify all 6 pass**

```
pytest tests/test_editing.py -v
```

Expected: 6 PASSED

- [ ] **Step 5: Commit**

```
git add src/ladcp/qa/editing.py tests/test_editing.py
git commit -m "feat: edit_sidelobes() surface and bottom sidelobe contamination masking"
```

---

### Task 2: Export, integration wiring, and xfail update

**Files:**
- Modify: `src/ladcp/qa/__init__.py` (add export)
- Modify: `tests/integration/test_inverse_p16n_cast003.py` (wire QC + update xfail strings)

**Interfaces:**
- Consumes: `edit_sidelobes` from `ladcp.qa.editing` (Task 1)
- Consumes: `EnsembleData`, `rdi.blen_m`, `THETA_DEG` (already present in integration test)

- [ ] **Step 1: Export `edit_sidelobes` from the qa package**

Replace the entire content of `src/ladcp/qa/__init__.py` with:

```python
from ladcp.qa.editing import edit_sidelobes

__all__ = ["edit_sidelobes"]
```

- [ ] **Step 2: Verify the export is importable**

```
python -c "from ladcp.qa import edit_sidelobes; print(edit_sidelobes)"
```

Expected: `<function edit_sidelobes at 0x...>`

- [ ] **Step 3: Wire `edit_sidelobes` into the integration test fixture**

In `tests/integration/test_inverse_p16n_cast003.py`, make these two changes:

**Change A — add import** (after the existing `from ladcp.transforms.beam2earth import beam2earth` line):

```python
from ladcp.qa.editing import edit_sidelobes
```

**Change B — call it** between `EnsembleData` construction and `prepare_superensembles` (i.e., between the closing `)` of `ens = EnsembleData(...)` and `se = prepare_superensembles(...)`):

```python
    ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)
```

The relevant section of the fixture should look like this after the change:

```python
    ens = EnsembleData(
        u=u_earth,
        v=v_earth,
        w=w_earth,
        weight=weight,
        izm=izm_neg,
        z=z_neg,
        time_jul=rdi.time_julian,
        bvel=bvel,
        bvels=bvels,
        hbot=hbot,
        izd=np.arange(rdi.nbin),
        izu=np.array([], dtype=int),
        slat=np.full(rdi.nens, np.nan),
        slon=np.full(rdi.nens, np.nan),
    )

    ens = edit_sidelobes(ens, theta_deg=THETA_DEG, cell_size_m=rdi.blen_m)

    se = prepare_superensembles(ens, dz=16.0)
    return compute_inverse(se)
```

- [ ] **Step 4: Update the xfail reason strings**

Both `test_inverse_u_rmse` and `test_inverse_v_rmse` have:

```python
@pytest.mark.xfail(strict=False, reason="Pipeline gaps: no edit_data.m QC pass, no GPS/SADCP constraint; remove once pipeline complete")
```

Replace with (both occurrences):

```python
@pytest.mark.xfail(strict=False, reason="Pipeline gap: no GPS/SADCP constraint; remove once pipeline complete")
```

- [ ] **Step 5: Run the full test suite (non-integration)**

```
pytest tests/ --ignore=tests/integration -v
```

Expected: all tests pass (including the 6 new editing tests)

- [ ] **Step 6: Commit**

```
git add src/ladcp/qa/__init__.py tests/integration/test_inverse_p16n_cast003.py
git commit -m "feat: wire edit_sidelobes into inverse pipeline; update xfail reason strings"
```
