# Shear Solution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `compute_shear()` — a function that takes Earth-frame ADCP velocities and bin depths, computes central-difference shear, averages into depth bins with 2σ outlier editing, and integrates bottom-up to produce a zero-mean relative velocity profile matching the LDEO `getshear2.m` reference.

**Architecture:** Single module `src/ladcp/solution/shear.py` replaces the existing `shear_solution` stub with a `ShearProfile` dataclass and a `compute_shear()` top-level function backed by three private helpers: `_central_diff_shear`, `_bin_average_shear`, and `_integrate_shear`. Unit tests use synthetic arrays with known analytical answers; one integration test compares depth range and velocity magnitude against the P16N cast 003 NetCDF reference.

**Tech Stack:** Python 3.11, numpy (array ops), pytest (tests), netCDF4 (integration test only — already a project dependency)

## Global Constraints

- Python ≥ 3.11; numpy is the only new dependency (already present)
- Array axis convention: `(nbin, nens)` for 2-D fields, matching `RDIData` and `assign_bin_depths` output
- Depths are **positive downward** in meters (Python convention differs from MATLAB's negative depths)
- Functions return plain numpy arrays / dataclasses — no mutation of input objects
- Integration tests must be gated with `@pytest.mark.integration` and skip when `TEST_DATA_DIR` env var is absent
- All commits must pass `uv run pytest -x -q` before being recorded
- Run tests with: `uv run pytest tests/ -x -q`

---

### Task 1: `ShearProfile` dataclass and `compute_shear()` API skeleton

Replace the narrow `shear_solution` stub with the full API surface. Write failing tests first.

**Files:**
- Modify: `src/ladcp/solution/shear.py` (replace entire stub)
- Create: `tests/test_shear.py`

**Interfaces:**
- Produces: `ShearProfile` dataclass (fields below), `compute_shear(u, v, w, izm, weight, dz, *, stdf, weight_min) -> ShearProfile`

- [ ] **Step 1: Write the failing test**

Create `tests/test_shear.py`:

```python
"""Unit tests for shear solution — src/ladcp/solution/shear.py."""
import numpy as np
import pytest
from ladcp.solution.shear import ShearProfile, compute_shear


def _constant_shear_inputs(
    nbin: int = 10,
    nens: int = 20,
    du_dz: float = 1e-3,
    dv_dz: float = -5e-4,
    dw_dz: float = 0.0,
    bin_spacing: float = 10.0,
    first_bin_depth: float = 10.0,
):
    """Synthetic field with exact linear shear du/dz everywhere."""
    izm = np.outer(
        np.arange(nbin) * bin_spacing + first_bin_depth,
        np.ones(nens),
    )  # (nbin, nens) positive depths
    u = izm * du_dz   # u = du_dz * z → shear = du_dz everywhere
    v = izm * dv_dz
    w = izm * dw_dz
    weight = np.ones((nbin, nens), dtype=np.float64)
    return u, v, w, izm, weight


def test_compute_shear_returns_shear_profile():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    assert isinstance(result, ShearProfile)


def test_shear_profile_fields_exist():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    for field in ("z", "u_shear", "v_shear", "w_shear",
                  "u_shear_err", "v_shear_err", "w_shear_err",
                  "n", "u_rel", "v_rel", "w_rel"):
        assert hasattr(result, field), f"ShearProfile missing field: {field}"


def test_shear_profile_arrays_same_length():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    nz = len(result.z)
    for arr in (result.u_shear, result.v_shear, result.w_shear,
                result.u_shear_err, result.v_shear_err, result.w_shear_err,
                result.n, result.u_rel, result.v_rel, result.w_rel):
        assert arr.shape == (nz,), f"Expected ({nz},), got {arr.shape}"


def test_z_axis_starts_at_half_dz():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    assert result.z[0] == pytest.approx(5.0)   # dz/2 = 5 m


def test_integrated_profile_is_zero_mean():
    u, v, w, izm, weight = _constant_shear_inputs()
    result = compute_shear(u, v, w, izm, weight, dz=10.0)
    assert np.mean(result.u_rel) == pytest.approx(0.0, abs=1e-12)
    assert np.mean(result.v_rel) == pytest.approx(0.0, abs=1e-12)
```

- [ ] **Step 2: Run test to verify it fails**

```
uv run pytest tests/test_shear.py -x -q
```

Expected: `ImportError` or `cannot import name 'ShearProfile'` — the stub only defines `shear_solution`.

- [ ] **Step 3: Replace the stub with the ShearProfile dataclass and a compute_shear skeleton**

Replace the entire contents of `src/ladcp/solution/shear.py`:

```python
"""Shear-based horizontal velocity solution.

Reference: docs/legacy/getshear2.m (Visbeck, LDEO 1997).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ladcp._typing import NDArray


@dataclass
class ShearProfile:
    """Depth-binned shear profile and integrated relative velocity.

    All arrays have shape (nz,) where nz = ceil(z_max / dz).
    z[k] is the centre of the k-th depth bin in metres (positive down).
    """

    z: NDArray           # (nz,) bin centre depths, m
    u_shear: NDArray     # (nz,) mean du/dz per bin, s⁻¹
    v_shear: NDArray     # (nz,) mean dv/dz per bin, s⁻¹
    w_shear: NDArray     # (nz,) mean dw/dz per bin, s⁻¹
    u_shear_err: NDArray # (nz,) 1-σ shear uncertainty, s⁻¹
    v_shear_err: NDArray # (nz,) 1-σ shear uncertainty, s⁻¹
    w_shear_err: NDArray # (nz,) 1-σ shear uncertainty, s⁻¹
    n: NDArray           # (nz,) int — estimates used per bin after editing
    u_rel: NDArray       # (nz,) relative eastward velocity, m/s (zero-mean)
    v_rel: NDArray       # (nz,) relative northward velocity, m/s (zero-mean)
    w_rel: NDArray       # (nz,) relative vertical velocity, m/s (zero-mean)


def compute_shear(
    u: NDArray,
    v: NDArray,
    w: NDArray,
    izm: NDArray,
    weight: NDArray,
    dz: float = 10.0,
    *,
    stdf: float = 2.0,
    weight_min: float = 0.1,
) -> ShearProfile:
    """Compute depth-binned shear profile and integrate to relative velocities.

    Parameters
    ----------
    u, v, w : ndarray, shape (nbin, nens)
        Earth-frame velocity components (East, North, Up) in m/s.
    izm : ndarray, shape (nbin, nens)
        Depth of each ADCP bin in metres (positive down).
    weight : ndarray, shape (nbin, nens)
        Quality weight in [0, 1]. Bins with weight ≤ weight_min are excluded.
    dz : float
        Output depth-bin size in metres. Default 10 m.
    stdf : float
        Outlier rejection threshold in units of std. Default 2 (2σ editing).
    weight_min : float
        Minimum weight for a bin to be included. Default 0.1.

    Returns
    -------
    ShearProfile
        Depth-binned shear and integrated relative velocity profile.
    """
    raise NotImplementedError
```

- [ ] **Step 4: Run tests to confirm the right failures**

```
uv run pytest tests/test_shear.py -x -q
```

Expected: first test `test_compute_shear_returns_shear_profile` fails with `NotImplementedError`. The dataclass and import succeed.

- [ ] **Step 5: Commit the skeleton**

```bash
git add src/ladcp/solution/shear.py tests/test_shear.py
git commit -m "feat: add ShearProfile dataclass and compute_shear() skeleton"
```

---

### Task 2: Implement `_central_diff_shear()` — Stage 1

Compute stride-2 central differences and apply the weight mask. The MATLAB `diff2(x)` computes `x[2:] - x[:-2]` (stride 2, not stride 1), and the result is padded with `NaN` at first and last row to maintain the original shape.

**Files:**
- Modify: `src/ladcp/solution/shear.py` (add private helper)
- Modify: `tests/test_shear.py` (add helper tests)

**Interfaces:**
- Produces: `_central_diff_shear(u, v, w, izm, weight_mask) -> tuple[NDArray, NDArray, NDArray]`
  - Returns `(shear_u, shear_v, shear_w)`, each shape `(nbin, nens)` with NaN at rows 0 and -1 and wherever weight_mask is NaN

- [ ] **Step 1: Write the failing tests for `_central_diff_shear`**

Add to `tests/test_shear.py`:

```python
from ladcp.solution.shear import _central_diff_shear


def test_central_diff_shear_exact_gradient():
    """Linear u = du_dz * z → shear = du_dz everywhere in interior."""
    du_dz = 1e-3
    nbin, nens = 6, 4
    bin_spacing = 8.0
    izm = np.outer(np.arange(nbin) * bin_spacing + 8.0, np.ones(nens))
    u = izm * du_dz
    v = np.zeros_like(u)
    w = np.zeros_like(u)
    weight_mask = np.ones((nbin, nens))
    su, sv, sw = _central_diff_shear(u, v, w, izm, weight_mask)
    # Interior bins (rows 1 to nbin-2) should equal du_dz exactly
    interior = su[1:-1, :]
    assert np.allclose(interior, du_dz, rtol=1e-10)


def test_central_diff_shear_boundary_nan():
    """First and last row must be NaN (no neighbour on one side)."""
    u, v, w, izm, weight = _constant_shear_inputs()
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    assert np.all(np.isnan(su[0, :]))
    assert np.all(np.isnan(su[-1, :]))


def test_central_diff_shear_weight_mask_applies():
    """Bins with weight ≤ weight_min become NaN shear (caller pre-masks)."""
    u, v, w, izm, _ = _constant_shear_inputs(nbin=6, nens=4)
    # weight_mask is 1 everywhere except column 2 → NaN
    mask = np.ones_like(u)
    mask[:, 2] = np.nan
    su, _, _ = _central_diff_shear(u, v, w, izm, mask)
    assert np.all(np.isnan(su[:, 2]))


def test_central_diff_shear_output_shape():
    u, v, w, izm, weight = _constant_shear_inputs(nbin=8, nens=5)
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    assert su.shape == (8, 5)
    assert sv.shape == (8, 5)
    assert sw.shape == (8, 5)
```

- [ ] **Step 2: Run tests to verify they fail**

```
uv run pytest tests/test_shear.py::test_central_diff_shear_exact_gradient -x -q
```

Expected: `ImportError: cannot import name '_central_diff_shear'`

- [ ] **Step 3: Implement `_central_diff_shear`**

Add to `src/ladcp/solution/shear.py` (before `compute_shear`):

```python
def _central_diff_shear(
    u: NDArray,
    v: NDArray,
    w: NDArray,
    izm: NDArray,
    weight_mask: NDArray,
) -> tuple[NDArray, NDArray, NDArray]:
    """Compute stride-2 central-difference shear, weighted and NaN-padded.

    Replicates MATLAB diff2(): x[2:,:] - x[:-2,:] divided by izm[2:] - izm[:-2].
    First and last rows are NaN because no two-sided neighbour exists.
    weight_mask is 1.0 where valid, NaN where excluded.
    """
    du = u[2:, :] - u[:-2, :]     # (nbin-2, nens)
    dv = v[2:, :] - v[:-2, :]
    dw = w[2:, :] - w[:-2, :]
    dz = izm[2:, :] - izm[:-2, :]  # depth increment between skip-one neighbours

    shear_u = np.full_like(u, np.nan)
    shear_v = np.full_like(v, np.nan)
    shear_w = np.full_like(w, np.nan)

    shear_u[1:-1, :] = du / dz
    shear_v[1:-1, :] = dv / dz
    shear_w[1:-1, :] = dw / dz

    shear_u *= weight_mask
    shear_v *= weight_mask
    shear_w *= weight_mask

    return shear_u, shear_v, shear_w
```

- [ ] **Step 4: Run the new tests**

```
uv run pytest tests/test_shear.py -k "central_diff" -v
```

Expected: all four `test_central_diff_shear_*` tests PASS.

- [ ] **Step 5: Run the full suite to verify no regressions**

```
uv run pytest -x -q
```

Expected: all previously passing tests still pass; shear integration tests still fail with `NotImplementedError`.

- [ ] **Step 6: Commit**

```bash
git add src/ladcp/solution/shear.py tests/test_shear.py
git commit -m "feat: implement _central_diff_shear() — stride-2 central difference shear"
```

---

### Task 3: Implement `_bin_average_shear()` — Stage 2

For each output depth bin, gather all shear estimates within a ±dz window, reject outliers beyond `stdf` standard deviations from the median, and return the bin mean and std.

**Files:**
- Modify: `src/ladcp/solution/shear.py`
- Modify: `tests/test_shear.py`

**Interfaces:**
- Consumes: `shear_u, shear_v, shear_w` each `(nbin, nens)` from `_central_diff_shear`; `izm (nbin, nens)`; `z_bins (nz,)`; `dz float`; `stdf float`
- Produces: `_bin_average_shear(shear_u, shear_v, shear_w, izm, z_bins, dz, stdf) -> tuple[NDArray×6, NDArray]`
  - Returns `(usm, vsm, wsm, use, vse, wse, nn)` each `(nz,)` — mean shear, std shear, count

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_shear.py`:

```python
from ladcp.solution.shear import _bin_average_shear


def test_bin_average_recovers_known_shear():
    """Constant shear field: bin mean should equal the true shear value."""
    du_dz = 2e-3
    u, v, w, izm, weight = _constant_shear_inputs(
        nbin=10, nens=30, du_dz=du_dz, dv_dz=0.0, dw_dz=0.0,
        bin_spacing=10.0, first_bin_depth=10.0,
    )
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, vsm, wsm, use, vse, wse, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    # Interior bins (have enough samples) should recover du_dz
    populated = nn > 2
    assert populated.sum() > 0, "No bins had more than 2 samples"
    assert np.allclose(usm[populated], du_dz, atol=1e-10)


def test_bin_average_outlier_rejection():
    """Inject a single large outlier; after 2σ editing it must not influence the mean."""
    nbin, nens = 6, 40
    du_dz = 1e-3
    u, v, w, izm, weight = _constant_shear_inputs(
        nbin=nbin, nens=nens, du_dz=du_dz, dv_dz=0.0, dw_dz=0.0,
        bin_spacing=10.0, first_bin_depth=10.0,
    )
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    # Inject spike in bin 3, first ensemble — 100× larger than normal
    su[3, 0] = 0.1
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, _, _, _, _, _, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    populated = nn > 2
    # After 2σ editing the mean in every bin should still be ≈ du_dz
    assert np.allclose(usm[populated], du_dz, atol=1e-4)


def test_bin_average_few_samples_returns_nan():
    """Bins with ≤ 2 valid samples return NaN (not enough to estimate std)."""
    nbin, nens = 4, 2   # only 2 ensembles → ≤ 2 finite shear estimates per bin
    u, v, w, izm, weight = _constant_shear_inputs(nbin=nbin, nens=nens)
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, _, _, _, _, _, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    assert np.all(np.isnan(usm[nn <= 2]))


def test_bin_average_output_shapes():
    u, v, w, izm, weight = _constant_shear_inputs(nbin=8, nens=20)
    mask = np.ones_like(u)
    su, sv, sw = _central_diff_shear(u, v, w, izm, mask)
    dz = 10.0
    z_max = np.nanmax(izm)
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)
    usm, vsm, wsm, use, vse, wse, nn = _bin_average_shear(su, sv, sw, izm, z_bins, dz)
    nz = len(z_bins)
    for arr in (usm, vsm, wsm, use, vse, wse, nn):
        assert arr.shape == (nz,)
```

- [ ] **Step 2: Run tests to verify they fail**

```
uv run pytest tests/test_shear.py -k "bin_average" -x -q
```

Expected: `ImportError: cannot import name '_bin_average_shear'`

- [ ] **Step 3: Implement `_bin_average_shear`**

Add to `src/ladcp/solution/shear.py` (after `_central_diff_shear`, before `compute_shear`):

```python
def _bin_average_shear(
    shear_u: NDArray,
    shear_v: NDArray,
    shear_w: NDArray,
    izm: NDArray,
    z_bins: NDArray,
    dz: float,
    stdf: float = 2.0,
) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray, NDArray]:
    """Average shear estimates into depth bins with 2σ outlier editing.

    For each output bin centred at z_bins[k]:
      1. Collect all (bin, ensemble) pairs where depth is within dz of z_bins[k]+dz/2.
      2. Discard NaN estimates.
      3. If n ≤ 2: leave mean/std as NaN.
      4. Reject estimates more than stdf*std from the median.
      5. If ≥ 2 remain: store mean and std.

    Window condition replicates MATLAB: abs(izm - (z + dz/2)) <= dz.
    """
    nz = len(z_bins)
    iz_flat = izm.ravel()
    su_flat = shear_u.ravel()
    sv_flat = shear_v.ravel()
    sw_flat = shear_w.ravel()

    usm = np.full(nz, np.nan)
    vsm = np.full(nz, np.nan)
    wsm = np.full(nz, np.nan)
    use = np.full(nz, np.nan)
    vse = np.full(nz, np.nan)
    wse = np.full(nz, np.nan)
    nn = np.zeros(nz, dtype=np.intp)

    for k, center in enumerate(z_bins):
        in_window = np.abs(iz_flat - (center + dz / 2)) <= dz
        finite = in_window & np.isfinite(su_flat + sv_flat)
        su = su_flat[finite]
        sv = sv_flat[finite]
        sw = sw_flat[finite & np.isfinite(sw_flat)]
        n = len(su)
        nn[k] = n
        if n <= 2:
            continue

        for arr, mean_out, std_out in [
            (su, usm, use),
            (sv, vsm, vse),
        ]:
            med = np.median(arr)
            std = np.std(arr)
            keep = np.abs(arr - med) < stdf * std
            if keep.sum() > 1:
                mean_out[k] = np.mean(arr[keep])
                std_out[k] = np.std(arr[keep])

        # w may have fewer finite values — filter independently
        if len(sw) > 2:
            med_w = np.median(sw)
            std_w = np.std(sw)
            keep_w = np.abs(sw - med_w) < stdf * std_w
            if keep_w.sum() > 1:
                wsm[k] = np.mean(sw[keep_w])
                wse[k] = np.std(sw[keep_w])

    return usm, vsm, wsm, use, vse, wse, nn
```

- [ ] **Step 4: Run the bin-average tests**

```
uv run pytest tests/test_shear.py -k "bin_average" -v
```

Expected: all four `test_bin_average_*` tests PASS.

- [ ] **Step 5: Run the full suite**

```
uv run pytest -x -q
```

Expected: all previously passing tests still pass.

- [ ] **Step 6: Commit**

```bash
git add src/ladcp/solution/shear.py tests/test_shear.py
git commit -m "feat: implement _bin_average_shear() — depth-bin averaging with 2σ outlier editing"
```

---

### Task 4: Implement `_integrate_shear()` and wire `compute_shear()` end-to-end

Bottom-up cumulative sum integration produces the relative velocity profile, which is then zero-meaned. This completes the three-stage pipeline and makes all unit tests pass.

**Files:**
- Modify: `src/ladcp/solution/shear.py`
- Modify: `tests/test_shear.py`

**Interfaces:**
- Consumes: `usm, vsm, wsm` each `(nz,)` from `_bin_average_shear`; `dz float`
- Produces: `_integrate_shear(usm, vsm, wsm, dz) -> tuple[NDArray, NDArray, NDArray]`
  - Returns `(u_rel, v_rel, w_rel)` each `(nz,)`, zero-mean

- [ ] **Step 1: Write the failing integration tests**

Add to `tests/test_shear.py`:

```python
from ladcp.solution.shear import _integrate_shear


def test_integrate_shear_constant_shear_linear_profile():
    """Constant shear s → integrated velocity grows linearly from deepest bin."""
    nz = 10
    dz = 10.0
    shear = 1e-3  # s^-1, uniform
    usm = np.full(nz, shear)
    vsm = np.zeros(nz)
    wsm = np.zeros(nz)
    ur, vr, wr = _integrate_shear(usm, vsm, wsm, dz)
    # After zero-mean, profile should be linear (slope = shear * dz)
    diff = np.diff(ur)
    assert np.allclose(diff, -shear * dz, atol=1e-12), \
        "Integrated profile should decrease by shear*dz per bin (bottom-up)"


def test_integrate_shear_zero_mean():
    """Integrated profiles must have exactly zero mean."""
    nz = 15
    dz = 8.0
    rng = np.random.default_rng(42)
    usm = rng.standard_normal(nz) * 1e-3
    vsm = rng.standard_normal(nz) * 1e-3
    wsm = np.zeros(nz)
    ur, vr, wr = _integrate_shear(usm, vsm, wsm, dz)
    assert np.mean(ur) == pytest.approx(0.0, abs=1e-12)
    assert np.mean(vr) == pytest.approx(0.0, abs=1e-12)


def test_integrate_shear_nan_gaps_filled_with_zero():
    """NaN shear bins are treated as zero shear (no contribution to integral)."""
    nz = 6
    dz = 10.0
    usm = np.array([np.nan, 1e-3, np.nan, 1e-3, np.nan, np.nan])
    vsm = np.zeros(nz)
    wsm = np.zeros(nz)
    ur, vr, wr = _integrate_shear(usm, vsm, wsm, dz)
    assert np.all(np.isfinite(ur)), "NaN gaps must not propagate to integrated profile"
    assert np.mean(ur) == pytest.approx(0.0, abs=1e-12)
```

- [ ] **Step 2: Run tests to confirm they fail**

```
uv run pytest tests/test_shear.py -k "integrate" -x -q
```

Expected: `ImportError: cannot import name '_integrate_shear'`

- [ ] **Step 3: Implement `_integrate_shear`**

Add to `src/ladcp/solution/shear.py` (after `_bin_average_shear`, before `compute_shear`):

```python
def _integrate_shear(
    usm: NDArray,
    vsm: NDArray,
    wsm: NDArray,
    dz: float,
) -> tuple[NDArray, NDArray, NDArray]:
    """Integrate shear from bottom up to produce zero-mean relative velocities.

    NaN bins are replaced with 0 (no shear contribution) before integrating.
    Replicates MATLAB: flipud(cumsum(flipud(usm))) * dz, then subtract mean.
    """
    u = np.where(np.isnan(usm), 0.0, usm)
    v = np.where(np.isnan(vsm), 0.0, vsm)
    w = np.where(np.isnan(wsm), 0.0, wsm)

    ur = np.flipud(np.cumsum(np.flipud(u))) * dz
    vr = np.flipud(np.cumsum(np.flipud(v))) * dz
    wr = np.flipud(np.cumsum(np.flipud(w))) * dz

    ur -= np.mean(ur)
    vr -= np.mean(vr)
    wr -= np.mean(wr)

    return ur, vr, wr
```

- [ ] **Step 4: Wire `compute_shear` end-to-end**

Replace the `raise NotImplementedError` body inside `compute_shear` with:

```python
    weight_mask = np.where(weight > weight_min, 1.0, np.nan)

    shear_u, shear_v, shear_w = _central_diff_shear(u, v, w, izm, weight_mask)

    z_max = float(np.nanmax(izm))
    z_bins = np.arange(dz / 2, z_max + dz / 2, dz)

    usm, vsm, wsm, use, vse, wse, nn = _bin_average_shear(
        shear_u, shear_v, shear_w, izm, z_bins, dz, stdf
    )

    ur, vr, wr = _integrate_shear(usm, vsm, wsm, dz)

    return ShearProfile(
        z=z_bins,
        u_shear=usm, v_shear=vsm, w_shear=wsm,
        u_shear_err=use, v_shear_err=vse, w_shear_err=wse,
        n=nn,
        u_rel=ur, v_rel=vr, w_rel=wr,
    )
```

- [ ] **Step 5: Run all shear unit tests**

```
uv run pytest tests/test_shear.py -v
```

Expected: all shear unit tests PASS (including the end-to-end tests `test_compute_shear_*` from Task 1).

- [ ] **Step 6: Run the full suite**

```
uv run pytest -x -q
```

Expected: 71 previously passing tests still pass; new shear tests also pass.

- [ ] **Step 7: Commit**

```bash
git add src/ladcp/solution/shear.py tests/test_shear.py
git commit -m "feat: implement _integrate_shear(); wire compute_shear() end-to-end"
```

---

### Task 5: Integration test against P16N cast 003 reference

Validate that `compute_shear()` produces a physically plausible profile when driven by real PD0 data and CTD-derived bin depths. Compare depth range and velocity magnitude against the LDEO reference `003.nc`.

**Files:**
- Create: `tests/integration/test_shear_p16n_cast003.py`

**Interfaces:**
- Consumes: `compute_shear` (from Task 4), `load_rdi`, `load_ctd`, `assign_bin_depths` (all previously implemented)
- Produces: integration test module (no new exports)

Note: The reference `003.nc` was produced by the full LDEO pipeline including weight functions and data editing that we have not yet implemented. Tolerances therefore check shape, depth coverage, and order-of-magnitude agreement rather than exact reproduction.

- [ ] **Step 1: Write the integration test**

Create `tests/integration/test_shear_p16n_cast003.py`:

```python
"""Integration tests: shear solution against P16N cast 003 LDEO reference.

Requires TEST_DATA_DIR env var pointing to a directory containing:
  2015_P16N/003DL000.000   — Downlooker PD0 binary
  2015_P16N/003_01.cnv     — CTD time-series (binary SBE)
  2015_P16N/003.nc         — LDEO_IX processed reference output
"""
import os
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.ctd import assign_bin_depths, load_ctd
from ladcp.ingestion.rdi import load_rdi
from ladcp.solution.shear import ShearProfile, compute_shear


@pytest.fixture
def test_data_dir() -> Path:
    path = Path(os.environ.get("TEST_DATA_DIR", "test_data"))
    if not path.exists():
        pytest.skip("TEST_DATA_DIR not populated — see test_data/sources.md")
    return path


@pytest.fixture
def dl_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003DL000.000"
    if not p.exists():
        pytest.skip(f"DL PD0 file not found: {p}")
    return p


@pytest.fixture
def cnv_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003_01.cnv"
    if not p.exists():
        pytest.skip(f"CTD file not found: {p}")
    return p


@pytest.fixture
def ref_path(test_data_dir: Path) -> Path:
    p = test_data_dir / "2015_P16N" / "003.nc"
    if not p.exists():
        pytest.skip(f"Reference NetCDF not found: {p}")
    return p


@pytest.fixture
def shear_result(dl_path: Path, cnv_path: Path) -> ShearProfile:
    rdi = load_rdi(dl_path)
    ctd = load_ctd(cnv_path)
    _, izm = assign_bin_depths(rdi, ctd, looker="down")
    weight = np.ones((rdi.nbin, rdi.nens), dtype=np.float64)
    return compute_shear(rdi.u, rdi.v, rdi.w, izm, weight, dz=10.0)


@pytest.mark.integration
def test_shear_profile_is_shear_profile(shear_result: ShearProfile):
    assert isinstance(shear_result, ShearProfile)


@pytest.mark.integration
def test_shear_depth_axis_starts_at_5m(shear_result: ShearProfile):
    assert shear_result.z[0] == pytest.approx(5.0, abs=0.1)


@pytest.mark.integration
def test_shear_depth_axis_covers_cast_depth(shear_result: ShearProfile, ref_path: Path):
    ds = netCDF4.Dataset(ref_path)
    ref_z_max = float(np.max(ds.variables["z"][:]))
    ds.close()
    # Our profile must reach at least 80 % of the reference maximum depth
    assert shear_result.z[-1] >= 0.8 * ref_z_max, (
        f"z_max={shear_result.z[-1]:.0f} m vs ref {ref_z_max:.0f} m"
    )


@pytest.mark.integration
def test_integrated_profile_zero_mean(shear_result: ShearProfile):
    assert np.mean(shear_result.u_rel) == pytest.approx(0.0, abs=1e-12)
    assert np.mean(shear_result.v_rel) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.integration
def test_shear_magnitude_plausible(shear_result: ShearProfile, ref_path: Path):
    """Integrated velocity should be the same order of magnitude as the reference."""
    ds = netCDF4.Dataset(ref_path)
    ref_u = np.array(ds.variables["u_shear_method"][:])
    ds.close()
    ref_rms = float(np.sqrt(np.nanmean(ref_u**2)))
    our_rms = float(np.sqrt(np.nanmean(shear_result.u_rel**2)))
    # Without full weight editing our RMS may differ by up to 5×, but same order
    assert our_rms > ref_rms / 10.0, f"Our u_rel RMS {our_rms:.4f} is too small vs ref {ref_rms:.4f}"
    assert our_rms < ref_rms * 10.0, f"Our u_rel RMS {our_rms:.4f} is too large vs ref {ref_rms:.4f}"


@pytest.mark.integration
def test_shear_profile_n_populated(shear_result: ShearProfile):
    """More than half the depth bins should have valid shear estimates."""
    assert (shear_result.n > 2).sum() > len(shear_result.n) // 2
```

- [ ] **Step 2: Run to confirm skip when TEST_DATA_DIR absent**

```
uv run pytest tests/integration/test_shear_p16n_cast003.py -v
```

Expected: all 6 tests show `SKIPPED (TEST_DATA_DIR not populated)`.

- [ ] **Step 3: If TEST_DATA_DIR is set, run integration tests**

```
TEST_DATA_DIR=test_data uv run pytest tests/integration/test_shear_p16n_cast003.py -v
```

Expected: all 6 tests PASS (or skip if files absent).

- [ ] **Step 4: Run full suite**

```
uv run pytest -x -q
```

Expected: all previously passing tests still pass; new integration tests skip cleanly.

- [ ] **Step 5: Commit**

```bash
git add tests/integration/test_shear_p16n_cast003.py
git commit -m "test: add integration tests for shear solution against P16N cast 003 reference"
```

---

## Self-Review

### Spec coverage

| Requirement | Task |
|---|---|
| Central-difference shear from Earth-frame velocities + bin depths | Task 2 (`_central_diff_shear`) |
| Depth-bin averaging with 2σ outlier editing | Task 3 (`_bin_average_shear`) |
| Bottom-up cumsum integration → zero-mean profile | Task 4 (`_integrate_shear`) |
| `ShearProfile` dataclass with all output fields | Task 1 |
| `compute_shear()` public API replaces too-narrow stub | Tasks 1 + 4 |
| Validation against `test_data/2015_P16N/003.nc` | Task 5 |
| Integration tests gated on `TEST_DATA_DIR` | Tasks 5 (fixture skips) |
| All new tests pass alongside existing 71 | All tasks |

### Type consistency check

| Name | Defined in | Used in |
|---|---|---|
| `ShearProfile` | Task 1 | Tasks 4, 5 |
| `_central_diff_shear(u, v, w, izm, weight_mask) -> tuple[NDArray, NDArray, NDArray]` | Task 2 | Task 3 tests, Task 4 |
| `_bin_average_shear(...) -> tuple[NDArray×6, NDArray]` | Task 3 | Task 4 |
| `_integrate_shear(usm, vsm, wsm, dz) -> tuple[NDArray, NDArray, NDArray]` | Task 4 | Task 4 (compute_shear body) |
| `compute_shear(u, v, w, izm, weight, dz, *, stdf, weight_min) -> ShearProfile` | Task 1 (stub), Task 4 (body) | Task 5 |

All cross-references are consistent.

### Placeholder scan

No TBD, TODO, or vague instructions detected. All steps contain complete code.
