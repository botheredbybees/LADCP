# LADCP Inverse Solver Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the LDEO LADCP velocity inverse solver in Python, replicating `prepinv.m` (super-ensemble formation) and `getinv.m` (constrained least-squares inversion) to produce `u(z)`, `v(z)` profiles matching the P16N cast 003 reference output within ±0.05 m/s RMS.

**Architecture:** Two-stage pipeline. Stage 1 (`prepare_superensembles`) depth-window-averages raw ADCP ensembles into super-ensembles (fewer, higher-SNR time steps). Stage 2 (`compute_inverse`) builds two sparse matrices — `A_ocean` mapping each measurement to a depth bin and `A_ctd` mapping it to a time bin — and solves `d = [A_ocean | A_ctd] * [u_ocean; -u_CTD]` with optional constraints (bottom-track, GPS/barotropic, SADCP, smoothness) appended as additional rows. U and V components share the same A matrix; MATLAB's complex-number trick is avoided by two identical solves.

**Tech Stack:** Python 3.11, NumPy, SciPy (`scipy.sparse`, `scipy.linalg.lstsq`), xarray (integration test I/O only).

## Global Constraints

- Depth convention: **negative below surface** (`izm[b, e] ≤ 0`). `izv = -izm` gives positive depths used for matrix column indexing.
- `izd` / `izu`: 0-indexed integer arrays of downlooker / uplooker bin rows in the first axis of `u/v/w/weight/izm`.
- `np.std(ddof=1)` everywhere — MATLAB `std()` uses `ddof=1`.
- `np.nanmedian` / `np.nanmean` for all averaging; never call `np.mean` on potentially-NaN arrays.
- Sparse matrices: build as `lil_matrix`, convert to `csr_matrix` before arithmetic.
- After solving `m = [u_ocean; u_ctd_neg]`: `u_CTD = -m[n_zbins:]` (sign flip — see math note in Task 5).
- Integration test in `tests/integration/test_inverse_p16n_cast003.py`, gated on `TEST_DATA_DIR` env var (same pattern as shear integration test).
- All source in `src/ladcp/solution/inverse.py`; tests in `tests/test_inverse.py`.
- One commit per completed task.

---

### Task 1: EnsembleData + SuperEnsemble dataclasses + depth-window averaging

**Files:**
- Create: `src/ladcp/solution/inverse.py`
- Create: `tests/test_inverse.py`

**Interfaces:**
- Produces:
  - `EnsembleData` — Earth-frame velocity arrays + metadata; input to `prepare_superensembles()`
  - `SuperEnsemble` — averaged super-ensembles; input to `compute_inverse()`
  - `prepare_superensembles(ens: EnsembleData, *, dz: float | None = None) -> SuperEnsemble`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_inverse.py
import numpy as np
import pytest
from ladcp.solution.inverse import EnsembleData, SuperEnsemble, prepare_superensembles


def _make_ens(*, n_bins: int = 5, n_ens: int = 40,
               u_val: float = 0.5, v_val: float = 0.1) -> EnsembleData:
    """Synthetic EnsembleData where CTD descends 5 m per ensemble."""
    z = np.linspace(-20.0, -20.0 - 5.0 * n_ens, n_ens)
    izm = np.zeros((n_bins, n_ens))
    for b in range(n_bins):
        izm[b] = z - (b + 0.5) * 8.0
    return EnsembleData(
        u=np.full((n_bins, n_ens), u_val),
        v=np.full((n_bins, n_ens), v_val),
        w=np.zeros((n_bins, n_ens)),
        weight=np.ones((n_bins, n_ens)),
        izm=izm,
        z=z,
        time_jul=np.linspace(2457000.0, 2457001.0, n_ens),
        bvel=np.full((n_ens, 3), np.nan),
        bvels=np.full((n_ens, 3), np.nan),
        hbot=np.full(n_ens, np.nan),
        izd=np.arange(n_bins),
        izu=np.array([], dtype=int),
        slat=np.full(n_ens, np.nan),
        slon=np.full(n_ens, np.nan),
    )


def test_prepare_superensembles_reduces_n_ens():
    """Super-ensemble averaging must collapse multiple raw ensembles."""
    ens = _make_ens(n_ens=40)
    se = prepare_superensembles(ens, dz=20.0)
    assert se.ru.shape[1] < 40


def test_prepare_superensembles_preserves_mean_velocity():
    """Constant velocity field should survive depth-window averaging unchanged."""
    ens = _make_ens(u_val=0.5, v_val=0.1)
    se = prepare_superensembles(ens, dz=20.0)
    assert np.allclose(np.nanmean(se.ru), 0.5, atol=0.02)
    assert np.allclose(np.nanmean(se.rv), 0.1, atol=0.02)


def test_prepare_superensembles_default_dz():
    """Default dz (inferred from izm spacing) must give the same shape as explicit dz."""
    ens = _make_ens(n_ens=40)
    se_default = prepare_superensembles(ens)
    # default dz = median(|diff(izm[0])|) ≈ 5 m
    se_explicit = prepare_superensembles(ens, dz=5.0)
    assert se_default.ru.shape == se_explicit.ru.shape


def test_prepare_superensembles_bvel_preserved():
    """Bottom-track velocity should be averaged into super-ensembles."""
    ens = _make_ens(n_ens=40)
    ens.bvel[:] = [0.1, -0.2, -1.0]
    se = prepare_superensembles(ens, dz=20.0)
    assert np.allclose(np.nanmean(se.bvel[0]), 0.1, atol=0.02)
```

- [ ] **Step 2: Run tests to confirm import error**

```
pytest tests/test_inverse.py -v
```
Expected: `ModuleNotFoundError` (file not created yet)

- [ ] **Step 3: Implement `inverse.py` with dataclasses and `prepare_superensembles`**

```python
# src/ladcp/solution/inverse.py
"""LADCP inverse velocity solver.

Replicates prepinv.m (super-ensemble formation) and getinv.m (constrained
least-squares inversion) from the LDEO_IX MATLAB reference implementation.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field


@dataclass
class EnsembleData:
    """Earth-frame ADCP data aligned with CTD depths, input to prepare_superensembles().

    Depth convention: negative = below surface. izd / izu are 0-indexed row indices
    into the first dimension of u/v/w/weight/izm.
    """
    u: np.ndarray         # (n_bins, n_ens) eastward velocity m/s
    v: np.ndarray         # (n_bins, n_ens) northward velocity m/s
    w: np.ndarray         # (n_bins, n_ens) vertical velocity m/s
    weight: np.ndarray    # (n_bins, n_ens) quality weight 0–1
    izm: np.ndarray       # (n_bins, n_ens) bin depth m (≤ 0)
    z: np.ndarray         # (n_ens,) CTD depth m (≤ 0)
    time_jul: np.ndarray  # (n_ens,) Julian day
    bvel: np.ndarray      # (n_ens, 3) bottom track u, v, w m/s (NaN = absent)
    bvels: np.ndarray     # (n_ens, 3) bottom track std m/s
    hbot: np.ndarray      # (n_ens,) height above bottom m (NaN = absent)
    izd: np.ndarray       # (n_dl_bins,) downlooker bin row indices (int)
    izu: np.ndarray       # (n_ul_bins,) uplooker bin row indices (int)
    slat: np.ndarray      # (n_ens,) latitude (NaN if unavailable)
    slon: np.ndarray      # (n_ens,) longitude (NaN if unavailable)


@dataclass
class SuperEnsemble:
    """Depth-window-averaged ADCP data, output of prepare_superensembles().

    Shape conventions: 2-D fields are (n_bins, n_se); 1-D fields are (n_se,).
    bvel / bvels are (3, n_se) — transposed from EnsembleData for column access.
    """
    ru: np.ndarray        # (n_bins, n_se) eastward velocity m/s
    rv: np.ndarray        # (n_bins, n_se) northward velocity m/s
    rw: np.ndarray        # (n_bins, n_se) vertical velocity m/s
    ruvs: np.ndarray      # (n_bins, n_se) combined U+V std (velocity uncertainty)
    weight: np.ndarray    # (n_bins, n_se) quality weight 0–1
    izm: np.ndarray       # (n_bins, n_se) bin depth m (≤ 0)
    z: np.ndarray         # (n_se,) CTD depth m (≤ 0)
    dt: np.ndarray        # (n_se,) time interval s
    time_jul: np.ndarray  # (n_se,) Julian day
    bvel: np.ndarray      # (3, n_se) bottom track u, v, w m/s
    bvels: np.ndarray     # (3, n_se) bottom track std m/s
    hbot: np.ndarray      # (n_se,) height above bottom m
    slat: np.ndarray      # (n_se,) latitude
    slon: np.ndarray      # (n_se,) longitude
    izd: np.ndarray       # downlooker bin row indices (unchanged from input)
    izu: np.ndarray       # uplooker bin row indices (unchanged from input)


def _window_boundaries(depth0: np.ndarray, avdz: float) -> list[np.ndarray]:
    """Partition ensemble indices into depth-triggered windows.

    Replicates the while-loop in prepinv.m lines 499–605. Each window spans
    consecutive ensembles until |depth0[t] - depth0[window_start]| > avdz.
    """
    n = len(depth0)
    windows: list[np.ndarray] = []
    ilast = 0
    while ilast < n:
        if ilast + 1 >= n:
            windows.append(np.array([ilast]))
            break
        depth_change = np.abs(depth0[ilast + 1:] - depth0[ilast])
        ii = np.where(depth_change > avdz)[0]
        end = int(ii[0]) + 1 if len(ii) > 0 else n - ilast
        i1 = np.arange(ilast, min(ilast + end, n))
        windows.append(i1)
        ilast = int(i1[-1]) + 1
    return windows


def prepare_superensembles(
    ens: EnsembleData,
    *,
    dz: float | None = None,
) -> SuperEnsemble:
    """Form super-ensembles by depth-window averaging (replicates prepinv.m).

    Parameters
    ----------
    ens:
        Earth-frame ADCP data from beam2earth + assign_bin_depths.
    dz:
        Depth window size m. Defaults to ``median(|diff(izm[0, :])|)``,
        matching MATLAB's ``medianan(abs(diff(d.izm(:,1))))``.
    """
    if dz is None:
        dz = float(np.nanmedian(np.abs(np.diff(ens.izm[0]))))

    # Reference bins: 2nd and 3rd downlooker bins (0-indexed: offset 1, 2 from min(izd))
    # Replicates: izr = [min(d.izd)+1, min(d.izd)+2]  (MATLAB 1-indexed → same offsets)
    izr: np.ndarray
    if len(ens.izd) > 2:
        base = int(ens.izd.min())
        izr = np.array([base + 1, base + 2], dtype=int)
    else:
        izr = ens.izd.copy()
    if len(ens.izu) > 2:
        ul_top = int(ens.izu.max())
        izr = np.concatenate([izr, [ul_top - 1, ul_top - 2]])

    windows = _window_boundaries(ens.z, dz)
    n_se = len(windows)
    n_bins = ens.u.shape[0]

    ru = np.full((n_bins, n_se), np.nan)
    rv = np.full((n_bins, n_se), np.nan)
    rw = np.full((n_bins, n_se), np.nan)
    ruvs = np.full((n_bins, n_se), np.nan)
    weight_se = np.full((n_bins, n_se), np.nan)
    izm_se = np.full((n_bins, n_se), np.nan)
    z_se = np.full(n_se, np.nan)
    time_se = np.full(n_se, np.nan)
    bvel_se = np.full((3, n_se), np.nan)
    bvels_se = np.full((3, n_se), np.nan)
    hbot_se = np.full(n_se, np.nan)
    slat_se = np.full(n_se, np.nan)
    slon_se = np.full(n_se, np.nan)

    for im, i1 in enumerate(windows):
        u_win = ens.u[:, i1]   # (n_bins, n_win) — no weight masking here
        v_win = ens.v[:, i1]   # (MATLAB uses w=d.weight*0+1 = all ones)
        w_win = ens.w[:, i1]

        # Per-ensemble reference velocity: median of reference bins
        # MATLAB: ur = medianan(d.ru(izr, i1)) — plain median, no weight applied
        izr_valid = izr[izr < n_bins]
        ur_t = np.nanmedian(u_win[izr_valid], axis=0)  # (n_win,)
        vr_t = np.nanmedian(v_win[izr_valid], axis=0)

        # Remove per-ensemble reference, take median, add back mean reference
        # MATLAB: di.ru(:,im) = medianan(d.ru - ur_broadcast)' + mean(ur)
        u_deref = u_win - ur_t[np.newaxis, :]
        v_deref = v_win - vr_t[np.newaxis, :]

        ru[:, im] = np.nanmedian(u_deref, axis=1) + np.nanmean(ur_t)
        rv[:, im] = np.nanmedian(v_deref, axis=1) + np.nanmean(vr_t)
        rw[:, im] = np.nanmedian(w_win, axis=1)

        # Velocity uncertainty: combined U+V std over window
        ruvs[:, im] = np.sqrt(
            np.nanstd(u_win, axis=1, ddof=1) ** 2
            + np.nanstd(v_win, axis=1, ddof=1) ** 2
        )

        weight_se[:, im] = np.nanmean(ens.weight[:, i1], axis=1)
        izm_se[:, im] = np.nanmean(ens.izm[:, i1], axis=1)
        z_se[im] = np.nanmean(ens.z[i1])
        time_se[im] = np.nanmean(ens.time_jul[i1])

        bvel_se[:, im] = np.nanmean(ens.bvel[i1], axis=0)
        bvels_se[:, im] = np.nanstd(ens.bvel[i1], axis=0, ddof=1)
        hbot_se[im] = np.nanmean(ens.hbot[i1])
        slat_se[im] = np.nanmedian(ens.slat[i1])
        slon_se[im] = np.nanmedian(ens.slon[i1])

    # Time interval between super-ensembles (seconds); mirror edge values
    dt_mid = np.diff(time_se) * 24.0 * 3600.0
    dt = np.concatenate([[dt_mid[0]], (dt_mid[:-1] + dt_mid[1:]) / 2.0, [dt_mid[-1]]])

    # Floor std at single_ping_err (≈0.01 m/s); propagate NaN from weight
    # Replicates prepinv.m's superens_std_min logic
    single_ping_err = 0.01
    zero_mask = ruvs == 0
    weight_se[zero_mask] = np.nan
    ruvs[ruvs < single_ping_err] = single_ping_err
    ruvs = ruvs + weight_se * 0  # NaN-propagate

    return SuperEnsemble(
        ru=ru, rv=rv, rw=rw, ruvs=ruvs, weight=weight_se, izm=izm_se,
        z=z_se, dt=dt, time_jul=time_se,
        bvel=bvel_se, bvels=bvels_se, hbot=hbot_se,
        slat=slat_se, slon=slon_se,
        izd=ens.izd, izu=ens.izu,
    )
```

- [ ] **Step 4: Run tests**

```
pytest tests/test_inverse.py -v
```
Expected: all 4 tests PASS

- [ ] **Step 5: Commit**

```
git add src/ladcp/solution/inverse.py tests/test_inverse.py
git commit -m "feat: EnsembleData + SuperEnsemble dataclasses + depth-window averaging"
```

---

### Task 2: Flatten data + build sparse observation matrices + apply data weights

The observation equation is:
```
u_ADCP(b, e) = u_ocean(z_b_e) - u_CTD(e)
```
Written as a sparse linear system `d = [A_ocean | A_ctd] * [u_ocean; u_ctd_neg]`
where `u_ctd_neg = -u_CTD`. After solving, CTD velocity = `-m[n_zbins:]`.

`A_ocean[k, j] = 1` if observation k maps to depth bin j (j = round(izv[k] / dz)).
`A_ctd[k, e] = 1` if observation k is from super-ensemble e.

**Files:**
- Modify: `src/ladcp/solution/inverse.py` (add `_flatten_obs`, `_build_obs_matrix`, `_build_ctd_matrix`, `_apply_weights`)
- Modify: `tests/test_inverse.py`

**Interfaces:**
- Consumes: `SuperEnsemble` from Task 1
- Produces (all module-private):
  - `_flatten_obs(se, velerr, weightmin) -> tuple[d_u, d_v, izv, jprof, wm]`
  - `_build_obs_matrix(izv, dz) -> csr_matrix` shape (n_obs, n_zbins)
  - `_build_ctd_matrix(jprof, n_se) -> csr_matrix` shape (n_obs, n_se)
  - `_apply_weights(A_ocean, A_ctd, d, wm) -> tuple[A_ocean_w, A_ctd_w, d_w, idx_down, idx_up]`

- [ ] **Step 1: Add tests for matrix construction**

Append to `tests/test_inverse.py`:

```python
from ladcp.solution.inverse import (
    _flatten_obs, _build_obs_matrix, _build_ctd_matrix, _apply_weights,
)
import scipy.sparse


def test_build_obs_matrix_shape():
    """obs matrix rows = n_obs, cols = number of unique depth bins."""
    izv = np.array([10.0, 20.0, 30.0, 10.0])  # 3 unique bins
    A = _build_obs_matrix(izv, dz=10.0)
    assert A.shape == (4, 3)


def test_build_obs_matrix_one_nonzero_per_row():
    """Each observation maps to exactly one depth bin."""
    izv = np.array([5.0, 15.0, 25.0, 35.0])
    A = _build_obs_matrix(izv, dz=10.0)
    assert np.allclose(np.asarray(A.sum(axis=1)).ravel(), 1.0)


def test_build_ctd_matrix_shape():
    """ctd matrix rows = n_obs, cols = n_se."""
    jprof = np.array([0, 0, 1, 1, 2, 2], dtype=int)
    A = _build_ctd_matrix(jprof, n_se=3)
    assert A.shape == (6, 3)
    assert np.allclose(np.asarray(A.sum(axis=1)).ravel(), 1.0)


def test_flatten_obs_removes_nan_weight():
    """Observations with NaN or zero weight must be excluded."""
    ens = _make_ens(n_bins=3, n_ens=20)
    se = prepare_superensembles(ens, dz=10.0)
    # Force some NaN weights
    se.weight[:, :2] = np.nan
    d_u, d_v, izv, jprof, wm = _flatten_obs(se, velerr=0.05, weightmin=0.05)
    assert np.all(np.isfinite(d_u))
    assert np.all(np.isfinite(wm))
    assert np.all(wm >= 0.05)
```

- [ ] **Step 2: Run to confirm failures**

```
pytest tests/test_inverse.py::test_build_obs_matrix_shape -v
```
Expected: `ImportError`

- [ ] **Step 3: Implement the four functions**

Add to `src/ladcp/solution/inverse.py`:

```python
import scipy.sparse as sp
from scipy.sparse import lil_matrix, csr_matrix


def _flatten_obs(
    se: SuperEnsemble,
    velerr: float = 0.05,
    weightmin: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Flatten super-ensemble data into observation vectors for matrix construction.

    Column-major (Fortran) order matches MATLAB reshape(x, nbin*nt, 1).

    Returns
    -------
    d_u, d_v : (n_obs,) observed velocities m/s
    izv      : (n_obs,) positive bin depths m  (= -izm, column-major flattened)
    jprof    : (n_obs,) super-ensemble index 0..n_se-1
    wm       : (n_obs,) data weights (velerr / ruvs, NaN-propagated)
    """
    n_bins, n_se = se.izm.shape

    # Observation depth: positive (izv = -izm), column-major flatten
    izv_full = (-se.izm).ravel(order="F")   # (n_bins * n_se,)

    # Profile index per observation: ensemble 0 repeated n_bins times, etc.
    jprof_full = np.repeat(np.arange(n_se), n_bins)

    d_u_full = se.ru.ravel(order="F")
    d_v_full = se.rv.ravel(order="F")

    # Data weight: velerr / std  (std-based weighting, MATLAB std_weight=1)
    # NaN in weight propagates to wm, excluding those observations
    wm_full = velerr / se.ruvs + se.weight * 0  # NaN-propagate
    wm_full = wm_full.ravel(order="F")

    # Keep only valid, well-weighted observations
    valid = (
        np.isfinite(d_u_full)
        & np.isfinite(d_v_full)
        & np.isfinite(wm_full)
        & (wm_full >= weightmin)
        & (izv_full > 0)
    )

    return (
        d_u_full[valid],
        d_v_full[valid],
        izv_full[valid],
        jprof_full[valid],
        wm_full[valid],
    )


def _build_obs_matrix(izv: np.ndarray, dz: float) -> csr_matrix:
    """Build A_ocean: maps each observation to a depth bin.

    Column j = round(izv[k] / dz) - 1 (0-indexed).
    Replicates lainseta(izv, dz) from getinv.m.
    """
    n_obs = len(izv)
    j = np.round(izv / dz).astype(int) - 1  # 0-indexed depth bin
    j = np.clip(j, 0, None)
    n_zbins = int(j.max()) + 1
    i = np.arange(n_obs)
    return csr_matrix((np.ones(n_obs), (i, j)), shape=(n_obs, n_zbins))


def _build_ctd_matrix(jprof: np.ndarray, n_se: int) -> csr_matrix:
    """Build A_ctd: maps each observation to its super-ensemble (time bin).

    Replicates lainseta(jprof, 1) from getinv.m.
    """
    n_obs = len(jprof)
    i = np.arange(n_obs)
    j = jprof.astype(int)
    return csr_matrix((np.ones(n_obs), (i, j)), shape=(n_obs, n_se))


def _apply_weights(
    A_ocean: csr_matrix,
    A_ctd: csr_matrix,
    d: np.ndarray,
    wm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Apply data weights to observation system; return dense arrays + split indices.

    Replicates lainweig() from getinv.m. Returns dense arrays (not sparse) because
    constraints are appended row-by-row in subsequent tasks.

    Returns
    -------
    A_o, A_c : dense float64 arrays with weights applied
    d_w      : weighted observation vector
    idx_down : row indices belonging to the downcast
    idx_up   : row indices belonging to the upcast
    """
    A_o = A_ocean.toarray() * wm[:, np.newaxis]
    A_c = A_ctd.toarray() * wm[:, np.newaxis]
    d_w = d * wm

    # Split down/up cast at the deepest super-ensemble
    # Column of A_ocean with the highest depth index = deepest depth bin
    col_sums = A_ocean.sum(axis=0).A1
    deepest_col = int(np.argmax(col_sums > 0) + col_sums[col_sums > 0].argmax())
    # Row belonging to deepest obs (proxy for cast bottom)
    rows_at_bottom = np.where(A_ocean.getcol(deepest_col).toarray().ravel() > 0)[0]
    if len(rows_at_bottom) > 0:
        split = int(np.median(rows_at_bottom))
    else:
        split = len(d) // 2
    idx_down = np.arange(0, split + 1)
    idx_up = np.arange(split, len(d))

    return A_o, A_c, d_w, idx_down, idx_up
```

- [ ] **Step 4: Run tests**

```
pytest tests/test_inverse.py -v
```
Expected: all tests PASS (including the 4 new ones)

- [ ] **Step 5: Commit**

```
git add src/ladcp/solution/inverse.py tests/test_inverse.py
git commit -m "feat: observation matrix construction and data weight application"
```

---

### Task 3: Smoothness + zero-mean fallback constraints

Both constraints append rows to `[A_ocean | A_ctd | d]`.

**Smoothness** (`lainsmoo`): adds curvature-penalty rows `[-1, 2, -1] * weight` for each interior column of A. Penalizes large second derivatives in the velocity profile. Applied to both A_ocean (smooth ocean velocity) and A_ctd (smooth CTD velocity) independently.

**Zero-mean** (`lainocean`): when no external velocity reference exists (no bottom track, no GPS, no SADCP), constrains `mean(u_ocean) = 0`. Adds one row `[1/n_zbins, 1/n_zbins, ...] * scale` to A_ocean with d=0.

**Files:**
- Modify: `src/ladcp/solution/inverse.py`
- Modify: `tests/test_inverse.py`

**Interfaces:**
- Produces (module-private):
  - `_add_smoothness(A_ocean, A_ctd, d, smoofac) -> (A_ocean, A_ctd, d)`
  - `_add_zero_mean(A_ocean, A_ctd, d) -> (A_ocean, A_ctd, d)`

- [ ] **Step 1: Add tests**

Append to `tests/test_inverse.py`:

```python
from ladcp.solution.inverse import _add_smoothness, _add_zero_mean


def test_add_smoothness_increases_rows():
    """Smoothness adds curvature rows for interior columns."""
    n_obs, n_zbins, n_se = 20, 8, 5
    A_o = np.random.rand(n_obs, n_zbins)
    A_c = np.random.rand(n_obs, n_se)
    d = np.random.rand(n_obs)
    A_o2, A_c2, d2 = _add_smoothness(A_o, A_c, d, smoofac=1.0)
    # Must add at least n_zbins - 2 curvature rows (interior bins)
    assert A_o2.shape[0] > n_obs
    assert A_c2.shape[0] == A_o2.shape[0]
    assert len(d2) == A_o2.shape[0]


def test_add_smoothness_zero_smoofac_still_runs():
    """smoofac=0 must not raise and must add at least boundary rows."""
    A_o = np.eye(6)
    A_c = np.zeros((6, 3))
    d = np.zeros(6)
    A_o2, A_c2, d2 = _add_smoothness(A_o, A_c, d, smoofac=0.0)
    assert A_o2.shape[0] >= 6


def test_add_zero_mean_appends_one_row():
    """Zero-mean adds exactly one constraint row."""
    A_o = np.eye(5)
    A_c = np.zeros((5, 3))
    d = np.ones(5)
    A_o2, A_c2, d2 = _add_zero_mean(A_o, A_c, d)
    assert A_o2.shape[0] == 6
    assert d2[-1] == 0.0  # RHS = 0 for zero-mean
```

- [ ] **Step 2: Verify failures**

```
pytest tests/test_inverse.py::test_add_smoothness_increases_rows -v
```
Expected: `ImportError`

- [ ] **Step 3: Implement**

Add to `src/ladcp/solution/inverse.py`:

```python
def _add_smoothness(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    smoofac: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Append curvature-penalty rows to the ocean-velocity block (lainsmoo).

    For each interior column j (1..n_cols-2), adds one row with stencil
    [-1, 2, -1] scaled by smoofac * (median_norm / col_norm[j]).
    Also smooths the CTD block symmetrically (MATLAB calls lainsmoo twice).
    Boundary columns get first-derivative (slope) rows: [2,-2] and [-2,2].
    """
    def _smoo_one(A_target: np.ndarray, A_other: np.ndarray,
                  d_in: np.ndarray, fs0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n_rows, n_cols = A_target.shape
        col_norms = np.sqrt(np.abs(A_target).sum(axis=0))
        pos = col_norms > 0
        if not pos.any():
            return A_target, A_other, d_in
        median_norm = max(float(np.median(col_norms[pos])), 0.01)
        clipped = np.maximum(col_norms, median_norm * 0.1)
        fs = (median_norm / clipped) * fs0

        # Interior: curvature stencil [-1, 2, -1]
        cur = np.array([-1.0, 2.0, -1.0])
        smoo_rows_t, smoo_rows_o = [], []
        for j in range(1, n_cols - 1):
            if fs[j] > 0:
                row = np.zeros(n_cols)
                row[j - 1 : j + 2] = cur * fs[j]
                smoo_rows_t.append(row)
                smoo_rows_o.append(np.zeros(A_other.shape[1]))

        # Boundaries: slope constraint [2,-2] / [-2,2]
        if fs[0] > 0:
            row = np.zeros(n_cols); row[0:2] = [2.0, -2.0] * fs[0]
            smoo_rows_t.append(row)
            smoo_rows_o.append(np.zeros(A_other.shape[1]))
        if fs[-1] > 0:
            row = np.zeros(n_cols); row[-2:] = [-2.0, 2.0] * fs[-1]
            smoo_rows_t.append(row)
            smoo_rows_o.append(np.zeros(A_other.shape[1]))

        if not smoo_rows_t:
            return A_target, A_other, d_in

        block_t = np.array(smoo_rows_t)
        block_o = np.array(smoo_rows_o)
        return (
            np.vstack([A_target, block_t]),
            np.vstack([A_other, block_o]),
            np.concatenate([d_in, np.zeros(len(smoo_rows_t))]),
        )

    # Smooth ocean velocity, then CTD velocity (MATLAB calls lainsmoo twice)
    A_ocean, A_ctd, d = _smoo_one(A_ocean, A_ctd, d, smoofac)
    A_ctd, A_ocean, d = _smoo_one(A_ctd, A_ocean, d, smoofac)
    return A_ocean, A_ctd, d


def _add_zero_mean(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    weight: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Constrain mean(u_ocean) = 0 when no external velocity reference (lainocean).

    Appends one row: sum(A_ocean columns) * weight / n_zbins = 0.
    """
    n_zbins = A_ocean.shape[1]
    scale = float(np.mean(np.abs(A_ocean).sum(axis=0)))
    row_o = np.ones(n_zbins) * weight * scale / n_zbins
    row_c = np.zeros(A_ctd.shape[1])
    return (
        np.vstack([A_ocean, row_o[np.newaxis, :]]),
        np.vstack([A_ctd, row_c[np.newaxis, :]]),
        np.concatenate([d, [0.0]]),
    )
```

- [ ] **Step 4: Run tests**

```
pytest tests/test_inverse.py -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```
git add src/ladcp/solution/inverse.py tests/test_inverse.py
git commit -m "feat: smoothness curvature and zero-mean fallback constraints"
```

---

### Task 4: Bottom-track + barotropic GPS constraints

**Bottom-track** (`lainbott`): each super-ensemble with valid bottom-track velocity gets one constraint row that ties `u_CTD(e) = u_bottom(e)`. The constraint weight is `botfac * velerr / bvels_e` (inversely proportional to bottom-track std).

**Barotropic** (`lainbaro`): one constraint row asserting `sum(u_CTD(e) * dt(e) / T) = u_ship`, where T = total cast duration. This ties the time-mean CTD velocity to the GPS-derived ship velocity.

**Files:**
- Modify: `src/ladcp/solution/inverse.py`
- Modify: `tests/test_inverse.py`

**Interfaces:**
- Produces (module-private):
  - `_add_bottom_track(A_ocean, A_ctd, d, bvel, bvels, *, botfac, velerr) -> (A_ocean, A_ctd, d)`
  - `_add_barotropic(A_ocean, A_ctd, d, u_ship, v_ship, dt, *, barofac) -> (A_ocean_u, A_ctd_u, d_u, A_ocean_v, A_ctd_v, d_v)`

Note: barotropic constraint has different RHS for U and V, so it returns two sets of arrays.

- [ ] **Step 1: Add tests**

Append to `tests/test_inverse.py`:

```python
from ladcp.solution.inverse import _add_bottom_track, _add_barotropic


def test_add_bottom_track_appends_rows():
    """One constraint row per ensemble with valid bottom track."""
    n_obs, n_zbins, n_se = 10, 5, 4
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.zeros((n_obs, n_se))
    d = np.zeros(n_obs)
    bvel = np.array([0.1, np.nan, -0.2, 0.0])   # 3 valid, 1 NaN
    bvels = np.array([0.01, 0.01, 0.01, 0.01])
    A_o2, A_c2, d2 = _add_bottom_track(A_o, A_c, d, bvel, bvels,
                                        botfac=1.0, velerr=0.05)
    # Should append one row per finite bvel
    n_finite = int(np.sum(np.isfinite(bvel)))
    assert A_c2.shape[0] == n_obs + n_finite


def test_add_barotropic_appends_one_row():
    """Barotropic constraint adds exactly one row to each component."""
    n_obs, n_zbins, n_se = 10, 5, 4
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.zeros((n_obs, n_se))
    d_u = np.zeros(n_obs)
    d_v = np.zeros(n_obs)
    dt = np.ones(n_se) * 100.0
    A_ou, A_cu, du, A_ov, A_cv, dv = _add_barotropic(
        A_o, A_c, d_u, A_o.copy(), A_c.copy(), d_v,
        u_ship=0.5, v_ship=0.1, dt=dt, barofac=1.0,
    )
    assert A_cu.shape[0] == n_obs + 1
    assert du[-1] != 0.0   # RHS = -u_ship * weight
```

- [ ] **Step 2: Verify failures**

```
pytest tests/test_inverse.py::test_add_bottom_track_appends_rows -v
```
Expected: `ImportError`

- [ ] **Step 3: Implement**

Add to `src/ladcp/solution/inverse.py`:

```python
def _add_bottom_track(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    bvel: np.ndarray,
    bvels: np.ndarray,
    *,
    botfac: float = 1.0,
    velerr: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Constrain CTD velocity to bottom-track velocity where available (lainbott).

    Each valid ensemble e gets one row: A_ctd[new_row, e] = weight_e
    with d[new_row] = bvel[e] * weight_e.
    Weight = botfac * velerr / bvels[e] scaled by sqrt(sum(|A_ctd columns|)).

    Parameters
    ----------
    bvel  : (n_se,) bottom-track velocity component (NaN = no measurement)
    bvels : (n_se,) bottom-track velocity std (m/s)
    """
    n_se = A_ctd.shape[1]
    valid = np.isfinite(bvel) & np.isfinite(bvels) & (bvels > 0)
    if not valid.any():
        return A_ocean, A_ctd, d

    col_scale = np.sqrt(np.abs(A_ctd).sum(axis=0))  # (n_se,)

    rows_o, rows_c, rhs = [], [], []
    for e in np.where(valid)[0]:
        weight_e = botfac * (velerr / bvels[e]) * col_scale[e]
        row_c = np.zeros(n_se)
        row_c[e] = weight_e
        rows_c.append(row_c)
        rows_o.append(np.zeros(A_ocean.shape[1]))
        rhs.append(bvel[e] * weight_e)

    return (
        np.vstack([A_ocean, rows_o]),
        np.vstack([A_ctd, rows_c]),
        np.concatenate([d, rhs]),
    )


def _add_barotropic(
    A_ocean_u: np.ndarray,
    A_ctd_u: np.ndarray,
    d_u: np.ndarray,
    A_ocean_v: np.ndarray,
    A_ctd_v: np.ndarray,
    d_v: np.ndarray,
    *,
    u_ship: float,
    v_ship: float,
    dt: np.ndarray,
    barofac: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Constrain time-mean CTD velocity = GPS-derived ship velocity (lainbaro).

    Appends one row: sum(A_ctd * dt / T) = u_ship.
    Weight scaled by sqrt(sum of CTD column norms) to balance the system.
    """
    T = float(dt.sum())
    n_se = A_ctd_u.shape[1]
    col_scale = np.sqrt(np.abs(A_ctd_u).sum(axis=0))  # (n_se,)
    fac = float(np.sqrt(np.sum(col_scale)))

    # Barotropic row: dt[e]/T per column, scaled
    row_c = dt / T * barofac * fac
    row_o = np.zeros(A_ocean_u.shape[1])

    A_ocean_u2 = np.vstack([A_ocean_u, row_o[np.newaxis, :]])
    A_ctd_u2 = np.vstack([A_ctd_u, row_c[np.newaxis, :]])
    d_u2 = np.concatenate([d_u, [-u_ship * barofac * fac]])

    A_ocean_v2 = np.vstack([A_ocean_v, row_o[np.newaxis, :]])
    A_ctd_v2 = np.vstack([A_ctd_v, row_c[np.newaxis, :]])
    d_v2 = np.concatenate([d_v, [-v_ship * barofac * fac]])

    return A_ocean_u2, A_ctd_u2, d_u2, A_ocean_v2, A_ctd_v2, d_v2
```

- [ ] **Step 4: Run tests**

```
pytest tests/test_inverse.py -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```
git add src/ladcp/solution/inverse.py tests/test_inverse.py
git commit -m "feat: bottom-track and barotropic GPS velocity constraints"
```

---

### Task 5: Least-squares solver (`_solve_lsq`)

Solves `d = A * m` using SciPy's `lstsq`. Also computes per-parameter error estimates via the normal equations `me = sqrt(diag(inv(A'A)) * ||d - Am||² / (n-p))`, matching MATLAB's `lesqfit()`.

**Files:**
- Modify: `src/ladcp/solution/inverse.py`
- Modify: `tests/test_inverse.py`

**Interfaces:**
- Produces: `_solve_lsq(A, d) -> tuple[np.ndarray, np.ndarray]` — (m, me)

- [ ] **Step 1: Add tests**

Append to `tests/test_inverse.py`:

```python
from ladcp.solution.inverse import _solve_lsq


def test_solve_lsq_identity():
    """Identity system must recover exact solution."""
    A = np.eye(4)
    d = np.array([1.0, 2.0, 3.0, 4.0])
    m, me = _solve_lsq(A, d)
    assert np.allclose(m, d)
    assert np.all(me >= 0)


def test_solve_lsq_overdetermined_consistent():
    """Consistent overdetermined system must find exact solution."""
    # 3 equations, 2 unknowns: x=1, y=2
    A = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    d = np.array([1.0, 2.0, 3.0])
    m, me = _solve_lsq(A, d)
    assert np.allclose(m, [1.0, 2.0], atol=1e-10)
    assert np.all(np.isfinite(me))


def test_solve_lsq_error_shape():
    """Error vector must have same length as solution vector."""
    A = np.random.rand(20, 5)
    d = np.random.rand(20)
    m, me = _solve_lsq(A, d)
    assert m.shape == me.shape == (5,)
```

- [ ] **Step 2: Verify failures**

```
pytest tests/test_inverse.py::test_solve_lsq_identity -v
```
Expected: `ImportError`

- [ ] **Step 3: Implement**

Add to `src/ladcp/solution/inverse.py`:

```python
import scipy.linalg


def _solve_lsq(
    A: np.ndarray,
    d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the least-squares system d = A*m (replicates lesqfit() + lainsolv()).

    Returns
    -------
    m  : (n_params,) solution vector
    me : (n_params,) 1-sigma parameter error estimates

    Error formula matches MATLAB lesqfit:
        me = sqrt(diag(inv(A'A)) * ||d - Am||² / (n - p))
    """
    m, _, _, _ = scipy.linalg.lstsq(A, d, check_finite=False)

    # Error estimate via normal equations
    dm = A @ m
    n, p = A.shape
    dof = max(n - p, 1)
    sigma2 = float(np.sum((d - dm) ** 2) / dof)
    try:
        AtA_inv = np.linalg.inv(A.T @ A)
        me = np.sqrt(np.abs(np.diag(AtA_inv)) * sigma2)
    except np.linalg.LinAlgError:
        me = np.full(p, np.nan)

    return m, me
```

- [ ] **Step 4: Run tests**

```
pytest tests/test_inverse.py -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```
git add src/ladcp/solution/inverse.py tests/test_inverse.py
git commit -m "feat: least-squares solver replicating MATLAB lesqfit()"
```

---

### Task 6: `InverseParams` + `InverseResult` + `compute_inverse()` entry point

Wires all prior tasks into a single callable. Also solves down-cast and up-cast separately with a zero-mean fallback constraint (baroclinic separation).

**Files:**
- Modify: `src/ladcp/solution/inverse.py`
- Modify: `src/ladcp/solution/__init__.py`
- Modify: `tests/test_inverse.py`

**Interfaces:**
- Produces:
  - `InverseParams` dataclass
  - `InverseResult` dataclass
  - `compute_inverse(se: SuperEnsemble, *, params: InverseParams | None = None, u_ship: float = 0.0, v_ship: float = 0.0) -> InverseResult`

- [ ] **Step 1: Add tests**

Append to `tests/test_inverse.py`:

```python
from ladcp.solution.inverse import InverseParams, InverseResult, compute_inverse


def test_compute_inverse_returns_result():
    """compute_inverse must return an InverseResult with z, u, v arrays."""
    ens = _make_ens(n_bins=5, n_ens=60)
    se = prepare_superensembles(ens, dz=10.0)
    result = compute_inverse(se)
    assert isinstance(result, InverseResult)
    assert result.z.shape == result.u.shape == result.v.shape
    assert len(result.z) > 0


def test_compute_inverse_zero_mean_no_constraint():
    """Without external constraints, result mean should be near zero."""
    ens = _make_ens(n_bins=5, n_ens=60, u_val=0.0, v_val=0.0)
    se = prepare_superensembles(ens, dz=10.0)
    params = InverseParams(barofac=0.0, botfac=0.0, sadcpfac=0.0)
    result = compute_inverse(se, params=params)
    assert abs(result.ubar) < 0.05
    assert abs(result.vbar) < 0.05


def test_compute_inverse_down_up_same_length():
    """Down-cast and up-cast profiles must have same length as full profile."""
    ens = _make_ens(n_bins=5, n_ens=60)
    se = prepare_superensembles(ens, dz=10.0)
    result = compute_inverse(se)
    assert result.u_do.shape == result.u.shape
    assert result.u_up.shape == result.u.shape
```

- [ ] **Step 2: Verify failures**

```
pytest tests/test_inverse.py::test_compute_inverse_returns_result -v
```
Expected: `ImportError`

- [ ] **Step 3: Implement `InverseParams`, `InverseResult`, and `compute_inverse`**

Add to `src/ladcp/solution/inverse.py`:

```python
@dataclass
class InverseParams:
    """Tuning parameters for compute_inverse() (getinv.m ps struct)."""
    dz: float = 10.0          # depth bin size m
    botfac: float = 1.0       # bottom-track constraint weight (0 = disable)
    sadcpfac: float = 1.0     # SADCP constraint weight (0 = disable)
    barofac: float = 1.0      # GPS barotropic constraint weight (0 = disable)
    smoofac: float = 0.0      # curvature smoothing weight (0 = minimal)
    velerr: float = 0.05      # nominal velocity error m/s
    weightmin: float = 0.05   # minimum observation weight threshold
    nav_error: float = 30.0   # navigation error m (for barvelerr computation)
    down_up: bool = True       # also solve down-cast and up-cast separately


@dataclass
class InverseResult:
    """Output of compute_inverse() (getinv.m dr struct)."""
    z: np.ndarray       # (n_zbins,) depth m (positive, increasing downward)
    u: np.ndarray       # (n_zbins,) eastward velocity m/s
    v: np.ndarray       # (n_zbins,) northward velocity m/s
    uerr: np.ndarray    # (n_zbins,) velocity error estimate m/s
    nvel: np.ndarray    # (n_zbins,) number of observations per depth bin
    u_do: np.ndarray    # (n_zbins,) downcast-only eastward velocity
    v_do: np.ndarray    # (n_zbins,) downcast-only northward velocity
    u_up: np.ndarray    # (n_zbins,) upcast-only eastward velocity
    v_up: np.ndarray    # (n_zbins,) upcast-only northward velocity
    u_ctd: np.ndarray   # (n_se,) CTD eastward velocity m/s
    v_ctd: np.ndarray   # (n_se,) CTD northward velocity m/s
    ubar: float         # depth-mean eastward velocity
    vbar: float         # depth-mean northward velocity
    zctd: np.ndarray    # (n_se,) CTD depth time series m
    wctd: np.ndarray    # (n_se,) CTD vertical velocity m/s


def compute_inverse(
    se: SuperEnsemble,
    *,
    params: InverseParams | None = None,
    u_ship: float = 0.0,
    v_ship: float = 0.0,
) -> InverseResult:
    """Solve the LADCP inverse velocity problem (replicates getinv.m).

    Parameters
    ----------
    se      : SuperEnsemble from prepare_superensembles().
    params  : Tuning parameters; defaults to InverseParams().
    u_ship  : Mean eastward ship velocity m/s over cast (from GPS start/end).
    v_ship  : Mean northward ship velocity m/s over cast.
    """
    if params is None:
        params = InverseParams()

    n_se = se.izm.shape[1]

    # --- Flatten observations and build A matrices ---
    d_u, d_v, izv, jprof, wm = _flatten_obs(se, params.velerr, params.weightmin)
    if len(d_u) < 10:
        raise ValueError("Too few valid observations for inversion")

    A_ocean_sp = _build_obs_matrix(izv, params.dz)
    A_ctd_sp = _build_ctd_matrix(jprof, n_se)
    n_zbins = A_ocean_sp.shape[1]

    A_o_u, A_c_u, dw_u, idx_down, idx_up = _apply_weights(
        A_ocean_sp, A_ctd_sp, d_u, wm
    )
    A_o_v, A_c_v, dw_v, _, _ = _apply_weights(A_ocean_sp, A_ctd_sp, d_v, wm)

    # --- Depth vector for output ---
    z = np.arange(1, n_zbins + 1) * params.dz  # positive, 1-indexed depth bins

    # --- Smoothness constraints (applied to both U and V identically) ---
    A_o_u, A_c_u, dw_u = _add_smoothness(A_o_u, A_c_u, dw_u, params.smoofac)
    A_o_v, A_c_v, dw_v = _add_smoothness(A_o_v, A_c_v, dw_v, params.smoofac)

    # --- Bottom-track constraint ---
    has_btrack = params.botfac > 0 and np.any(np.isfinite(se.bvel[0]))
    if has_btrack:
        A_o_u, A_c_u, dw_u = _add_bottom_track(
            A_o_u, A_c_u, dw_u,
            se.bvel[0], se.bvels[0], botfac=params.botfac, velerr=params.velerr,
        )
        A_o_v, A_c_v, dw_v = _add_bottom_track(
            A_o_v, A_c_v, dw_v,
            se.bvel[1], se.bvels[1], botfac=params.botfac, velerr=params.velerr,
        )

    # --- Barotropic (GPS) constraint ---
    has_baro = params.barofac > 0 and (u_ship != 0.0 or v_ship != 0.0)
    if has_baro:
        A_o_u, A_c_u, dw_u, A_o_v, A_c_v, dw_v = _add_barotropic(
            A_o_u, A_c_u, dw_u, A_o_v, A_c_v, dw_v,
            u_ship=u_ship, v_ship=v_ship, dt=se.dt, barofac=params.barofac,
        )

    # --- Zero-mean fallback when no external constraint ---
    if not has_btrack and not has_baro:
        A_o_u, A_c_u, dw_u = _add_zero_mean(A_o_u, A_c_u, dw_u)
        A_o_v, A_c_v, dw_v = _add_zero_mean(A_o_v, A_c_v, dw_v)

    # --- Solve full-cast system ---
    A_full_u = np.hstack([A_o_u, A_c_u])
    A_full_v = np.hstack([A_o_v, A_c_v])
    m_u, me_u = _solve_lsq(A_full_u, dw_u)
    m_v, me_v = _solve_lsq(A_full_v, dw_v)

    u_ocean = m_u[:n_zbins]
    v_ocean = m_v[:n_zbins]
    u_ctd_neg = m_u[n_zbins:]
    v_ctd_neg = m_v[n_zbins:]
    # Sign convention: solved u_ctd_neg = -u_CTD  (see MATLAB dr.uctd = -real(uctd))
    u_ctd = -u_ctd_neg[:n_se] if len(u_ctd_neg) >= n_se else np.full(n_se, np.nan)
    v_ctd = -v_ctd_neg[:n_se] if len(v_ctd_neg) >= n_se else np.full(n_se, np.nan)

    uerr = np.sqrt(me_u[:n_zbins] ** 2 + me_v[:n_zbins] ** 2)
    nvel = np.asarray(A_ocean_sp.sum(axis=0)).ravel()

    # --- Down/up cast separately (ps.down_up=1) ---
    _BAROCLINIC_FAC = 10.0  # MATLAB: baroclinfac = 10 (large = forces zero mean)

    def _solve_subset(idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        A_os = A_o_u[idx, :n_zbins]
        A_cs = A_c_u[idx, :]
        ds_u = dw_u[idx]
        A_os2 = A_o_v[idx, :n_zbins]
        A_cs2 = A_c_v[idx, :]
        ds_v = dw_v[idx]

        # Remove zero-constrained CTD columns
        active = np.where(np.abs(A_cs).sum(axis=0) > 0)[0]
        A_cs = A_cs[:, active]
        A_cs2 = A_cs2[:, active]

        # Add zero-mean baroclinic constraint
        zero_row = np.zeros(n_zbins)
        A_os = np.vstack([A_os, (_BAROCLINIC_FAC * np.ones(n_zbins))[np.newaxis, :]])
        A_cs = np.vstack([A_cs, np.zeros((1, A_cs.shape[1]))])
        ds_u = np.concatenate([ds_u, [0.0]])
        A_os2 = np.vstack([A_os2, (_BAROCLINIC_FAC * np.ones(n_zbins))[np.newaxis, :]])
        A_cs2 = np.vstack([A_cs2, np.zeros((1, A_cs2.shape[1]))])
        ds_v = np.concatenate([ds_v, [0.0]])

        if A_os.shape[0] < A_os.shape[1] + 2:
            return np.full(n_zbins, np.nan), np.full(n_zbins, np.nan)

        mu, _ = _solve_lsq(np.hstack([A_os, A_cs]), ds_u)
        mv, _ = _solve_lsq(np.hstack([A_os2, A_cs2]), ds_v)
        u_sub = mu[:n_zbins]
        v_sub = mv[:n_zbins]
        # Clip unreasonably large values
        u_sub[np.abs(u_sub) > 5.0] = np.nan
        v_sub[np.abs(v_sub) > 5.0] = np.nan
        return u_sub, v_sub

    if params.down_up and len(idx_down) > 5 and len(idx_up) > 5:
        u_do, v_do = _solve_subset(idx_down)
        u_up, v_up = _solve_subset(idx_up)
    else:
        u_do = u_ocean.copy()
        v_do = v_ocean.copy()
        u_up = u_ocean.copy()
        v_up = v_ocean.copy()

    return InverseResult(
        z=z,
        u=u_ocean, v=v_ocean, uerr=uerr, nvel=nvel,
        u_do=u_do, v_do=v_do, u_up=u_up, v_up=v_up,
        u_ctd=u_ctd, v_ctd=v_ctd,
        ubar=float(np.nanmean(u_ocean)),
        vbar=float(np.nanmean(v_ocean)),
        zctd=se.z,
        wctd=-np.nanmean(se.rw, axis=0),
    )
```

- [ ] **Step 4: Export from `src/ladcp/solution/__init__.py`**

Read the current `__init__.py` first, then add:

```python
from ladcp.solution.inverse import (
    EnsembleData,
    SuperEnsemble,
    InverseParams,
    InverseResult,
    prepare_superensembles,
    compute_inverse,
)

__all__ = [
    # existing shear exports...
    "EnsembleData",
    "SuperEnsemble",
    "InverseParams",
    "InverseResult",
    "prepare_superensembles",
    "compute_inverse",
]
```

- [ ] **Step 5: Run all tests**

```
pytest tests/test_inverse.py -v
```
Expected: all PASS (≥16 tests)

- [ ] **Step 6: Commit**

```
git add src/ladcp/solution/inverse.py src/ladcp/solution/__init__.py tests/test_inverse.py
git commit -m "feat: compute_inverse() entry point with InverseParams/InverseResult"
```

---

### Task 7: SADCP constraint (`_add_sadcp`)

Ship ADCP velocity is an independent in-situ ocean velocity profile that can directly constrain `u_ocean(z)`. Each SADCP depth bin where the instrument measured velocity gets one row in A_ocean (not A_ctd) with weight `sadcpfac * velerr / sadcp_err[j]`.

**Files:**
- Modify: `src/ladcp/solution/inverse.py`
- Modify: `tests/test_inverse.py`

**Interfaces:**
- Produces: `_add_sadcp(A_ocean, A_ctd, d, sadcp_z, sadcp_vel, sadcp_err, dz, *, sadcpfac, velerr) -> (A_ocean, A_ctd, d)`
- `compute_inverse` gains an optional `sadcp` parameter.

- [ ] **Step 1: Add tests**

Append to `tests/test_inverse.py`:

```python
from ladcp.solution.inverse import _add_sadcp


def test_add_sadcp_appends_rows():
    """Each finite SADCP measurement adds one row to A_ocean."""
    n_obs, n_zbins, n_se = 10, 8, 4
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.zeros((n_obs, n_se))
    d = np.zeros(n_obs)
    sadcp_z = np.array([15.0, 25.0, 35.0, np.nan])  # 3 valid
    sadcp_vel = np.array([0.1, 0.2, 0.3, np.nan])
    sadcp_err = np.array([0.02, 0.02, 0.02, 0.02])
    A_o2, A_c2, d2 = _add_sadcp(
        A_o, A_c, d,
        sadcp_z=sadcp_z, sadcp_vel=sadcp_vel, sadcp_err=sadcp_err,
        dz=10.0, sadcpfac=1.0, velerr=0.05,
    )
    assert A_o2.shape[0] == n_obs + 3  # 3 finite measurements


def test_add_sadcp_zeros_in_A_ctd():
    """SADCP constraint rows must have no A_ctd contribution."""
    n_obs, n_zbins, n_se = 5, 4, 3
    A_o = np.zeros((n_obs, n_zbins))
    A_c = np.zeros((n_obs, n_se))
    d = np.zeros(n_obs)
    A_o2, A_c2, d2 = _add_sadcp(
        A_o, A_c, d,
        sadcp_z=np.array([10.0]), sadcp_vel=np.array([0.5]),
        sadcp_err=np.array([0.02]), dz=10.0, sadcpfac=1.0, velerr=0.05,
    )
    assert A_c2[-1].sum() == 0.0  # last row of A_ctd = zeros
```

- [ ] **Step 2: Verify failures**

```
pytest tests/test_inverse.py::test_add_sadcp_appends_rows -v
```
Expected: `ImportError`

- [ ] **Step 3: Implement `_add_sadcp` and wire into `compute_inverse`**

Add to `src/ladcp/solution/inverse.py`:

```python
def _add_sadcp(
    A_ocean: np.ndarray,
    A_ctd: np.ndarray,
    d: np.ndarray,
    *,
    sadcp_z: np.ndarray,
    sadcp_vel: np.ndarray,
    sadcp_err: np.ndarray,
    dz: float,
    sadcpfac: float = 1.0,
    velerr: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Constrain u_ocean at SADCP depth bins (lainsadcp from getinv.m).

    Each valid SADCP measurement at depth z_j gets one row in A_ocean at
    column round(z_j / dz) - 1 with weight sadcpfac * velerr / sadcp_err[j].

    Parameters
    ----------
    sadcp_z   : (n_sadcp,) positive depth m
    sadcp_vel : (n_sadcp,) velocity component m/s
    sadcp_err : (n_sadcp,) velocity std m/s
    """
    n_zbins = A_ocean.shape[1]
    valid = np.isfinite(sadcp_z) & np.isfinite(sadcp_vel) & np.isfinite(sadcp_err)
    if not valid.any():
        return A_ocean, A_ctd, d

    col_scale = np.sqrt(np.abs(A_ocean).sum(axis=0))  # (n_zbins,)

    rows_o, rows_c, rhs = [], [], []
    for k in np.where(valid)[0]:
        j = int(np.round(sadcp_z[k] / dz)) - 1
        j = min(max(j, 0), n_zbins - 1)
        w = sadcpfac * (velerr / max(sadcp_err[k], 1e-6)) * col_scale[j]
        row_o = np.zeros(n_zbins)
        row_o[j] = w
        rows_o.append(row_o)
        rows_c.append(np.zeros(A_ctd.shape[1]))
        rhs.append(sadcp_vel[k] * w)

    return (
        np.vstack([A_ocean, rows_o]),
        np.vstack([A_ctd, rows_c]),
        np.concatenate([d, rhs]),
    )
```

Update `compute_inverse` signature and body to accept optional SADCP data:

```python
# In compute_inverse signature, add:
#   sadcp_z: np.ndarray | None = None,
#   sadcp_u: np.ndarray | None = None,
#   sadcp_v: np.ndarray | None = None,
#   sadcp_err: np.ndarray | None = None,

# After bottom-track block and before barotropic, add:
has_sadcp = (params.sadcpfac > 0 and sadcp_z is not None
             and np.any(np.isfinite(sadcp_z)))
if has_sadcp:
    A_o_u, A_c_u, dw_u = _add_sadcp(
        A_o_u, A_c_u, dw_u,
        sadcp_z=sadcp_z, sadcp_vel=sadcp_u, sadcp_err=sadcp_err,
        dz=params.dz, sadcpfac=params.sadcpfac, velerr=params.velerr,
    )
    A_o_v, A_c_v, dw_v = _add_sadcp(
        A_o_v, A_c_v, dw_v,
        sadcp_z=sadcp_z, sadcp_vel=sadcp_v, sadcp_err=sadcp_err,
        dz=params.dz, sadcpfac=params.sadcpfac, velerr=params.velerr,
    )

# Update zero-mean condition:
if not has_btrack and not has_baro and not has_sadcp:
    ...
```

- [ ] **Step 4: Run all tests**

```
pytest tests/test_inverse.py -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```
git add src/ladcp/solution/inverse.py tests/test_inverse.py
git commit -m "feat: SADCP velocity constraint; wire sadcp args into compute_inverse"
```

---

### Task 8: Integration test vs P16N cast 003 reference

Validates the full pipeline against the LDEO_IX-processed reference output `003.nc`. The test is gated on `TEST_DATA_DIR` (same pattern as the shear integration test).

**Files:**
- Create: `tests/integration/test_inverse_p16n_cast003.py`

**Interfaces:**
- Consumes: all public functions from Tasks 1 + 6.
- Reference data: `$TEST_DATA_DIR/003.nc` (LDEO_IX NetCDF output, P16N cast 003).
- Inputs: `$TEST_DATA_DIR/003DL000.000` (RDI PD0), `$TEST_DATA_DIR/003_01.cnv` (CTD).

**Tolerance:** RMS difference in u, v < 0.05 m/s over the depth range where the reference has ≥ 3 observations per bin.

- [ ] **Step 1: Write the integration test**

```python
# tests/integration/test_inverse_p16n_cast003.py
"""Integration test: inverse solver vs P16N cast 003 LDEO_IX reference."""
import os
import numpy as np
import pytest

DATA_DIR = os.environ.get("TEST_DATA_DIR", "")
SKIP = not bool(DATA_DIR)
SKIP_REASON = "TEST_DATA_DIR not set"

pytestmark = pytest.mark.skipif(SKIP, reason=SKIP_REASON)


@pytest.fixture(scope="module")
def reference():
    """Load LDEO_IX reference velocities from 003.nc."""
    import xarray as xr
    ds = xr.open_dataset(os.path.join(DATA_DIR, "003.nc"))
    return ds


@pytest.fixture(scope="module")
def inverse_result():
    """Run full pipeline on P16N cast 003 raw data."""
    import xarray as xr
    from ladcp.ingestion.rdi import load_rdi
    from ladcp.ingestion.ctd import load_ctd, assign_bin_depths
    from ladcp.transforms.beam2earth import beam2earth
    from ladcp.solution.inverse import EnsembleData, prepare_superensembles, compute_inverse

    dl_path = os.path.join(DATA_DIR, "003DL000.000")
    ctd_path = os.path.join(DATA_DIR, "003_01.cnv")

    rdi = load_rdi(dl_path)
    ctd = load_ctd(ctd_path)

    # Earth-frame velocities (u, v, w): shape (n_bins, n_ens)
    u, v, w = beam2earth(rdi)

    # Bin depths and quality weight from correlation
    z_m, izm = assign_bin_depths(rdi, ctd, looker="down")
    weight = np.nanmean(rdi.correlation, axis=2) / 128.0  # 0–1

    # CTD depth interpolated to each ADCP ensemble time
    from scipy.interpolate import interp1d
    f_z = interp1d(ctd.time_julian, ctd.pressure_dbar * (-1.0),
                   bounds_error=False, fill_value=np.nan)
    z_ctd = f_z(rdi.timestamp)

    # Bottom-track (u, v, w, std) — shape (n_ens, 3)/(n_ens, 3)
    bt_u = rdi.bottom_track[:, 0].astype(float) / 1000.0  # mm/s → m/s
    bt_v = rdi.bottom_track[:, 1].astype(float) / 1000.0
    bt_w = rdi.bottom_track[:, 2].astype(float) / 1000.0
    bvel = np.stack([bt_u, bt_v, bt_w], axis=1)
    bvels = np.full_like(bvel, 0.02)  # 2 cm/s nominal std
    hbot = rdi.bottom_range if hasattr(rdi, "bottom_range") else np.full(len(z_ctd), np.nan)

    ens = EnsembleData(
        u=u, v=v, w=w,
        weight=weight,
        izm=izm,
        z=z_ctd,
        time_jul=rdi.timestamp / 86400.0,   # s → Julian day (relative)
        bvel=bvel, bvels=bvels, hbot=hbot,
        izd=np.arange(u.shape[0]),
        izu=np.array([], dtype=int),
        slat=np.full(u.shape[1], np.nan),
        slon=np.full(u.shape[1], np.nan),
    )
    se = prepare_superensembles(ens, dz=16.0)
    return compute_inverse(se)


def test_inverse_u_rmse(inverse_result, reference):
    """RMS error in u vs LDEO_IX reference must be < 0.05 m/s."""
    ref_z = reference["depth"].values  # positive m
    ref_u = reference["u"].values      # m/s (adjust variable name to match 003.nc)

    result_u = np.interp(ref_z, inverse_result.z, inverse_result.u,
                         left=np.nan, right=np.nan)
    valid = np.isfinite(ref_u) & np.isfinite(result_u)
    assert valid.sum() > 10, "Too few overlapping depth bins to compare"

    rmse = float(np.sqrt(np.mean((result_u[valid] - ref_u[valid]) ** 2)))
    assert rmse < 0.05, f"u RMSE {rmse:.4f} m/s exceeds 0.05 m/s tolerance"


def test_inverse_v_rmse(inverse_result, reference):
    """RMS error in v vs LDEO_IX reference must be < 0.05 m/s."""
    ref_z = reference["depth"].values
    ref_v = reference["v"].values

    result_v = np.interp(ref_z, inverse_result.z, inverse_result.v,
                         left=np.nan, right=np.nan)
    valid = np.isfinite(ref_v) & np.isfinite(result_v)
    assert valid.sum() > 10

    rmse = float(np.sqrt(np.mean((result_v[valid] - ref_v[valid]) ** 2)))
    assert rmse < 0.05, f"v RMSE {rmse:.4f} m/s exceeds 0.05 m/s tolerance"


def test_inverse_profile_has_reasonable_depth_range(inverse_result, reference):
    """Computed profile depth range should cover most of the reference range."""
    ref_z = reference["depth"].values
    ref_valid = np.isfinite(reference["u"].values)
    ref_max_depth = float(ref_z[ref_valid].max())

    result_max_depth = float(inverse_result.z.max())
    assert result_max_depth > 0.8 * ref_max_depth, (
        f"Profile only reaches {result_max_depth:.0f} m vs {ref_max_depth:.0f} m reference"
    )
```

- [ ] **Step 2: Run unit tests to confirm no regressions**

```
pytest tests/test_inverse.py -v
```
Expected: all PASS

- [ ] **Step 3: Run integration tests (requires TEST_DATA_DIR)**

```
set TEST_DATA_DIR=test_data/2015_P16N && pytest tests/integration/test_inverse_p16n_cast003.py -v
```

Expected: either PASS or SKIP (if `TEST_DATA_DIR` not set). The variable names in `003.nc` (`depth`, `u`, `v`) may need adjustment — run `python -c "import xarray as xr; print(xr.open_dataset('test_data/2015_P16N/003.nc'))"` to confirm actual variable names and update the test accordingly.

- [ ] **Step 4: Commit**

```
git add tests/integration/test_inverse_p16n_cast003.py
git commit -m "test: integration test for inverse solver vs P16N cast 003 reference"
```

---

## Self-Review

**Spec coverage:**
- ✅ `prepinv.m` depth-window averaging → Task 1
- ✅ `lainseta` / `lainweig` matrix construction → Task 2
- ✅ `lainsmoo` / `lainocean` → Task 3
- ✅ `lainbott` / `lainbaro` → Task 4
- ✅ `lesqfit` / `lainsolv` → Task 5
- ✅ Full `getinv.m` orchestration + down/up cast separation → Task 6
- ✅ `lainsadcp` → Task 7
- ✅ Integration test vs reference → Task 8
- ⚠️ Drag model (`laindrag`, `dragfac > 0`) — not included; `dragfac=0` is the default and the common case. Add as a future task if needed.
- ⚠️ Up/down looker heading rotation (`rotup2down`, `offsetup2down` in prepinv.m) — deferred; single-instrument casts or pre-aligned dual-instrument casts don't need it.

**Placeholder scan:** All steps contain code. No "TBD" or "implement later" phrases.

**Type consistency check:**
- `EnsembleData.bvel` shape `(n_ens, 3)` → `SuperEnsemble.bvel` shape `(3, n_se)` — transpose happens in Task 1, used as `se.bvel[0]` (U) and `se.bvel[1]` (V) in Task 6. ✓
- `_flatten_obs` returns `jprof` as 0-indexed integers; `_build_ctd_matrix` expects the same. ✓
- `_apply_weights` returns dense `np.ndarray`; all downstream functions (`_add_smoothness`, `_add_bottom_track`, etc.) accept `np.ndarray`. ✓
- `compute_inverse` passes `se.bvel[0]` (shape `(n_se,)`) to `_add_bottom_track(bvel=...)` which iterates `np.where(np.isfinite(bvel))[0]`. ✓
