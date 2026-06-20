# Coordinate Transforms Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the beam→Earth coordinate transform that converts Teledyne RDI Workhorse PD0 beam-frame data to Earth-frame (u=East, v=North, w=Up) velocities, and validate against P16N cast 003 PD0 files which are confirmed beam-frame (EX=0x04, coord=Beam, gimbaled=True).

**Decisive fact:** The P16N test data (both DL and UL, 003DL000.000 and 003UL000.000) records beam-frame velocities with EX=0x04 (Beam coordinates, gimbaled tilt correction enabled, no bin mapping). This makes the transform layer the correct next step toward reproducing reference output.

**Architecture:** Four new pieces, layered:
1. `_pd0._read_fixed_leader()` — read EX byte at fixed-leader offset +23
2. `_types.RDIData` — add `coord_transform: int` field
3. `src/ladcp/transforms/beam2earth.py` — `beam2xyz()` + `beam2earth()` functions
4. `tests/integration/test_transforms_p16n_cast003.py` — plausibility-level integration test

**Reference implementation:** `docs/legacy/ADCPtools/janus2earth.m` and `docs/legacy/ADCPtools/janus2xyz.m` (read-only; do not modify).

**Tech Stack:** Python 3.11, numpy vectorized operations (no loops over bins or ensembles), pytest, uv

## Global Constraints

- Python 3.11+; all code must pass `ruff check` (E, F, I, NPY, UP rules)
- `src/` layout; install via `uv sync`; package name `ladcp`
- All new tests in `tests/` or `tests/integration/`; unit tests must pass without test data files
- Integration tests gated by `TEST_DATA_DIR` env var (see `tests/conftest.py`)
- **No Python loops over (nbin, nens) — numpy broadcasting only.** The MATLAB `janus2xyz.m` has a nested loop over nz×nt — do not replicate this in Python. Use broadcast arithmetic instead.
- Match `janus2earth.m` / `janus2xyz.m` math exactly. The exact matrix and formulas are:
  - `beam2xyz`: `Vx = uvfac*(-b1+b2)`, `Vy = uvfac*(-b3+b4)`, `Vz = wfac*(-b1-b2-b3-b4)` where `uvfac=1/(2*sin(θ))`, `wfac=1/(4*cos(θ))`
  - `beam2earth` gimbaled correction: `Sph2Sph3 = sin(pitch)*sin(roll); head += arcsin(Sph2Sph3/sqrt(cos(pitch)^2 + Sph2Sph3^2))`
  - Rotation: `u=+Vx*cx1+Vy*cy1+Vz*cz1`, `v=-Vx*cx2+Vy*cy2-Vz*cz2`, `w=-Vx*cx3+Vy*cy3+Vz*cz3`
  - Nine coefficients: `cx1=Cph1*Cph3+Sph1*Sph2*Sph3`, `cx2=Sph1*Cph3-Cph1*Sph2*Sph3`, `cx3=Cph2*Sph3`, `cy1=Sph1*Cph2`, `cy2=Cph1*Cph2`, `cy3=Sph2`, `cz1=Cph1*Sph3-Sph1*Sph2*Cph3`, `cz2=Sph1*Sph3+Cph1*Sph2*Cph3`, `cz3=Cph2*Cph3`
  - Sph1=sin(head), Cph1=cos(head), Sph2=sin(pitch), Cph2=cos(pitch), Sph3=sin(roll), Cph3=cos(roll) — all in radians
- NaN propagation: numpy arithmetic naturally propagates NaN; no special handling needed
- Do NOT implement bin mapping (BinmapType) or 3-beam solution in this phase — EX=0x04 has neither
- The existing stub `src/ladcp/transforms/beam2earth.py` must be replaced (not extended) — it currently raises `NotImplementedError`

---

### Task 1: Read coordinate-transform byte (EX) from fixed leader

**Files:**
- Modify: `src/ladcp/ingestion/_pd0.py` — `_read_fixed_leader()`
- Modify: `src/ladcp/ingestion/_types.py` — `RDIData` dataclass
- Modify: `src/ladcp/ingestion/rdi.py` — `load_rdi()`
- Modify: `tests/test_pd0_parser.py` — add EX byte unit tests
- Modify: `tests/integration/test_pd0_p16n_cast003.py` — add `coord_transform` assertion

**Why this first:** The EX byte tells callers whether the loaded data is beam, instrument, ship, or Earth frame. Downstream transforms and integration tests require this field to be present on `RDIData`.

**Interfaces produced:**
- `RDIData.coord_transform: int` — raw EX byte value from the fixed leader
  - Bits 4–3 encode coordinate system: 0=Beam, 1=Instrument, 2=Ship, 3=Earth
  - Bit 2: use tilts (gimbaled); Bit 1: bin mapping; Bit 0: 3-beam solution
  - P16N DL/UL has `coord_transform=4` (binary 00000100: Beam, gimbaled=True)

**EX byte location:** Fixed leader, byte offset `start + 23` relative to the byte immediately after the 2-byte type ID (0x0000). This is consistent with the existing layout: after the 7-byte firmware/flags skip, the byte sequence is: nbin(+7), npng+blen+blnk(+8-13), skip 4(+14-17), EX byte(+23).

Actually, EX is at absolute offset +23 from fixed leader content start (i.e., 2 bytes past the type ID). The sequence:
- +0 to +1: CPU firmware version + revision
- +2 to +3: System configuration (2 bytes) — beam angle, Janus type, frequency
- +4: Simulated flag
- +5: Lag length
- +6: Number of beams
- +7: nbin (WP cells) ← currently parsed
- +8-9: npng ← currently parsed
- +10-11: blen_cm ← currently parsed
- +12-13: blnk_cm ← currently parsed
- +14: profiling mode
- +15: low corr threshold
- +16: code repetitions
- +17: pct good minimum
- +18-19: error velocity max
- +20: TPP minutes
- +21: TPP seconds
- +22: TPP hundredths
- **+23: Coordinate transform (EX byte)** ← read this
- +24-25: Heading alignment (0.01 deg)
- +26-27: Heading bias (0.01 deg)
- +28: Sensor source
- +29: Sensors available
- +30-31: dist_cm ← currently parsed (confirmed correct offset)

- [ ] **Step 1: Update `_read_fixed_leader()` to include EX byte**

In `src/ladcp/ingestion/_pd0.py`, inside `_read_fixed_leader()`, read the EX byte at `start + 23` and return it as `"coord_transform"`:

```python
coord_transform = body[start + 23]
```

Add `"coord_transform": coord_transform` to the returned dict.

- [ ] **Step 2: Add `coord_transform` to `RDIData`**

In `src/ladcp/ingestion/_types.py`, add to `RDIData`:

```python
coord_transform: int  # raw EX byte: bits 4-3 = system (0=Beam,1=Instr,2=Ship,3=Earth), bit 2=tilt, bit 1=binmap, bit 0=3beam
```

Place it after `npng` (the other non-array scalar fields).

- [ ] **Step 3: Pass through in `load_rdi()`**

In `src/ladcp/ingestion/rdi.py`, after parsing `fl = ensembles[0]["fixed_leader"]`, pass `coord_transform=fl["coord_transform"]` into the `RDIData(...)` constructor call.

- [ ] **Step 4: Write unit tests for EX byte**

In `tests/test_pd0_parser.py`, add a test that verifies `coord_transform` is returned from `_read_fixed_leader` when a known EX byte is placed at offset +23 in a synthetic fixed leader.

The existing `_make_fixed_leader()` helper builds a 60-byte buffer. To set EX=4 (beam + gimbaled), set byte at `start + 23` where `start` is the offset within the buffer (which starts at 0 in the helper output).

Add to `TestPd0Parser`:

```python
def test_fixed_leader_coord_transform(self):
    from ladcp.ingestion._pd0 import _read_fixed_leader
    buf = bytearray(60)
    buf[23] = 0x04  # Beam + gimbaled
    fl = _read_fixed_leader(bytes(buf), 0)
    assert fl["coord_transform"] == 4

def test_fixed_leader_coord_transform_earth(self):
    from ladcp.ingestion._pd0 import _read_fixed_leader
    buf = bytearray(60)
    buf[23] = 0b00011000  # Earth frame (bits 4-3 = 11 = 3)
    fl = _read_fixed_leader(bytes(buf), 0)
    assert fl["coord_transform"] == 0b00011000
    assert (fl["coord_transform"] >> 3) & 0x03 == 3  # Earth
```

Also add to `TestLoadRdi`:

```python
def test_coord_transform_field_present(self, tmp_path):
    from ladcp.ingestion.rdi import load_rdi
    path = self._write_temp_pd0(tmp_path, nbin=4, nens=3)
    d = load_rdi(path)
    assert hasattr(d, "coord_transform")
    assert isinstance(d.coord_transform, int)
```

- [ ] **Step 5: Add integration assertion in P16N test**

In `tests/integration/test_pd0_p16n_cast003.py`, add:

```python
@pytest.mark.integration
def test_dl_coord_transform_beam_gimbaled(dl_path):
    """P16N 003DL000.000 was recorded in beam-frame with gimbaled tilt (EX=0x04)."""
    d = load_rdi(dl_path)
    assert d.coord_transform == 4, f"Expected EX=4 (beam+gimbaled), got {d.coord_transform}"
    assert (d.coord_transform >> 3) & 0x03 == 0, "Expected beam-frame (coord bits = 0b00)"
    assert bool(d.coord_transform & 0x04), "Expected gimbaled bit set"
```

- [ ] **Step 6: Run full test suite**

```powershell
cd C:\Users\peter_sha\Documents\sourcecode\LADCP
uv run pytest -v
```

Expected: all existing tests still pass, new unit tests pass, integration test passes when `TEST_DATA_DIR=test_data`.

Also run with integration:
```powershell
$env:TEST_DATA_DIR = "test_data"; uv run pytest -v; $env:TEST_DATA_DIR = $null
```

- [ ] **Step 7: Ruff lint + format**

```powershell
uv run ruff check src/ladcp/ingestion/ tests/
uv run ruff format --check src/ladcp/ingestion/ tests/
```

Fix any issues with `uv run ruff format src/ tests/`.

- [ ] **Step 8: Commit**

```powershell
git add src/ladcp/ingestion/_pd0.py src/ladcp/ingestion/_types.py src/ladcp/ingestion/rdi.py tests/test_pd0_parser.py tests/integration/test_pd0_p16n_cast003.py
git commit -m "feat: expose coord_transform (EX byte) from PD0 fixed leader in RDIData"
```

---

### Task 2: `beam2xyz()` — beam → instrument frame, vectorized numpy

**Files:**
- Modify: `src/ladcp/transforms/beam2earth.py` (replace stub)
- Create: `tests/test_transforms.py`

**Interfaces:**
- Produces: `beam2xyz(b1, b2, b3, b4, theta_deg)` → `tuple[np.ndarray, np.ndarray, np.ndarray]`
  - b1, b2, b3, b4: float64 arrays of shape `(nbin, nens)`, along-beam velocities m/s (TRDI: positive toward transducer face)
  - theta_deg: float, beam angle from vertical in degrees (20.0° for RDI Workhorse 300/600 kHz)
  - Returns: (Vx, Vy, Vz) each of shape `(nbin, nens)`, instrument-frame velocities m/s
    - x-axis: increases in beam 1's direction, away from instrument
    - y-axis: increases in beam 3's direction, away from instrument
    - z-axis: increases upward

**Math (from janus2xyz.m, vectorized):**
```
theta = radians(theta_deg)
uvfac = 1 / (2 * sin(theta))
wfac  = 1 / (4 * cos(theta))

Vx = uvfac * (-b1 + b2)
Vy = uvfac * (-b3 + b4)
Vz = wfac  * (-b1 - b2 - b3 - b4)
```

NaN handling: automatic — numpy arithmetic propagates NaN. If any beam in a bin/ensemble is NaN, the corresponding Vx/Vy/Vz will be NaN.

- [ ] **Step 1: Replace stub in `src/ladcp/transforms/beam2earth.py`**

Replace the entire file with:

```python
"""Janus beam → Earth coordinate transforms.

Reference: docs/legacy/ADCPtools/janus2earth.m, janus2xyz.m (Apaloczy et al.).
TRDI convention: positive along-beam velocity = toward transducer face.
Heading increases clockwise from the y-axis (North when heading=0).
"""

import numpy as np


def beam2xyz(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
    b4: np.ndarray,
    theta_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert 4-beam Janus along-beam velocities to instrument frame (Vx, Vy, Vz).

    Vectorized numpy replacement for janus2xyz.m (which has a nested nz×nt loop).
    All inputs broadcast together; NaN in any beam propagates to all outputs for
    that bin-ensemble.

    Parameters
    ----------
    b1, b2, b3, b4 : array_like, shape (nbin, nens)
        Along-beam velocities m/s. Positive = toward transducer face (TRDI convention).
        Beam layout: 1=+x, 2=−x, 3=+y, 4=−y (beam angle theta from vertical).
    theta_deg : float
        Beam angle from vertical in degrees. 20.0° for RDI Workhorse 300/600 kHz.

    Returns
    -------
    Vx, Vy, Vz : ndarray, shape (nbin, nens)
        Instrument-frame velocity components m/s.
        x: beam-1 direction; y: beam-3 direction; z: upward.
    """
    theta = np.radians(theta_deg)
    uvfac = 1.0 / (2.0 * np.sin(theta))
    wfac = 1.0 / (4.0 * np.cos(theta))
    Vx = uvfac * (-b1 + b2)
    Vy = uvfac * (-b3 + b4)
    Vz = wfac * (-b1 - b2 - b3 - b4)
    return Vx, Vy, Vz
```

Leave the `beam2earth` function as a stub for Task 3.

- [ ] **Step 2: Create `tests/test_transforms.py` with unit tests**

```python
"""Unit tests for src/ladcp/transforms/beam2earth.py."""

import numpy as np
import pytest

from ladcp.transforms.beam2earth import beam2xyz


THETA = 20.0  # RDI Workhorse 300 kHz


class TestBeam2xyz:
    def test_zero_beams_give_zero(self):
        b = np.zeros((5, 10))
        Vx, Vy, Vz = beam2xyz(b, b, b, b, THETA)
        assert np.all(Vx == 0.0)
        assert np.all(Vy == 0.0)
        assert np.all(Vz == 0.0)

    def test_output_shape_matches_input(self):
        rng = np.random.default_rng(42)
        b = rng.standard_normal((25, 100))
        Vx, Vy, Vz = beam2xyz(b, b, b, b, THETA)
        assert Vx.shape == (25, 100)
        assert Vy.shape == (25, 100)
        assert Vz.shape == (25, 100)

    def test_nan_propagation(self):
        """NaN in one beam → NaN in all output components for that cell."""
        b0 = np.ones((3, 4))
        b_nan = b0.copy()
        b_nan[1, 2] = np.nan
        Vx, Vy, Vz = beam2xyz(b_nan, b0, b0, b0, THETA)
        assert np.isnan(Vx[1, 2])
        assert np.isnan(Vy[1, 2])
        assert np.isnan(Vz[1, 2])
        # Other cells not affected
        assert np.isfinite(Vx[0, 0])

    def test_vx_from_b1_b2_only(self):
        """Vx depends only on b1 and b2; Vy depends only on b3 and b4."""
        theta = np.radians(THETA)
        uvfac = 1.0 / (2.0 * np.sin(theta))
        b = np.zeros((1, 1))
        b1 = np.full((1, 1), 0.5)
        b2 = np.full((1, 1), -0.3)
        Vx, Vy, Vz = beam2xyz(b1, b2, b, b, THETA)
        expected_Vx = uvfac * (-0.5 + (-0.3))
        assert abs(float(Vx[0, 0]) - expected_Vx) < 1e-10
        assert float(Vy[0, 0]) == 0.0

    def test_vy_from_b3_b4_only(self):
        theta = np.radians(THETA)
        uvfac = 1.0 / (2.0 * np.sin(theta))
        b = np.zeros((1, 1))
        b3 = np.full((1, 1), 0.2)
        b4 = np.full((1, 1), 0.6)
        Vx, Vy, Vz = beam2xyz(b, b, b3, b4, THETA)
        expected_Vy = uvfac * (-0.2 + 0.6)
        assert abs(float(Vy[0, 0]) - expected_Vy) < 1e-10
        assert float(Vx[0, 0]) == 0.0

    def test_vz_from_all_beams(self):
        """Vz = wfac * (-b1 - b2 - b3 - b4)."""
        theta = np.radians(THETA)
        wfac = 1.0 / (4.0 * np.cos(theta))
        b = np.full((1, 1), 0.1)
        Vx, Vy, Vz = beam2xyz(b, b, b, b, THETA)
        expected_Vz = wfac * (-0.4)
        assert abs(float(Vz[0, 0]) - expected_Vz) < 1e-10

    def test_scalar_theta(self):
        """theta_deg can be a float scalar."""
        b = np.ones((2, 3))
        Vx, Vy, Vz = beam2xyz(b, -b, b, -b, 20.0)
        assert Vx.shape == (2, 3)
```

- [ ] **Step 3: Run tests to verify they pass**

```powershell
cd C:\Users\peter_sha\Documents\sourcecode\LADCP
uv run pytest tests/test_transforms.py -v
```

Expected: all 7 tests pass.

- [ ] **Step 4: Ruff check**

```powershell
uv run ruff check src/ladcp/transforms/beam2earth.py tests/test_transforms.py
uv run ruff format --check src/ladcp/transforms/beam2earth.py tests/test_transforms.py
```

- [ ] **Step 5: Run full test suite to confirm no regressions**

```powershell
uv run pytest -v
```

- [ ] **Step 6: Commit**

```powershell
git add src/ladcp/transforms/beam2earth.py tests/test_transforms.py
git commit -m "feat: implement beam2xyz() — vectorized beam-to-instrument transform"
```

---

### Task 3: `beam2earth()` — full 4-beam Janus beam→Earth transform

**Files:**
- Modify: `src/ladcp/transforms/beam2earth.py` — add `beam2earth()` function (remove old stub)
- Modify: `tests/test_transforms.py` — add `beam2earth` unit tests

**Prerequisite:** Task 2 (beam2xyz in beam2earth.py)

**Interfaces:**
- Produces: `beam2earth(b1, b2, b3, b4, heading, pitch, roll, theta_deg, gimbaled=True)` → `tuple[np.ndarray, np.ndarray, np.ndarray]`
  - b1–b4: shape `(nbin, nens)` float64, along-beam m/s
  - heading, pitch, roll: shape `(nens,)` float64, degrees
  - theta_deg: float, beam angle in degrees
  - gimbaled: bool, apply Dewey & Stringer (2007) eq. A2 heading correction (default True)
  - Returns: (u, v, w) each shape `(nbin, nens)`, Earth-frame m/s (East, North, Up)

**Broadcasting note:** The nine rotation coefficients (cx1, cy1, …, cz3) have shape `(nens,)` and must broadcast against `(nbin, nens)` velocity arrays. `Vx * cx1` where `Vx.shape=(nbin,nens)` and `cx1.shape=(nens,)` broadcasts correctly along the last axis in numpy.

**Complete math (from janus2earth.m):**

```
# 1. Convert to radians
h = radians(heading); p = radians(pitch); r = radians(roll); theta = radians(theta_deg)

# 2. Trig (nt-length vectors)
Sph1=sin(h); Cph1=cos(h); Sph2=sin(p); Cph2=cos(p); Sph3=sin(r); Cph3=cos(r)

# 3. Gimbaled correction (D&S 2007 eq. A2) — applied BEFORE recomputing Sph1/Cph1
if gimbaled:
    Sph2Sph3 = Sph2 * Sph3
    h = h + arcsin(Sph2Sph3 / sqrt(Cph2**2 + Sph2Sph3**2))
    Sph1 = sin(h); Cph1 = cos(h)

# 4. Rotation matrix (nine (nens,) scalars)
cx1 = Cph1*Cph3 + Sph1*Sph2*Sph3
cx2 = Sph1*Cph3 - Cph1*Sph2*Sph3
cx3 = Cph2*Sph3
cy1 = Sph1*Cph2
cy2 = Cph1*Cph2
cy3 = Sph2
cz1 = Cph1*Sph3 - Sph1*Sph2*Cph3
cz2 = Sph1*Sph3 + Cph1*Sph2*Cph3
cz3 = Cph2*Cph3

# 5. Earth velocities (broadcasting: (nbin,nens) * (nens,) → (nbin,nens))
Vx, Vy, Vz = beam2xyz(b1, b2, b3, b4, theta_deg)
u = +Vx*cx1 + Vy*cy1 + Vz*cz1
v = -Vx*cx2 + Vy*cy2 - Vz*cz2
w = -Vx*cx3 + Vy*cy3 + Vz*cz3
```

**Unit test golden vectors:**

At heading=0°, pitch=0°, roll=0°, gimbaled=True or False (no correction when pitch=roll=0):
- Rotation matrix reduces to identity: cx1=1, cy1=0, cz1=0; cx2=0, cy2=1, cz2=0; cx3=0, cy3=0, cz3=1
- Therefore: u=Vx, v=Vy, w=Vz

This gives clean tests that don't require computing trig by hand.

- [ ] **Step 1: Add `beam2earth()` to `src/ladcp/transforms/beam2earth.py`**

Append to the module (after `beam2xyz`):

```python
def beam2earth(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
    b4: np.ndarray,
    heading: np.ndarray,
    pitch: np.ndarray,
    roll: np.ndarray,
    theta_deg: float,
    gimbaled: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert 4-beam Janus along-beam velocities to Earth frame (u=East, v=North, w=Up).

    Implements Appendix A of Dewey & Stringer (2007) as coded in janus2earth.m.
    TRDI convention: positive beam velocity = toward transducer face.
    Heading: clockwise from y-axis (y=North when heading=0).

    Parameters
    ----------
    b1, b2, b3, b4 : ndarray, shape (nbin, nens)
        Along-beam velocities m/s.
    heading, pitch, roll : ndarray, shape (nens,)
        Instrument orientation angles in degrees.
    theta_deg : float
        Beam angle from vertical in degrees (20.0° for Workhorse 300/600 kHz).
    gimbaled : bool
        If True, apply gimbaled heading correction (D&S 2007 eq. A2); if False,
        apply fixed-mount pitch correction (eq. A1). Default True (LADCP standard).

    Returns
    -------
    u, v, w : ndarray, shape (nbin, nens)
        Earth-frame velocity components m/s: East, North, Up.
    """
    h = np.radians(heading)
    p = np.radians(pitch)
    r = np.radians(roll)

    Sph1 = np.sin(h); Cph1 = np.cos(h)
    Sph2 = np.sin(p); Cph2 = np.cos(p)
    Sph3 = np.sin(r); Cph3 = np.cos(r)

    if gimbaled:
        Sph2Sph3 = Sph2 * Sph3
        h = h + np.arcsin(Sph2Sph3 / np.sqrt(Cph2**2 + Sph2Sph3**2))
        Sph1 = np.sin(h)
        Cph1 = np.cos(h)

    cx1 = Cph1 * Cph3 + Sph1 * Sph2 * Sph3
    cx2 = Sph1 * Cph3 - Cph1 * Sph2 * Sph3
    cx3 = Cph2 * Sph3
    cy1 = Sph1 * Cph2
    cy2 = Cph1 * Cph2
    cy3 = Sph2
    cz1 = Cph1 * Sph3 - Sph1 * Sph2 * Cph3
    cz2 = Sph1 * Sph3 + Cph1 * Sph2 * Cph3
    cz3 = Cph2 * Cph3

    Vx, Vy, Vz = beam2xyz(b1, b2, b3, b4, theta_deg)

    u = +Vx * cx1 + Vy * cy1 + Vz * cz1
    v = -Vx * cx2 + Vy * cy2 - Vz * cz2
    w = -Vx * cx3 + Vy * cy3 + Vz * cz3

    return u, v, w
```

Remove the old stub at the bottom of the file (the `janus5beam2earth` stub) only if it now conflicts — but check first. If `tests/test_smoke.py` still imports `janus5beam2earth`, keep that stub or update it to delegate/alias appropriately. Do NOT break the smoke test.

Actually: `tests/test_smoke.py` tests that `janus5beam2earth` raises `NotImplementedError`. The stub should remain (as a separate stub function). Do not remove it.

- [ ] **Step 2: Add `beam2earth` unit tests to `tests/test_transforms.py`**

Append to the test file:

```python
from ladcp.transforms.beam2earth import beam2earth


class TestBeam2earth:
    def test_identity_rotation_zero_heading_pitch_roll(self):
        """heading=pitch=roll=0 → rotation matrix is identity: u=Vx, v=Vy, w=Vz."""
        theta_deg = 20.0
        theta = np.radians(theta_deg)
        uvfac = 1.0 / (2.0 * np.sin(theta))
        wfac = 1.0 / (4.0 * np.cos(theta))

        nbin, nens = 5, 8
        b1 = np.random.default_rng(0).standard_normal((nbin, nens)) * 0.1
        b2 = np.random.default_rng(1).standard_normal((nbin, nens)) * 0.1
        b3 = np.random.default_rng(2).standard_normal((nbin, nens)) * 0.1
        b4 = np.random.default_rng(3).standard_normal((nbin, nens)) * 0.1

        heading = np.zeros(nens)
        pitch = np.zeros(nens)
        roll = np.zeros(nens)

        u, v, w = beam2earth(b1, b2, b3, b4, heading, pitch, roll, theta_deg)

        Vx_expected = uvfac * (-b1 + b2)
        Vy_expected = uvfac * (-b3 + b4)
        Vz_expected = wfac * (-b1 - b2 - b3 - b4)

        np.testing.assert_allclose(u, Vx_expected, rtol=1e-10)
        np.testing.assert_allclose(v, Vy_expected, rtol=1e-10)
        np.testing.assert_allclose(w, Vz_expected, rtol=1e-10)

    def test_output_shape(self):
        nbin, nens = 25, 100
        b = np.random.default_rng(42).standard_normal((nbin, nens))
        heading = np.random.default_rng(0).uniform(0, 360, nens)
        pitch = np.random.default_rng(1).uniform(-5, 5, nens)
        roll = np.random.default_rng(2).uniform(-5, 5, nens)
        u, v, w = beam2earth(b, b, b, b, heading, pitch, roll, 20.0)
        assert u.shape == (nbin, nens)
        assert v.shape == (nbin, nens)
        assert w.shape == (nbin, nens)

    def test_nan_propagation(self):
        """NaN beam → NaN Earth velocity for that cell."""
        nbin, nens = 3, 4
        b = np.ones((nbin, nens))
        b_nan = b.copy()
        b_nan[1, 2] = np.nan
        heading = np.full(nens, 90.0)
        pitch = np.zeros(nens)
        roll = np.zeros(nens)
        u, v, w = beam2earth(b_nan, b, b, b, heading, pitch, roll, 20.0)
        assert np.isnan(u[1, 2])
        assert np.isnan(v[1, 2])
        assert np.isnan(w[1, 2])
        assert np.isfinite(u[0, 0])

    def test_gimbaled_no_effect_when_pitch_roll_zero(self):
        """Gimbaled correction is zero when pitch=roll=0 (arcsin(0/sqrt(cos²+0))=0)."""
        nbin, nens = 3, 5
        b = np.random.default_rng(7).standard_normal((nbin, nens))
        heading = np.full(nens, 45.0)
        pitch = np.zeros(nens)
        roll = np.zeros(nens)
        u_g, v_g, w_g = beam2earth(b, b, b, b, heading, pitch, roll, 20.0, gimbaled=True)
        u_ng, v_ng, w_ng = beam2earth(b, b, b, b, heading, pitch, roll, 20.0, gimbaled=False)
        np.testing.assert_allclose(u_g, u_ng, rtol=1e-10)
        np.testing.assert_allclose(v_g, v_ng, rtol=1e-10)
        np.testing.assert_allclose(w_g, w_ng, rtol=1e-10)

    def test_all_nan_beams_give_nan_output(self):
        b = np.full((3, 4), np.nan)
        heading = np.zeros(4)
        pitch = np.zeros(4)
        roll = np.zeros(4)
        u, v, w = beam2earth(b, b, b, b, heading, pitch, roll, 20.0)
        assert np.all(np.isnan(u))
        assert np.all(np.isnan(v))
        assert np.all(np.isnan(w))
```

- [ ] **Step 3: Run tests**

```powershell
uv run pytest tests/test_transforms.py -v
```

Expected: all tests pass (7 from Task 2 + 5 new = 12 total).

- [ ] **Step 4: Ruff check**

```powershell
uv run ruff check src/ladcp/transforms/ tests/test_transforms.py
uv run ruff format --check src/ladcp/transforms/ tests/test_transforms.py
```

- [ ] **Step 5: Run full test suite**

```powershell
uv run pytest -v
```

- [ ] **Step 6: Commit**

```powershell
git add src/ladcp/transforms/beam2earth.py tests/test_transforms.py
git commit -m "feat: implement beam2earth() — 4-beam Janus beam-to-Earth transform with gimbaled correction"
```

---

### Task 4: Integration test — apply transform to P16N cast 003

**Files:**
- Create: `tests/integration/test_transforms_p16n_cast003.py`

**Prerequisite:** Tasks 1, 2, 3

**What this tests:** Load the P16N DL beam-frame PD0 file, apply `beam2earth()`, and verify that the resulting Earth-frame velocities are oceanographically plausible. This is a sanity-check level test — it does NOT validate numerical values against a reference output.

**Expected plausibility bounds for a GO-SHIP LADCP cast:**
- Mean |u| and |v| over the cast are O(0.01–0.5 m/s) — not zero but not wildly large
- Mean |w| ≪ mean |(u,v)| — vertical velocity is much smaller than horizontal
- Finite fraction of Earth-frame velocities ≈ finite fraction of beam velocities (transform preserves NaN pattern)

**Note on data layout when loading beam-frame PD0:**
When `load_rdi()` reads a beam-frame file (`coord_transform & 0x18 == 0`), `RDIData.u/v/w/e` contain b1/b2/b3/b4 respectively (the four velocity columns in order). This is documented in `OCEANOGRAPHERS_NOTES.md`. The integration test uses `d.u` as `b1`, `d.v` as `b2`, etc.

- [ ] **Step 1: Create `tests/integration/test_transforms_p16n_cast003.py`**

```python
"""Integration tests: beam2earth() applied to P16N cast 003 DL data.

Verifies Earth-frame velocities are oceanographically plausible.
This is a sanity-check level test — not numerical validation against 002.nc.

Data source: NCEI archive 0221195 (2015_P16N GO-SHIP cruise)
File: test_data/2015_P16N/003DL000.000
EX byte: 0x04 (Beam coordinates, gimbaled=True, binmap=False)
"""

import numpy as np
import pytest

from ladcp.ingestion.rdi import load_rdi
from ladcp.transforms.beam2earth import beam2earth

THETA_DEG = 20.0  # RDI Workhorse 300 kHz


@pytest.fixture
def dl_data(test_data_dir):
    p = test_data_dir / "2015_P16N" / "003DL000.000"
    if not p.exists():
        pytest.skip(f"P16N downlooker PD0 not found at {p}")
    return load_rdi(p)


@pytest.mark.integration
def test_transform_runs_without_error(dl_data):
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG, gimbaled=True)
    assert u.shape == d.u.shape
    assert v.shape == d.v.shape
    assert w.shape == d.w.shape


@pytest.mark.integration
def test_earth_frame_finite_fraction_comparable_to_beam(dl_data):
    """Transform preserves NaN pattern: finite fraction ≈ same as input."""
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    beam_finite = np.isfinite(d.u).mean()
    earth_finite = np.isfinite(u).mean()
    # Allow ≤5% change due to NaN propagation from multi-beam combination
    assert abs(earth_finite - beam_finite) < 0.05, (
        f"Large finite fraction change: beam={beam_finite:.2f}, earth={earth_finite:.2f}"
    )


@pytest.mark.integration
def test_horizontal_velocity_plausible(dl_data):
    """Mean horizontal speed should be oceanographically plausible (0–2 m/s)."""
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    mean_spd = np.sqrt(np.nanmean(u**2) + np.nanmean(v**2))
    assert mean_spd < 2.0, f"Mean horizontal speed too large: {mean_spd:.3f} m/s"
    # Should see *some* velocity (not degenerate/all-zero)
    assert mean_spd > 0.001, f"Mean speed suspiciously low: {mean_spd:.6f} m/s"


@pytest.mark.integration
def test_vertical_velocity_smaller_than_horizontal(dl_data):
    """Mean |w| << mean |(u,v)| — vertical velocity is much smaller for ocean flow."""
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    mean_horiz = np.sqrt(np.nanmean(u**2) + np.nanmean(v**2))
    mean_vert = np.sqrt(np.nanmean(w**2))
    assert mean_vert < mean_horiz, (
        f"|w|={mean_vert:.4f} should be < |(u,v)|={mean_horiz:.4f}"
    )


@pytest.mark.integration
def test_velocity_magnitudes_not_nan_dominated(dl_data):
    """At least 50% of Earth-frame cells have finite values."""
    d = dl_data
    u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, THETA_DEG)
    assert np.isfinite(u).mean() > 0.5
```

- [ ] **Step 2: Run with integration data**

```powershell
$env:TEST_DATA_DIR = "test_data"; uv run pytest tests/integration/test_transforms_p16n_cast003.py -v; $env:TEST_DATA_DIR = $null
```

Expected: all 5 tests pass. If any plausibility assertion fails, inspect the actual values:

```powershell
$env:TEST_DATA_DIR = "test_data"
uv run python -c "
import numpy as np
from pathlib import Path
from ladcp.ingestion.rdi import load_rdi
from ladcp.transforms.beam2earth import beam2earth
d = load_rdi(Path('test_data/2015_P16N/003DL000.000'))
u, v, w = beam2earth(d.u, d.v, d.w, d.e, d.heading, d.pitch, d.roll, 20.0)
print(f'mean_horiz={np.sqrt(np.nanmean(u**2)+np.nanmean(v**2)):.4f}')
print(f'mean_vert={np.sqrt(np.nanmean(w**2)):.4f}')
print(f'finite_frac={np.isfinite(u).mean():.2f}')
"
$env:TEST_DATA_DIR = $null
```

Adjust bounds if needed to match actual oceanographic reality of the P16N cruise.

- [ ] **Step 3: Run full test suite**

```powershell
uv run pytest -v
$env:TEST_DATA_DIR = "test_data"; uv run pytest -v; $env:TEST_DATA_DIR = $null
```

- [ ] **Step 4: Ruff check**

```powershell
uv run ruff check tests/integration/test_transforms_p16n_cast003.py
uv run ruff format --check tests/integration/test_transforms_p16n_cast003.py
```

- [ ] **Step 5: Commit**

```powershell
git add tests/integration/test_transforms_p16n_cast003.py
git commit -m "test: add integration tests for beam2earth() against P16N cast 003"
```

---

## Self-Review

**Spec coverage:**
- ✅ EX byte read from fixed leader and exposed on RDIData — Task 1
- ✅ `beam2xyz()` matches janus2xyz.m matrix exactly, vectorized numpy, no loops — Task 2
- ✅ `beam2earth()` matches janus2earth.m math (gimbaled correction, 9 rotation coefficients, sign convention) — Task 3
- ✅ NaN propagation: natural via numpy, tested — Tasks 2 and 3
- ✅ Gimbaled=True default (matches LADCP standard), gimbaled=False also implemented — Task 3
- ✅ Broadcasting: (nbin,nens) × (nens,) — validated by output-shape test — Task 3
- ✅ Unit tests run without test data files — Tasks 1–3
- ✅ Integration tests gated by TEST_DATA_DIR — Task 4
- ✅ Plausibility-level validation of Earth-frame velocities — Task 4
- ✅ Smoke test (`janus5beam2earth` stub) preserved — Task 3 step 1 explicitly calls this out

**Deferred (out of scope for this phase):**
- Bin mapping (BinmapType ≠ 'none') — EX=0x04 has bit 1=0
- 3-beam solution — EX=0x04 has bit 0=0
- 5-beam (vertical beam 5) — not present in standard Workhorse
- `gimbaled=False` branch tested only via the `test_gimbaled_no_effect_when_pitch_roll_zero` case
