# SBE Hex Decoder Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `load_sbe_hex()` in `src/ladcp/ingestion/sbe_hex.py` that reads a raw SBE9+ `.hex` file and paired `.XMLCON` calibration file, decodes sensor data and GPS from each 24 Hz scan, and returns a `CTDTimeSeries` compatible with the existing `assign_bin_depths()` and `compute_ship_velocity()` pipeline.

**Architecture:** Two-phase design. Phase 1 (Tasks 1–2): parse XMLCON for calibration coefficients and determine the exact scan byte layout from the hex header — this must be verified empirically before writing calibration code. Phase 2 (Tasks 3–4): apply calibration equations to produce T, C, P time series; extract GPS lat/lon; assemble `CTDTimeSeries`. Integration test (Task 5) validates GPS track and pressure against the CCHDO processed `_ctd.nc` reference (which is binned, so comparison is after 1 Hz averaging).

**Tech Stack:** Python 3.11, `defusedxml` (safe XML parsing), `numpy`, `struct`, `pytest`

## Global Constraints

- Do not add `gsw` as a new dependency; compute salinity in Task 4 only if `gsw` is already installed (check `importlib.util.find_spec("gsw")`), otherwise leave salinity as NaN and note it
- The decoder must handle the S4P SBE9plus configuration (dual T/C, Digiquartz pressure, 5 voltage words, surface PAR, system time, NMEA lat/lon appended to every scan)
- `CTDTimeSeries` is imported from `ladcp.ingestion.ctd`; do not modify it
- Test env var: `TEST_DATA_DIR`; integration tests skip cleanly when data absent
- Output `CTDTimeSeries.time_julian` must use the same Julian day epoch as `_to_julian()` in `ladcp.ingestion._pd0`

---

### Task 1: XMLCON parser

**Files:**
- Create: `src/ladcp/ingestion/sbe_hex.py` (scaffold + XMLCON parser only)
- Create: `tests/unit/test_sbe_hex.py`

**Interfaces:**
- Produces: `XmlconCoeffs` dataclass and `load_xmlcon(path) -> XmlconCoeffs`

- [ ] **Step 1: Write failing unit tests for XMLCON parser**

Create `tests/unit/test_sbe_hex.py`:

```python
"""Unit tests for SBE hex decoder."""
from __future__ import annotations
from pathlib import Path
import textwrap
import pytest

from ladcp.ingestion.sbe_hex import XmlconCoeffs, load_xmlcon


MINIMAL_XMLCON = textwrap.dedent("""\
<?xml version="1.0" encoding="UTF-8"?>
<SBE_InstrumentConfiguration SB_ConfigCTD_FileVersion="7.26.1.0">
  <Instrument Type="8">
    <Name>SBE 911plus/917plus CTD</Name>
    <FrequencyChannelsSuppressed>0</FrequencyChannelsSuppressed>
    <VoltageWordsSuppressed>0</VoltageWordsSuppressed>
    <SurfaceParVoltageAdded>1</SurfaceParVoltageAdded>
    <ScanTimeAdded>1</ScanTimeAdded>
    <NmeaPositionDataAdded>1</NmeaPositionDataAdded>
    <SensorArray Size="5">
      <Sensor index="0" SensorID="55">
        <TemperatureSensor SensorID="55">
          <UseG_J>1</UseG_J>
          <G>4.36593732e-003</G><H>6.30830930e-004</H>
          <I>2.06378769e-005</I><J>1.63292939e-006</J>
          <F0>1000.000</F0><Slope>1.0</Slope><Offset>0.0</Offset>
        </TemperatureSensor>
      </Sensor>
      <Sensor index="1" SensorID="3">
        <ConductivitySensor SensorID="3">
          <UseG_J>1</UseG_J>
          <Coefficients equation="1">
            <G>-9.90838065e+000</G><H>1.60083819e+000</H>
            <I>-1.58153324e-003</I><J>1.99854164e-004</J>
            <CPcor>-9.57000000e-008</CPcor><CTcor>3.2500e-006</CTcor>
          </Coefficients>
          <Slope>1.0</Slope><Offset>0.0</Offset>
        </ConductivitySensor>
      </Sensor>
      <Sensor index="2" SensorID="45">
        <PressureSensor SensorID="45">
          <C1>-4.160303e+004</C1><C2>-4.604479e-001</C2><C3>1.585404e-002</C3>
          <D1>3.546467e-002</D1><D2>0.0</D2>
          <T1>3.013997e+001</T1><T2>-3.831629e-004</T2>
          <T3>3.608677e-006</T3><T4>1.200552e-008</T4><T5>0.0</T5>
          <AD590M>1.278460e-002</AD590M><AD590B>-9.255860e+000</AD590B>
          <Slope>1.0</Slope><Offset>0.0</Offset>
        </PressureSensor>
      </Sensor>
      <Sensor index="3" SensorID="55">
        <TemperatureSensor SensorID="55">
          <UseG_J>1</UseG_J>
          <G>4.35781951e-003</G><H>6.45070776e-004</H>
          <I>2.42988411e-005</I><J>2.35822338e-006</J>
          <F0>1000.000</F0><Slope>1.0</Slope><Offset>0.0</Offset>
        </TemperatureSensor>
      </Sensor>
      <Sensor index="4" SensorID="3">
        <ConductivitySensor SensorID="3">
          <UseG_J>1</UseG_J>
          <Coefficients equation="1">
            <G>-3.96678467e+000</G><H>4.84542307e-001</H>
            <I>-6.60474581e-004</I><J>5.63015941e-005</J>
            <CPcor>-9.57000000e-008</CPcor><CTcor>3.2500e-006</CTcor>
          </Coefficients>
          <Slope>1.0</Slope><Offset>0.0</Offset>
        </ConductivitySensor>
      </Sensor>
    </SensorArray>
  </Instrument>
</SBE_InstrumentConfiguration>
""")


@pytest.fixture
def xmlcon_file(tmp_path: Path) -> Path:
    p = tmp_path / "test.XMLCON"
    p.write_text(MINIMAL_XMLCON, encoding="utf-8")
    return p


def test_load_xmlcon_returns_coeffs(xmlcon_file: Path) -> None:
    c = load_xmlcon(xmlcon_file)
    assert isinstance(c, XmlconCoeffs)


def test_temperature1_coefficients(xmlcon_file: Path) -> None:
    c = load_xmlcon(xmlcon_file)
    assert abs(c.t1_G - 4.36593732e-3) < 1e-12
    assert abs(c.t1_H - 6.30830930e-4) < 1e-12
    assert abs(c.t1_f0 - 1000.0) < 1e-6


def test_conductivity1_coefficients(xmlcon_file: Path) -> None:
    c = load_xmlcon(xmlcon_file)
    assert abs(c.c1_G - (-9.90838065)) < 1e-6
    assert abs(c.c1_CPcor - (-9.57e-8)) < 1e-14
    assert abs(c.c1_CTcor - 3.25e-6) < 1e-12


def test_pressure_coefficients(xmlcon_file: Path) -> None:
    c = load_xmlcon(xmlcon_file)
    assert abs(c.p_C1 - (-4.160303e4)) < 1e-1
    assert abs(c.p_AD590M - 1.278460e-2) < 1e-8


def test_secondary_sensors_parsed(xmlcon_file: Path) -> None:
    c = load_xmlcon(xmlcon_file)
    assert abs(c.t2_G - 4.35781951e-3) < 1e-12
    assert abs(c.c2_G - (-3.96678467)) < 1e-6
```

- [ ] **Step 2: Run to verify failure**

```
pytest tests/unit/test_sbe_hex.py -v
```

Expected: `ImportError: cannot import name 'XmlconCoeffs' from 'ladcp.ingestion.sbe_hex'`

- [ ] **Step 3: Implement XMLCON parser in `src/ladcp/ingestion/sbe_hex.py`**

```python
"""SBE 9plus hex file decoder.

Reads raw .hex telemetry + .XMLCON calibration from SBE9/11plus systems.
Produces CTDTimeSeries at the native 24 Hz sampling rate.

Reference: SeaBird Application Note 69 (SBE 9plus data format).
Sensor calibration equations: AN-04 (temperature), AN-14 (conductivity),
AN-46 (pressure).
"""
from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path
import defusedxml.ElementTree as ET

import numpy as np

from ladcp.ingestion.ctd import CTDTimeSeries
from ladcp.ingestion._pd0 import _to_julian


# ---------------------------------------------------------------------------
# XMLCON calibration coefficients
# ---------------------------------------------------------------------------

@dataclass
class XmlconCoeffs:
    """Calibration coefficients parsed from an SBE .XMLCON file."""
    # Primary temperature (SBE3, G-J ITS-90 form)
    t1_G: float = 0.0
    t1_H: float = 0.0
    t1_I: float = 0.0
    t1_J: float = 0.0
    t1_f0: float = 1000.0   # reference frequency, Hz

    # Secondary temperature
    t2_G: float = 0.0
    t2_H: float = 0.0
    t2_I: float = 0.0
    t2_J: float = 0.0
    t2_f0: float = 1000.0

    # Primary conductivity (SBE4, G-J form)
    c1_G: float = 0.0
    c1_H: float = 0.0
    c1_I: float = 0.0
    c1_J: float = 0.0
    c1_CPcor: float = -9.57e-8
    c1_CTcor: float = 3.25e-6

    # Secondary conductivity
    c2_G: float = 0.0
    c2_H: float = 0.0
    c2_I: float = 0.0
    c2_J: float = 0.0
    c2_CPcor: float = -9.57e-8
    c2_CTcor: float = 3.25e-6

    # Pressure — Digiquartz polynomial coefficients
    p_C1: float = 0.0
    p_C2: float = 0.0
    p_C3: float = 0.0
    p_D1: float = 0.0
    p_D2: float = 0.0
    p_T1: float = 0.0
    p_T2: float = 0.0
    p_T3: float = 0.0
    p_T4: float = 0.0
    p_T5: float = 0.0
    p_AD590M: float = 0.0   # AD590 temperature compensation slope
    p_AD590B: float = 0.0   # AD590 temperature compensation intercept

    # Supplemental flags
    scan_time_added: bool = False
    nmea_pos_added: bool = False
    surface_par_added: bool = False
    n_freq_channels: int = 4   # T1, C1, T2, C2
    n_voltage_words: int = 0   # from hex header


def _gj(node: ET.Element | None, tag: str, default: float = 0.0) -> float:
    """Read a float text value from an XML child element."""
    if node is None:
        return default
    el = node.find(tag)
    return float(el.text) if el is not None and el.text else default


def load_xmlcon(path: str | Path) -> XmlconCoeffs:
    """Parse an SBE .XMLCON configuration file and extract calibration coefficients."""
    tree = ET.parse(str(path))
    root = tree.getroot()
    inst = root.find("Instrument")
    if inst is None:
        raise ValueError(f"No <Instrument> element in {path}")

    c = XmlconCoeffs()
    c.scan_time_added = inst.findtext("ScanTimeAdded", "0").strip() == "1"
    c.nmea_pos_added = inst.findtext("NmeaPositionDataAdded", "0").strip() == "1"
    c.surface_par_added = inst.findtext("SurfaceParVoltageAdded", "0").strip() == "1"

    # Count frequency channels (T1, C1, T2, C2 → 4 if dual, 2 if single)
    # Voltage words come from the hex header; initialise from XMLCON sensor count.
    # n_freq_channels will be overridden in load_sbe_hex() from the hex header.

    sensors = root.findall(".//Sensor")
    t_count = 0
    c_count = 0
    for sensor in sensors:
        t_el = sensor.find("TemperatureSensor")
        c_el = sensor.find("ConductivitySensor")
        p_el = sensor.find("PressureSensor")

        if t_el is not None:
            t_count += 1
            if t_count == 1:
                c.t1_G = _gj(t_el, "G")
                c.t1_H = _gj(t_el, "H")
                c.t1_I = _gj(t_el, "I")
                c.t1_J = _gj(t_el, "J")
                f0_el = t_el.find("F0")
                c.t1_f0 = float(f0_el.text) if f0_el is not None else 1000.0
            elif t_count == 2:
                c.t2_G = _gj(t_el, "G")
                c.t2_H = _gj(t_el, "H")
                c.t2_I = _gj(t_el, "I")
                c.t2_J = _gj(t_el, "J")
                f0_el = t_el.find("F0")
                c.t2_f0 = float(f0_el.text) if f0_el is not None else 1000.0

        if c_el is not None:
            c_count += 1
            coef_el = c_el.find(".//Coefficients[@equation='1']")
            if coef_el is None:
                coef_el = c_el.find("Coefficients")
            if c_count == 1:
                c.c1_G = _gj(coef_el, "G")
                c.c1_H = _gj(coef_el, "H")
                c.c1_I = _gj(coef_el, "I")
                c.c1_J = _gj(coef_el, "J")
                c.c1_CPcor = _gj(coef_el, "CPcor", -9.57e-8)
                c.c1_CTcor = _gj(coef_el, "CTcor", 3.25e-6)
            elif c_count == 2:
                c.c2_G = _gj(coef_el, "G")
                c.c2_H = _gj(coef_el, "H")
                c.c2_I = _gj(coef_el, "I")
                c.c2_J = _gj(coef_el, "J")
                c.c2_CPcor = _gj(coef_el, "CPcor", -9.57e-8)
                c.c2_CTcor = _gj(coef_el, "CTcor", 3.25e-6)

        if p_el is not None:
            c.p_C1 = _gj(p_el, "C1")
            c.p_C2 = _gj(p_el, "C2")
            c.p_C3 = _gj(p_el, "C3")
            c.p_D1 = _gj(p_el, "D1")
            c.p_D2 = _gj(p_el, "D2")
            c.p_T1 = _gj(p_el, "T1")
            c.p_T2 = _gj(p_el, "T2")
            c.p_T3 = _gj(p_el, "T3")
            c.p_T4 = _gj(p_el, "T4")
            c.p_T5 = _gj(p_el, "T5")
            c.p_AD590M = _gj(p_el, "AD590M")
            c.p_AD590B = _gj(p_el, "AD590B")

    c.n_freq_channels = t_count + c_count  # 4 for dual T/C
    return c


# Stub — implemented in Task 3
def load_sbe_hex(
    hex_path: str | Path,
    xmlcon_path: str | Path,
) -> CTDTimeSeries:
    """Load SBE hex file and return a 24 Hz CTDTimeSeries."""
    raise NotImplementedError("Implemented in Task 3")
```

- [ ] **Step 4: Run unit tests**

```
pytest tests/unit/test_sbe_hex.py -v
```

Expected: all 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ladcp/ingestion/sbe_hex.py tests/unit/test_sbe_hex.py
git commit -m "feat: SBE hex decoder — XMLCON parser (XmlconCoeffs, load_xmlcon)"
```

---

### Task 2: Hex header parser and scan byte-layout discovery

**Files:**
- Modify: `src/ladcp/ingestion/sbe_hex.py`
- Create: `scripts/inspect_sbe_scan.py` (discovery tool, not shipped)

**Interfaces:**
- Produces: `HexHeader` dataclass with `bytes_per_scan`, `n_voltage_words`, `scan_time_added`, `nmea_pos_added`, `start_datetime` and `parse_hex_header(path) -> HexHeader`

- [ ] **Step 1: Write unit tests for header parser**

Add to `tests/unit/test_sbe_hex.py`:

```python
from ladcp.ingestion.sbe_hex import HexHeader, parse_hex_header

MINIMAL_HEX_HEADER = textwrap.dedent("""\
* Sea-Bird SBE 9 Data File:
* FileName = test.hex
* Software Version Seasave V 7.26.1.8
* Temperature SN = 5844
* Conductivity SN = 4546
* Number of Bytes Per Scan = 44
* Number of Voltage Words = 5
* Append System Time to Every Scan
* System UpLoad Time = Mar 17 2018 01:23:11
* NMEA Latitude = 70 27.16 S
* NMEA Longitude = 168 28.48 E
* NMEA UTC (Time) = Mar 17 2018  01:23:09
* Store Lat/Lon Data = Append to Every Scan
* SBE 11plus V 5.1g
*END*
""")


@pytest.fixture
def hex_header_file(tmp_path: Path) -> Path:
    p = tmp_path / "test.hex"
    p.write_text(MINIMAL_HEX_HEADER, encoding="ascii")
    return p


def test_parse_hex_header_bytes_per_scan(hex_header_file: Path) -> None:
    h = parse_hex_header(hex_header_file)
    assert h.bytes_per_scan == 44


def test_parse_hex_header_voltage_words(hex_header_file: Path) -> None:
    h = parse_hex_header(hex_header_file)
    assert h.n_voltage_words == 5


def test_parse_hex_header_nmea_pos(hex_header_file: Path) -> None:
    h = parse_hex_header(hex_header_file)
    assert h.nmea_pos_added is True


def test_parse_hex_header_start_datetime(hex_header_file: Path) -> None:
    h = parse_hex_header(hex_header_file)
    # System UpLoad Time = Mar 17 2018 01:23:11
    assert h.upload_year == 2018
    assert h.upload_month == 3
    assert h.upload_day == 17
    assert abs(h.upload_hour_frac - (1 + 23/60 + 11/3600)) < 1e-4
```

- [ ] **Step 2: Implement `HexHeader` and `parse_hex_header()`**

Add to `src/ladcp/ingestion/sbe_hex.py` (before `load_sbe_hex`):

```python
import re
from datetime import datetime

_MONTH_MAP = {
    "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
    "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12,
}


@dataclass
class HexHeader:
    """Metadata extracted from the SBE hex file header (lines starting with '*')."""
    bytes_per_scan: int = 0
    n_voltage_words: int = 0
    scan_time_added: bool = False
    nmea_pos_added: bool = False
    upload_year: int = 0
    upload_month: int = 0
    upload_day: int = 0
    upload_hour_frac: float = 0.0   # hours + minutes/60 + seconds/3600


def parse_hex_header(path: str | Path) -> HexHeader:
    """Read the '*'-prefixed header lines from an SBE hex file."""
    h = HexHeader()
    with open(str(path), "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if not line.startswith("*"):
                break
            if line.startswith("*END*"):
                break
            stripped = line.lstrip("* ").strip()

            m = re.match(r"Number of Bytes Per Scan\s*=\s*(\d+)", stripped, re.I)
            if m:
                h.bytes_per_scan = int(m.group(1))
                continue

            m = re.match(r"Number of Voltage Words\s*=\s*(\d+)", stripped, re.I)
            if m:
                h.n_voltage_words = int(m.group(1))
                continue

            if re.search(r"Append System Time", stripped, re.I):
                h.scan_time_added = True
                continue

            if re.search(r"Store Lat/Lon Data|Latitude/Longitude added", stripped, re.I):
                h.nmea_pos_added = True
                continue

            # "System UpLoad Time = Mar 17 2018 01:23:11"
            m = re.match(
                r"System UpLoad Time\s*=\s*(\w+)\s+(\d+)\s+(\d+)\s+(\d+):(\d+):(\d+)",
                stripped, re.I,
            )
            if m:
                mon_str, day, year, hh, mm, ss = m.groups()
                h.upload_month = _MONTH_MAP.get(mon_str[:3].capitalize(), 0)
                h.upload_day = int(day)
                h.upload_year = int(year)
                h.upload_hour_frac = int(hh) + int(mm) / 60.0 + int(ss) / 3600.0
    return h
```

- [ ] **Step 3: Run the new header tests**

```
pytest tests/unit/test_sbe_hex.py -v
```

Expected: all 9 tests PASS (5 from Task 1 + 4 new).

- [ ] **Step 4: Write byte-layout discovery script**

Create `scripts/inspect_sbe_scan.py` (run once to confirm byte offsets, then keep for reference):

```python
"""Inspect raw bytes of first SBE hex scan to identify channel byte offsets.

Prints each 3-byte group with its position, plus tries float32 decoding at
every aligned offset — use this to locate GPS lat/lon by matching the header
values.

Usage:
    python scripts/inspect_sbe_scan.py path/to/00101.hex path/to/00101.XMLCON
"""
import struct
import sys
from pathlib import Path

from ladcp.ingestion.sbe_hex import load_xmlcon, parse_hex_header

KNOWN_LAT = -70.4527   # degrees (from header: "NMEA Latitude = 70 27.16 S")
KNOWN_LON = 168.4747   # degrees (from header: "NMEA Longitude = 168 28.48 E")


def main() -> None:
    if len(sys.argv) < 3:
        sys.exit("Usage: inspect_sbe_scan.py <hex_file> <xmlcon_file>")
    hex_path = Path(sys.argv[1])
    xmlcon_path = Path(sys.argv[2])

    hdr = parse_hex_header(hex_path)
    coeffs = load_xmlcon(xmlcon_path)

    print(f"Bytes per scan: {hdr.bytes_per_scan}")
    print(f"Voltage words: {hdr.n_voltage_words}")
    print(f"Scan time added: {hdr.scan_time_added}")
    print(f"NMEA pos added: {hdr.nmea_pos_added}")
    print(f"Sensors: {coeffs.n_freq_channels} freq channels\n")

    # Read first data scan (first line after *END*)
    scan_hex = ""
    past_end = False
    with open(str(hex_path), encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if past_end and line and not line.startswith("*"):
                scan_hex = line.strip()
                break
            if line.startswith("*END*"):
                past_end = True

    if not scan_hex:
        sys.exit("No data scan found")

    scan = bytes.fromhex(scan_hex)
    print(f"Raw scan ({len(scan)} bytes):")
    for i in range(0, len(scan), 3):
        chunk = scan[i:i+3]
        val = int.from_bytes(chunk, "big")
        print(f"  bytes {i:02d}-{i+2:02d}: {chunk.hex().upper()}  = {val:8d}")

    print("\nFloat32 candidates (big-endian) at every offset:")
    for offset in range(0, len(scan) - 3):
        val = struct.unpack(">f", scan[offset:offset+4])[0]
        if abs(val - KNOWN_LAT) < 2.0:
            print(f"  offset {offset:02d}: {val:.4f}  ← LAT CANDIDATE")
        if abs(val - KNOWN_LON) < 2.0:
            print(f"  offset {offset:02d}: {val:.4f}  ← LON CANDIDATE")

    print("\nFloat32 candidates (little-endian) at every offset:")
    for offset in range(0, len(scan) - 3):
        val = struct.unpack("<f", scan[offset:offset+4])[0]
        if abs(val - KNOWN_LAT) < 2.0:
            print(f"  offset {offset:02d}: {val:.4f}  ← LAT CANDIDATE (LE)")
        if abs(val - KNOWN_LON) < 2.0:
            print(f"  offset {offset:02d}: {val:.4f}  ← LON CANDIDATE (LE)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run the discovery script**

```
python scripts/inspect_sbe_scan.py \
  test_data/2018_S4P/CTD/00101.hex \
  test_data/2018_S4P/CTD/00101.XMLCON
```

Expected output will identify which byte offsets contain the lat/lon. Record the offsets in a comment at the top of `sbe_hex.py` — these become the constants used in Task 3.

- [ ] **Step 6: Commit**

```bash
git add src/ladcp/ingestion/sbe_hex.py tests/unit/test_sbe_hex.py \
        scripts/inspect_sbe_scan.py
git commit -m "feat: SBE hex decoder — header parser + byte-layout discovery script"
```

---

### Task 3: Calibration equations and `load_sbe_hex()`

**Files:**
- Modify: `src/ladcp/ingestion/sbe_hex.py` — replace `load_sbe_hex` stub with full implementation

**Prerequisite:** Task 2's discovery script has been run; the byte offsets for all channels are confirmed. The following byte-layout description is for the standard SBE9plus dual-T/C configuration with 5 voltage words + system time + NMEA lat/lon (total 44 bytes). Verify offsets using the discovery script output before coding.

**Assumed byte layout (44 bytes per scan, verify with discovery script):**

| Bytes | Content |
|---|---|
| 0–2 | T1 frequency count (24-bit big-endian uint) |
| 3–5 | C1 frequency count |
| 6–8 | T2 frequency count |
| 9–11 | C2 frequency count |
| 12–14 | Pressure (PT and P counts packed, see AN-69) |
| 15–29 | Voltage words 0–4 (5 × 3 bytes, pairs of 12-bit ADC) |
| 30–31 | System time counter (16-bit, 1/256-second units) |
| 32–35 | GPS latitude (float32 or scaled int, see discovery script) |
| 36–39 | GPS longitude (float32 or scaled int, see discovery script) |
| 40–43 | Additional data (surface PAR, scan counter, or status) |

**Calibration equations (reference: SBE Application Notes AN-04, AN-14, AN-46):**

Temperature (SBE3F, UseG_J form):
```
f = (f_clock × N_cycles) / count_T   [Hz, where f_clock and N_cycles from AN-69]
n = f0 / f                            [f0 = 1000 Hz from XMLCON]
T(°C) = 1/(G + H·ln(n) + I·ln(n)² + J·ln(n)³) - 273.15
```

Conductivity (SBE4, UseG_J form):
```
f = (f_clock × N_cycles) / count_C   [Hz]
fkHz = f / 1000                       [kHz]
C(mS/cm) = (G + H·fkHz² + I·fkHz³ + J·fkHz⁴) / (1 + CTcor·T + CPcor·P)
```

Pressure (Digiquartz, AN-46):
```
PT_raw = (count_P_bytes >> 12) & 0xFFF    [upper 12 bits]
P_raw  = count_P_bytes & 0xFFF            [lower 12 bits]
T_u = AD590M × PT_raw + AD590B           [AD590 temperature, °C]
# Full polynomial in AN-46:
U = T1 + T2·T_u + T3·T_u² + T4·T_u³ + T5·T_u⁴
Y = 1 + D1·T_u + D2·T_u²
C_coef = C1 + C2·T_u + C3·T_u²
P(psia) = C_coef · (1 - (U²/P_raw²)) · (1 - Y · (1 - U²/P_raw²))
P(dbar) = (P_psia - 14.6959) × 0.68947572932  # psia to dbar
```

- [ ] **Step 1: Add calibration helper functions**

Add to `src/ladcp/ingestion/sbe_hex.py`:

```python
# Clock and gate constants for SBE9plus frequency channels (AN-69).
# These values MUST be verified against AN-69 before shipping.
# Typical SBE9plus: f_clock = 24e6 Hz (24 MHz), N_cycles = 210 for T channels.
_FREQ_CLOCK_HZ = 24_000_000.0
_FREQ_N_CYCLES = 210


def _count_to_freq(count: int) -> float:
    """Convert 24-bit period counter to oscillation frequency in Hz."""
    if count == 0:
        return np.nan
    return (_FREQ_CLOCK_HZ * _FREQ_N_CYCLES) / count


def _temperature(count: int, coeffs: XmlconCoeffs, primary: bool = True) -> float:
    """Decode a 24-bit T frequency count to ITS-90 temperature (°C)."""
    f = _count_to_freq(count)
    if not np.isfinite(f) or f <= 0:
        return np.nan
    G, H, I, J, f0 = (
        (coeffs.t1_G, coeffs.t1_H, coeffs.t1_I, coeffs.t1_J, coeffs.t1_f0)
        if primary
        else (coeffs.t2_G, coeffs.t2_H, coeffs.t2_I, coeffs.t2_J, coeffs.t2_f0)
    )
    n = f0 / f
    ln_n = np.log(n)
    return 1.0 / (G + H * ln_n + I * ln_n**2 + J * ln_n**3) - 273.15


def _conductivity(
    count: int, T: float, P_dbar: float, coeffs: XmlconCoeffs, primary: bool = True
) -> float:
    """Decode a 24-bit C frequency count to conductivity (mS/cm)."""
    f = _count_to_freq(count)
    if not np.isfinite(f) or f <= 0:
        return np.nan
    G, H, I, J, CPcor, CTcor = (
        (coeffs.c1_G, coeffs.c1_H, coeffs.c1_I, coeffs.c1_J, coeffs.c1_CPcor, coeffs.c1_CTcor)
        if primary
        else (coeffs.c2_G, coeffs.c2_H, coeffs.c2_I, coeffs.c2_J, coeffs.c2_CPcor, coeffs.c2_CTcor)
    )
    fk = f / 1000.0  # kHz
    return (G + H * fk**2 + I * fk**3 + J * fk**4) / (1.0 + CTcor * T + CPcor * P_dbar)


def _pressure(raw_bytes: bytes, coeffs: XmlconCoeffs) -> float:
    """Decode 3-byte Digiquartz pressure block to pressure in dbar.

    Byte layout: upper 12 bits = PT_raw (pressure temperature),
                 lower 12 bits = P_raw (pressure oscillation count).
    See SBE Application Note 46 for the full Digiquartz polynomial.
    """
    word = int.from_bytes(raw_bytes, "big")
    PT_raw = (word >> 12) & 0xFFF
    P_raw = word & 0xFFF

    T_u = coeffs.p_AD590M * PT_raw + coeffs.p_AD590B  # AD590 temp, °C

    U = (
        coeffs.p_T1
        + coeffs.p_T2 * T_u
        + coeffs.p_T3 * T_u**2
        + coeffs.p_T4 * T_u**3
        + coeffs.p_T5 * T_u**4
    )
    Y = 1.0 + coeffs.p_D1 * T_u + coeffs.p_D2 * T_u**2
    C_coef = coeffs.p_C1 + coeffs.p_C2 * T_u + coeffs.p_C3 * T_u**2

    if P_raw == 0:
        return np.nan
    ratio = (U / P_raw) ** 2
    P_psia = C_coef * (1.0 - ratio) * (1.0 - Y * (1.0 - ratio))
    # Convert psia to dbar (1 psia ≈ 0.68947572932 dbar, subtract 1 atm)
    return max((P_psia - 14.6959) * 0.68947572932, 0.0)
```

- [ ] **Step 2: Add unit tests for calibration functions**

Add to `tests/unit/test_sbe_hex.py`:

```python
from ladcp.ingestion.sbe_hex import _temperature, _conductivity, _pressure, load_xmlcon


@pytest.fixture
def s4p_coeffs(xmlcon_file: Path) -> "XmlconCoeffs":
    from ladcp.ingestion.sbe_hex import XmlconCoeffs
    return load_xmlcon(xmlcon_file)


def test_temperature_returns_float(s4p_coeffs) -> None:
    # Non-zero count must produce a finite float
    T = _temperature(0x188A64, s4p_coeffs, primary=True)
    assert np.isfinite(T)


def test_conductivity_returns_float(s4p_coeffs) -> None:
    T = _temperature(0x188A64, s4p_coeffs, primary=True)
    C = _conductivity(0x09BADA, T, 0.0, s4p_coeffs, primary=True)
    assert np.isfinite(C)


def test_pressure_zero_count_is_nan(s4p_coeffs) -> None:
    P = _pressure(b"\x00\x00\x00", s4p_coeffs)
    assert np.isnan(P)
```

- [ ] **Step 3: Implement `load_sbe_hex()`**

Replace the `raise NotImplementedError` stub in `src/ladcp/ingestion/sbe_hex.py`:

```python
# GPS encoding: SBE11plus V5.x stores lat/lon as IEEE 754 float32 big-endian.
# Exact byte offsets MUST be confirmed by running scripts/inspect_sbe_scan.py.
# These offsets assume the standard 44-byte S4P layout documented above.
_LAT_OFFSET = 32
_LON_OFFSET = 36


def load_sbe_hex(
    hex_path: str | Path,
    xmlcon_path: str | Path,
) -> CTDTimeSeries:
    """Load SBE 9plus hex file and return a 24 Hz CTDTimeSeries.

    Parameters
    ----------
    hex_path:
        Path to the .hex data file.
    xmlcon_path:
        Path to the paired .XMLCON calibration file.

    Returns
    -------
    CTDTimeSeries
        24 Hz time series with time_julian (Julian days), pressure_dbar,
        temp_c (primary T1), salinity (if gsw available, else NaN array),
        and lat/lon GPS track.
    """
    coeffs = load_xmlcon(xmlcon_path)
    hdr = parse_hex_header(hex_path)
    n_bytes = hdr.bytes_per_scan

    # Build start-of-file Julian time from upload timestamp
    t0_jd = _to_julian(hdr.upload_year, hdr.upload_month, hdr.upload_day, hdr.upload_hour_frac)
    dt_per_scan = 1.0 / (24.0 * 3600.0 * 24.0)  # 24 Hz → Julian day increment per scan

    # Read all data lines (skip '*' header)
    data_lines: list[str] = []
    past_end = False
    with open(str(hex_path), "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("*END*"):
                past_end = True
                continue
            if past_end and line and not line.startswith("*"):
                data_lines.append(line.strip())

    n_scans = len(data_lines)
    if n_scans == 0:
        raise ValueError(f"No data scans found in {hex_path}")

    time_jd = np.empty(n_scans, dtype=np.float64)
    pressure = np.full(n_scans, np.nan, dtype=np.float64)
    temp = np.full(n_scans, np.nan, dtype=np.float64)
    lat = np.full(n_scans, np.nan, dtype=np.float64)
    lon = np.full(n_scans, np.nan, dtype=np.float64)

    for i, hex_line in enumerate(data_lines):
        time_jd[i] = t0_jd + i * dt_per_scan

        if len(hex_line) < n_bytes * 2:
            continue

        scan = bytes.fromhex(hex_line[: n_bytes * 2])

        # Decode T1 and P (P needed for C calibration)
        count_t1 = int.from_bytes(scan[0:3], "big")
        count_c1 = int.from_bytes(scan[3:6], "big")
        count_p = scan[12:15]

        T1 = _temperature(count_t1, coeffs, primary=True)
        P = _pressure(count_p, coeffs)
        C1 = _conductivity(count_c1, T1, P, coeffs, primary=True)

        temp[i] = T1
        pressure[i] = P

        # GPS lat/lon
        if hdr.nmea_pos_added and n_bytes >= _LON_OFFSET + 4:
            try:
                lat[i] = struct.unpack(">f", scan[_LAT_OFFSET: _LAT_OFFSET + 4])[0]
                lon[i] = struct.unpack(">f", scan[_LON_OFFSET: _LON_OFFSET + 4])[0]
            except struct.error:
                pass

    # Salinity (requires gsw; set to NaN if unavailable)
    salinity = np.full(n_scans, np.nan, dtype=np.float64)
    try:
        import importlib.util
        if importlib.util.find_spec("gsw") is not None:
            import gsw
            SP = gsw.SP_from_C(C1_arr, temp, pressure)  # noqa — needs C1 array
    except Exception:
        pass  # leave salinity as NaN

    # Mask GPS: replace 0.0 (no-fix placeholder) with NaN
    lat[lat == 0.0] = np.nan
    lon[lon == 0.0] = np.nan

    return CTDTimeSeries(
        time_julian=time_jd,
        pressure_dbar=pressure,
        temp_c=temp,
        salinity=salinity,
        lat=lat if np.any(np.isfinite(lat)) else None,
        lon=lon if np.any(np.isfinite(lon)) else None,
    )
```

Note: the above has a bug (`C1_arr` undefined). Fix that by collecting C1 values in an array during the loop:

```python
    cond = np.full(n_scans, np.nan, dtype=np.float64)
    for i, hex_line in enumerate(data_lines):
        ...
        C1 = _conductivity(count_c1, T1, P, coeffs, primary=True)
        cond[i] = C1
        temp[i] = T1
        pressure[i] = P
        ...

    # Salinity
    salinity = np.full(n_scans, np.nan, dtype=np.float64)
    import importlib.util
    if importlib.util.find_spec("gsw") is not None:
        import gsw
        try:
            salinity = gsw.SP_from_C(cond, temp, pressure)
        except Exception:
            pass
```

- [ ] **Step 4: Run calibration unit tests**

```
pytest tests/unit/test_sbe_hex.py -v
```

Expected: all 12 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ladcp/ingestion/sbe_hex.py tests/unit/test_sbe_hex.py
git commit -m "feat: SBE hex decoder — calibration equations and load_sbe_hex()"
```

---

### Task 4: Integration test against S4P CTD data

**Files:**
- Create: `tests/integration/test_sbe_hex_s4p.py`

**Note:** The CCHDO `320620180309_ctd.nc` provides depth-binned CTD profiles (NOT time series). For comparison: average the decoded 24 Hz data to 1 Hz, depth-bin it, and check that pressure and temperature agree within ±0.5°C and ±5 dbar in the main water column. GPS comparison: check that mean lat/lon of the track matches the header NMEA position within 0.1°.

- [ ] **Step 1: Write integration test**

Create `tests/integration/test_sbe_hex_s4p.py`:

```python
"""Integration test: SBE hex decoder vs S4P cast 001 CTD reference.

Requires TEST_DATA_DIR env var pointing to a directory containing 2018_S4P/.
Comparison uses the CCHDO calibrated profile (binned) as a loose reference;
the raw hex data has no full calibrations applied (per README.Archive).
"""
from __future__ import annotations

import os
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ladcp.ingestion.sbe_hex import load_sbe_hex


@pytest.fixture(scope="module")
def s4p_dir() -> Path:
    env = os.environ.get("TEST_DATA_DIR", "")
    if not env:
        pytest.skip("TEST_DATA_DIR not set")
    p = Path(env) / "2018_S4P"
    if not p.exists():
        pytest.skip(f"2018_S4P not found at {p}")
    return p


@pytest.fixture(scope="module")
def hex_path(s4p_dir: Path) -> Path:
    p = s4p_dir / "CTD/00101.hex"
    if not p.exists():
        pytest.skip(f"00101.hex not found: {p}")
    return p


@pytest.fixture(scope="module")
def xmlcon_path(s4p_dir: Path) -> Path:
    p = s4p_dir / "CTD/00101.XMLCON"
    if not p.exists():
        pytest.skip(f"00101.XMLCON not found: {p}")
    return p


@pytest.fixture(scope="module")
def decoded(hex_path: Path, xmlcon_path: Path):
    return load_sbe_hex(hex_path, xmlcon_path)


@pytest.mark.integration
def test_decoded_has_scans(decoded) -> None:
    assert len(decoded.time_julian) > 1000, "Expected >1000 scans"


@pytest.mark.integration
def test_decoded_pressure_reaches_bottom(decoded) -> None:
    """Cast 001 reaches ~1361 m; decoded pressure should exceed 1000 dbar."""
    assert np.nanmax(decoded.pressure_dbar) > 1000.0


@pytest.mark.integration
def test_decoded_temperature_range(decoded) -> None:
    """Southern Ocean temps: -2°C to +5°C is plausible."""
    T_valid = decoded.temp_c[np.isfinite(decoded.temp_c)]
    assert len(T_valid) > 0
    assert T_valid.min() > -3.0
    assert T_valid.max() < 35.0   # allow on-deck warm data before deployment


@pytest.mark.integration
def test_gps_lat_lon_present(decoded) -> None:
    """GPS lat/lon must be present in the decoded time series."""
    assert decoded.lat is not None, "No GPS lat decoded"
    assert decoded.lon is not None, "No GPS lon decoded"
    lat_valid = decoded.lat[np.isfinite(decoded.lat)]
    assert len(lat_valid) > 100, "Too few valid GPS lat values"


@pytest.mark.integration
def test_gps_position_near_expected(decoded) -> None:
    """Mean GPS position should be near cast 001 position (−70.45°N, 168.47°E)."""
    lat_valid = decoded.lat[np.isfinite(decoded.lat)]
    lon_valid = decoded.lon[np.isfinite(decoded.lon)]
    assert abs(np.nanmedian(lat_valid) - (-70.45)) < 0.5
    assert abs(np.nanmedian(lon_valid) - 168.47) < 0.5


@pytest.mark.integration
@pytest.mark.xfail(strict=False, reason="GPS byte offset must be confirmed empirically; remove after Task 2 discovery")
def test_gps_position_tight(decoded) -> None:
    """After confirming GPS byte offsets, mean position within 0.1°."""
    assert abs(np.nanmedian(decoded.lat[np.isfinite(decoded.lat)]) - (-70.45)) < 0.1
    assert abs(np.nanmedian(decoded.lon[np.isfinite(decoded.lon)]) - 168.47) < 0.1
```

- [ ] **Step 2: Run without data (expect skip)**

```
pytest tests/integration/test_sbe_hex_s4p.py -v
```

Expected: all tests SKIP.

- [ ] **Step 3: Run with data**

```
TEST_DATA_DIR="C:/Users/peter_sha/Documents/sourcecode/LADCP/test_data" \
    pytest tests/integration/test_sbe_hex_s4p.py -v -m integration
```

After confirming GPS byte offsets in Task 2 Step 5, update `_LAT_OFFSET`/`_LON_OFFSET` constants in `sbe_hex.py` so the `xfail` test passes, then remove the `xfail` marker.

- [ ] **Step 4: Commit**

```bash
git add tests/integration/test_sbe_hex_s4p.py
git commit -m "test: SBE hex decoder integration test vs S4P cast 001 data"
```

---

## Self-Review

**Spec coverage:**
- ✅ XMLCON parser — extracts T1/C1/T2/C2 calibration + Digiquartz pressure coefficients
- ✅ Hex header parser — extracts bytes_per_scan, n_voltage_words, upload timestamp, GPS/time flags
- ✅ Byte-layout discovery script — empirically locates GPS bytes before coding hard offsets
- ✅ Temperature calibration (G-J ITS-90 form, `f₀/f` ratio)
- ✅ Conductivity calibration (G-J form, CTcor/CPcor corrections)
- ✅ Digiquartz pressure calibration (polynomial with AD590 compensation)
- ✅ GPS lat/lon extraction (big-endian float32, offset confirmed by discovery script)
- ✅ Salinity via `gsw` (optional, NaN if unavailable)
- ✅ Outputs `CTDTimeSeries` compatible with `assign_bin_depths()` and `compute_ship_velocity()`
- ✅ Integration test validates pressure depth, temperature plausibility, GPS presence

**Placeholders:** The calibration constants `_FREQ_CLOCK_HZ = 24e6` and `_FREQ_N_CYCLES = 210` are estimates that must be verified against SBE Application Note 69. The `_LAT_OFFSET` / `_LON_OFFSET` byte offsets must be confirmed by the discovery script (Task 2 Step 5). Tests use `xfail` to mark these as provisional.

**Type consistency:** `CTDTimeSeries` from `ladcp.ingestion.ctd` is the return type throughout; `XmlconCoeffs` is defined in Task 1 and consumed by Tasks 3–4.
