"""SBE 9plus hex file decoder.

Reads raw .hex telemetry + .XMLCON calibration from SBE9/11plus systems.
Produces CTDTimeSeries at the native 24 Hz sampling rate.

Reference: SeaBird Application Note 69 (SBE 9plus data format).
Sensor calibration equations: AN-04 (temperature), AN-14 (conductivity),
AN-46 (pressure).
"""
from __future__ import annotations

import re
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


# ---------------------------------------------------------------------------
# Hex header parsing
# ---------------------------------------------------------------------------

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


# Stub — implemented in Task 3
def load_sbe_hex(
    hex_path: str | Path,
    xmlcon_path: str | Path,
) -> CTDTimeSeries:
    """Load SBE hex file and return a 24 Hz CTDTimeSeries."""
    raise NotImplementedError("Implemented in Task 3")
