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


# ---------------------------------------------------------------------------
# Calibration equations (AN-04 temperature, AN-14 conductivity, AN-46 pressure)
# ---------------------------------------------------------------------------

# Clock and gate constants for SBE9plus frequency channels (AN-69).
# SBE9plus: f_clock = 24 MHz, N_cycles = 210 for T/C channels.
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
    # Convert psia to dbar (1 psia ≈ 0.68947572932 dbar, subtract 1 atm = 14.6959 psia)
    return max((P_psia - 14.6959) * 0.68947572932, 0.0)


# ---------------------------------------------------------------------------
# GPS byte offsets and load_sbe_hex()
# ---------------------------------------------------------------------------

# GPS encoding: SBE11plus V5.x stores lat/lon as signed int32 scaled by 1e-7
# (confirmed via Task 2 discovery: float32 interpretation produced garbage).
# Validated at runtime — values outside [-90,90] / [-180,180] become NaN.
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
    cond = np.full(n_scans, np.nan, dtype=np.float64)
    lat = np.full(n_scans, np.nan, dtype=np.float64)
    lon = np.full(n_scans, np.nan, dtype=np.float64)

    for i, hex_line in enumerate(data_lines):
        time_jd[i] = t0_jd + i * dt_per_scan

        if len(hex_line) < n_bytes * 2:
            continue

        scan = bytes.fromhex(hex_line[: n_bytes * 2])

        # Decode T1, C1, and P (P needed for C calibration)
        count_t1 = int.from_bytes(scan[0:3], "big")
        count_c1 = int.from_bytes(scan[3:6], "big")
        count_p = scan[12:15]

        T1 = _temperature(count_t1, coeffs, primary=True)
        P = _pressure(count_p, coeffs)
        C1 = _conductivity(count_c1, T1, P, coeffs, primary=True)

        temp[i] = T1
        pressure[i] = P
        cond[i] = C1

        # GPS lat/lon: try int32 scaled by 1e-7, validate range
        if hdr.nmea_pos_added and n_bytes >= _LON_OFFSET + 4:
            try:
                lat_raw = struct.unpack(">i", scan[_LAT_OFFSET: _LAT_OFFSET + 4])[0]
                lon_raw = struct.unpack(">i", scan[_LON_OFFSET: _LON_OFFSET + 4])[0]
                lat_val = lat_raw * 1e-7
                lon_val = lon_raw * 1e-7
                if -90.0 <= lat_val <= 90.0:
                    lat[i] = lat_val
                if -180.0 <= lon_val <= 180.0:
                    lon[i] = lon_val
            except struct.error:
                pass

    # Salinity (requires gsw; set to NaN array if unavailable)
    salinity = np.full(n_scans, np.nan, dtype=np.float64)
    import importlib.util
    if importlib.util.find_spec("gsw") is not None:
        import gsw
        try:
            salinity = gsw.SP_from_C(cond, temp, pressure)
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
