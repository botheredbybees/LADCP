"""Low-level Teledyne RDI PD0 binary parser.

Reference: docs/legacy/loadrdi.m::rdread.
"""

import struct

import numpy as np

_VEL_BAD = -32768
_VEL_SCALE = 0.001  # int16 LSB → m/s
_BT_RANGE_SCALE = 0.01  # uint16 LSB → m
_BT_VEL_SCALE = 0.001  # int16 LSB → m/s
_FL_LEN_SCALE = 0.01  # uint16 cm → m

# Data type IDs from loadrdi.m varid array
_ID_FIXED = 0x0000
_ID_VARIABLE = 0x0080
_ID_VELOCITY = 0x0100
_ID_CORR = 0x0200
_ID_ECHO = 0x0300
_ID_PG = 0x0400
_ID_BTRACK = 0x0600


def parse_pd0(data: bytes) -> list[dict]:
    """Parse a PD0 binary blob; return one dict per valid ensemble.

    Follows the structure of rdread()/rdhead()/rdflead()/rdvlead()/rdbtrack()
    in docs/legacy/loadrdi.m.  Bad velocity values (-32768) are replaced with NaN.
    Ensembles with a failed checksum are silently dropped.
    """
    ensembles = []
    offset = 0
    n = len(data)

    while offset < n - 5:
        # Find 0x7F 0x7F sync
        if data[offset] != 0x7F or data[offset + 1] != 0x7F:
            offset += 1
            continue

        if offset + 4 > n:
            break
        nbytes = struct.unpack_from("<H", data, offset + 2)[0]
        end = offset + nbytes + 2  # +2 for checksum
        if end > n:
            offset += 1
            continue

        body = data[offset : offset + nbytes]
        checksum = struct.unpack_from("<H", data, offset + nbytes)[0]
        if sum(body) % 65536 != checksum:
            offset += nbytes + 2
            continue

        ens = _parse_ensemble(body)
        if ens is not None:
            ensembles.append(ens)
        offset += nbytes + 2

    return ensembles


def _parse_ensemble(body: bytes) -> dict | None:
    """Parse a single ensemble body (without checksum)."""
    if len(body) < 6:
        return None

    # Header: magic(2) + nbytes(2) + spare(1) + ndt(1) + offsets(2*ndt)
    ndt = body[5]
    if ndt > 16 or 6 + 2 * ndt > len(body):
        return None
    offsets = list(struct.unpack_from(f"<{ndt}H", body, 6))

    # Read data type IDs at each offset
    blocks: dict[int, int] = {}  # type_id → offset
    for off in offsets:
        if off + 2 > len(body):
            continue
        type_id = struct.unpack_from("<H", body, off)[0]
        blocks[type_id] = off

    required = (_ID_FIXED, _ID_VARIABLE, _ID_VELOCITY)
    if any(k not in blocks for k in required):
        return None

    fl = _read_fixed_leader(body, blocks[_ID_FIXED] + 2)
    nbin = fl["nbin"]
    vl = _read_variable_leader(body, blocks[_ID_VARIABLE] + 2)
    vel = _read_matrix(
        body, blocks[_ID_VELOCITY] + 2, nbin, "<i2", _VEL_SCALE, _VEL_BAD
    )
    corr = (
        _read_matrix_uint8(body, blocks[_ID_CORR] + 2, nbin)
        if _ID_CORR in blocks
        else None
    )
    echo = (
        _read_matrix_uint8(body, blocks[_ID_ECHO] + 2, nbin)
        if _ID_ECHO in blocks
        else None
    )
    pg = (
        _read_matrix_uint8(body, blocks[_ID_PG] + 2, nbin) if _ID_PG in blocks else None
    )
    bt = (
        _read_bottom_track(body, blocks[_ID_BTRACK] + 2)
        if _ID_BTRACK in blocks
        else None
    )

    return {
        "fixed_leader": fl,
        "variable_leader": vl,
        "velocity": vel,
        "correlation": corr,
        "echo": echo,
        "percent_good": pg,
        "bottom_track": bt,
    }


def _read_fixed_leader(body: bytes, start: int) -> dict:
    """Reference: rdflead() in loadrdi.m.

    Byte layout (offsets relative to start, i.e. after the 2-byte ID):
      skip 7 bytes (CPU firmware/feature flags)
      nbin        uint8   @ +7
      npng        uint16  @ +8
      blen_cm     uint16  @ +10
      blnk_cm     uint16  @ +12
      skip 16 bytes       @ +14
      coord_transform uint8 @ +23
      dist_cm     uint16  @ +30
      plen_cm     uint16  @ +32
      skip 6 bytes        @ +34
      serial      8×uint8 @ +40
    """
    p = start + 7
    nbin = body[p]
    npng, blen_cm, blnk_cm = struct.unpack_from("<HHH", body, p + 1)
    # skip 16 bytes (water profiling mode, correlation threshold, etc.)
    p2 = p + 1 + 6 + 16
    coord_transform = body[start + 23]
    dist_cm, plen_cm = struct.unpack_from("<HH", body, p2)
    # skip 6 bytes (ref layer, false target, spare, bandwidth)
    p3 = p2 + 4 + 6
    serial = list(body[p3 : p3 + 8])
    return {
        "nbin": nbin,
        "npng": npng,
        "blen_m": blen_cm * _FL_LEN_SCALE,
        "blnk_m": blnk_cm * _FL_LEN_SCALE,
        "coord_transform": coord_transform,
        "dist_m": dist_cm * _FL_LEN_SCALE,
        "plen_m": plen_cm * _FL_LEN_SCALE,
        "serial": serial,
    }


def _read_variable_leader(body: bytes, start: int) -> dict:
    """Reference: rdvlead() in loadrdi.m.

    Byte layout (offsets relative to start, i.e. after the 2-byte ID):
      skip 2 bytes (ensemble number low+high) @ +0
      7-byte time (yy,mm,dd,hh,mm,ss,cc)      @ +2
      skip 3 bytes (RTC century + ens MSB)    @ +9
      sound_vel   uint16                       @ +12
      skip 2 bytes (depth of transducer)       @ +14
      heading     uint16                       @ +16
      pitch       int16                        @ +18
      roll        int16                        @ +20
      salinity    uint16                       @ +22
      temp        int16                        @ +24
      skip 6 bytes (MPT min/sec/hun + stdevs)  @ +26
      xmt_current uint8                        @ +32
      xmt_volt    uint8                        @ +33
      int_temp    uint8                        @ +34
    """
    # skip 2 bytes (ensemble number low + high)
    p = start + 2
    yy, mo, dd, hh, mm, ss, cc = struct.unpack_from("BBBBBBB", body, p)
    year = 2000 + yy if yy < 80 else 1900 + yy
    time_julian = _to_julian(year, mo, dd, hh + mm / 60 + ss / 3600 + cc / 360000)
    # skip 3 bytes (real-time clock century + ensemble MSB, bit)
    p2 = p + 7 + 3
    sound_vel = struct.unpack_from("<H", body, p2)[0]
    # skip 2 (depth of transducer)
    p3 = p2 + 2 + 2
    heading_01 = struct.unpack_from("<H", body, p3)[0]
    pitch_01, roll_01 = struct.unpack_from("<hh", body, p3 + 2)
    salinity_ppt = struct.unpack_from("<H", body, p3 + 6)[0] * 0.001
    temp_01 = struct.unpack_from("<h", body, p3 + 8)[0]
    # skip 6 bytes (MPT minutes/seconds/hundredths, heading/pitch/roll std)
    p4 = p3 + 10 + 6
    xmt_cur, xmt_volt, int_temp = struct.unpack_from("BBB", body, p4)
    return {
        "time_julian": time_julian,
        "heading_deg": heading_01 * 0.01,
        "pitch_deg": pitch_01 * 0.01,
        "roll_deg": roll_01 * 0.01,
        "temp_c": temp_01 * 0.01,
        "salinity_psu": salinity_ppt,
        "sound_vel_ms": float(sound_vel),
        "xmt_current": xmt_cur,
        "xmt_volt": xmt_volt,
        "int_temp": int_temp,
    }


def _read_matrix(
    body: bytes, start: int, nbin: int, dtype: str, scale: float, bad: int
) -> np.ndarray:
    """Read (nbin × 4) matrix of int16 values, apply scaling, replace bad→NaN."""
    nbytes = nbin * 4 * 2
    if start + nbytes > len(body):
        return np.full((nbin, 4), np.nan)
    raw = (
        np.frombuffer(body[start : start + nbytes], dtype=dtype)
        .reshape(nbin, 4)
        .astype(np.float64)
    )
    raw[raw == bad] = np.nan
    return raw * scale


def _read_matrix_uint8(body: bytes, start: int, nbin: int) -> np.ndarray:
    """Read (nbin × 4) matrix of uint8 values."""
    nbytes = nbin * 4
    if start + nbytes > len(body):
        return np.zeros((nbin, 4), dtype=np.uint8)
    arr = np.frombuffer(body[start : start + nbytes], dtype=np.uint8)
    return arr.reshape(nbin, 4).copy()


def _read_bottom_track(body: bytes, start: int) -> dict:
    """Reference: rdbtrack() in loadrdi.m. Skip 14 bytes, then range + velocity."""
    p = start + 14
    if p + 16 > len(body):
        return {
            "range_m": np.full(4, np.nan),
            "vel_ms": np.full(4, np.nan),
            "corr": np.zeros(4, np.uint8),
            "pg": np.zeros(4, np.uint8),
        }
    range_raw = np.array(struct.unpack_from("<4H", body, p), dtype=np.float64)
    vel_raw = np.array(struct.unpack_from("<4h", body, p + 8), dtype=np.float64)
    corr_pg = np.frombuffer(body[p + 16 : p + 24], dtype=np.uint8)
    range_m = range_raw * _BT_RANGE_SCALE
    range_m[range_m == 0] = np.nan
    vel_ms = vel_raw * _BT_VEL_SCALE
    vel_ms[vel_raw == _VEL_BAD] = np.nan
    return {
        "range_m": range_m,
        "vel_ms": vel_ms,
        "corr": corr_pg[:4].copy(),
        "pg": corr_pg[4:8].copy() if len(corr_pg) >= 8 else np.zeros(4, np.uint8),
    }


def _to_julian(year: int, month: int, day: int, hour_frac: float) -> float:
    """Midnight-based Julian day matching docs/legacy/julian.m (Fliegel/Van Flandern).

    In this convention JD 2440000 began at 0000 hours, May 23, 1968.
    Unlike the formal astronomical definition (noon-to-noon), this counts
    from midnight — matching the MATLAB reference and the LDEO processing chain.
    """
    mo = month + 9
    yr = year - 1
    if month > 2:
        mo = month - 3
        yr = year
    c = yr // 100
    yr = yr - c * 100
    j = (146097 * c) // 4 + (1461 * yr) // 4 + (153 * mo + 2) // 5 + day + 1721119
    return float(j) + hour_frac / 24.0
