"""Unit tests for the low-level PD0 binary parser."""

import struct

import numpy as np


def _make_fixed_leader(
    nbin=10, npng=1, blen_cm=800, blnk_cm=176, dist_cm=1024, serial=None
):
    """Build a minimal fixed leader block (without the 2-byte ID prefix).

    Byte layout matches rdflead() in docs/legacy/loadrdi.m:
      skip 7 bytes (CPU firmware/feature flags)
      nbin        uint8   @ offset  7
      npng        uint16  @ offset  8
      blen_cm     uint16  @ offset 10
      blnk_cm     uint16  @ offset 12
      skip 16 bytes       @ offset 14-29
      dist_cm     uint16  @ offset 30
      plen_cm     uint16  @ offset 32
      skip 6 bytes        @ offset 34-39
      serial      8×uint8 @ offset 40
    """
    if serial is None:
        serial = [0] * 8
    buf = bytearray(60)
    buf[7] = nbin
    struct.pack_into("<HHH", buf, 8, npng, blen_cm, blnk_cm)
    struct.pack_into("<HH", buf, 30, dist_cm, 800)  # offset 30, not 32
    for i, b in enumerate(serial[:8]):
        buf[40 + i] = b  # offset 40, not 42
    return bytes(buf)


def _make_variable_leader(
    year=2018,
    month=11,
    day=5,
    hour=12,
    minute=0,
    second=0,
    hundredths=0,
    pitch_01deg=100,
    roll_01deg=200,
    heading_01deg=18000,
    temp_01c=200,
    salinity_ppt=35,
    sound_vel_ms=1500,
):
    """Build a minimal variable leader block (without the 2-byte ID prefix).

    Byte layout matches rdvlead() in docs/legacy/loadrdi.m:
      skip 2 bytes (ensemble number low+high)    @ offset  0-1
      7-byte time (yy,mm,dd,hh,mm,ss,cc uint8)  @ offset  2-8
      skip 3 bytes (RTC century + ens MSB + bit) @ offset  9-11
      sound_vel   uint16                         @ offset 12
      skip 2 bytes (depth of transducer)         @ offset 14-15
      heading     uint16                         @ offset 16
      pitch       int16                          @ offset 18
      roll        int16                          @ offset 20
      salinity    uint16                         @ offset 22
      temp        int16                          @ offset 24
      skip 6 bytes (MPT min/sec/hun + std devs)  @ offset 26-31
      xmt_current uint8                          @ offset 32
      xmt_volt    uint8                          @ offset 33
      int_temp    uint8                          @ offset 34
    """
    buf = bytearray(65)
    year_2d = year % 100
    struct.pack_into(
        "BBBBBBB", buf, 2, year_2d, month, day, hour, minute, second, hundredths
    )
    struct.pack_into("<H", buf, 12, sound_vel_ms)  # offset 12, not 14
    struct.pack_into("<H", buf, 16, heading_01deg)  # offset 16, not 18
    struct.pack_into("<hh", buf, 18, pitch_01deg, roll_01deg)  # offset 18, not 20
    struct.pack_into("<H", buf, 22, salinity_ppt * 1000)  # offset 22, not 24
    struct.pack_into("<h", buf, 24, temp_01c)  # offset 24, not 26
    buf[32] = 50  # xmt_current  (offset 32, not 33)
    buf[33] = 48  # xmt_volt     (offset 33, not 34)
    buf[34] = 22  # int_temp     (offset 34, not 35)
    return bytes(buf)


def _make_velocity_block(nbin=3, values=None):
    """Build velocity data block (without 2-byte ID). int16 LE, 4 beams, nbin bins."""
    if values is None:
        values = np.zeros((nbin, 4), dtype=np.int16)
    # Layout: (nbin * 4) int16 values, row-major, 4 beams per bin
    buf = values.astype("<i2").tobytes()
    return buf


def _make_minimal_ensemble(nbin=3):
    """Build a complete minimal PD0 ensemble with known values."""
    # IDs as used by rdread: fixed=0x0000, variable=0x0080, velocity=0x0100
    # correlation=0x0200, echo=0x0300, percent_good=0x0400, bottom_track=0x0600
    id_fixed = struct.pack("<H", 0x0000)
    id_var = struct.pack("<H", 0x0080)
    id_vel = struct.pack("<H", 0x0100)
    id_corr = struct.pack("<H", 0x0200)
    id_echo = struct.pack("<H", 0x0300)
    id_pg = struct.pack("<H", 0x0400)

    fl_body = _make_fixed_leader(nbin=nbin)
    vl_body = _make_variable_leader()
    vel_vals = np.array(
        [[100, -200, 50, 10], [150, -250, 60, 15], [0, 0, 0, -32768]], dtype=np.int16
    )  # last row bad
    vel_body = _make_velocity_block(nbin=nbin, values=vel_vals)
    corr_body = bytes([80] * (nbin * 4))
    echo_body = bytes([90] * (nbin * 4))
    pg_body = bytes([100] * (nbin * 4))

    # Header: 0x7F 0x7F, then build offset table
    blocks = [
        id_fixed + fl_body,
        id_var + vl_body,
        id_vel + vel_body,
        id_corr + corr_body,
        id_echo + echo_body,
        id_pg + pg_body,
    ]
    n_types = len(blocks)
    header_size = 6 + 2 * n_types  # magic(2) + nbytes(2) + spare(1) + ndt(1) + offsets

    offsets = []
    pos = header_size
    for b in blocks:
        offsets.append(pos)
        pos += len(b)

    nbytes = pos  # total bytes before checksum
    header = struct.pack(
        "<BBHBb" + "H" * n_types, 0x7F, 0x7F, nbytes, 0, n_types, *offsets
    )
    body = header + b"".join(blocks)
    checksum = struct.pack("<H", sum(body) % 65536)
    return body + checksum


class TestPd0Parser:
    def test_parse_returns_list(self):
        from ladcp.ingestion._pd0 import parse_pd0

        data = _make_minimal_ensemble(nbin=3)
        ensembles = parse_pd0(data)
        assert isinstance(ensembles, list)
        assert len(ensembles) == 1

    def test_fixed_leader_fields(self):
        from ladcp.ingestion._pd0 import parse_pd0

        data = _make_minimal_ensemble(nbin=3)
        fl = parse_pd0(data)[0]["fixed_leader"]
        assert fl["nbin"] == 3
        assert fl["npng"] == 1
        assert abs(fl["blen_m"] - 8.0) < 0.001  # 800 cm → 8.0 m
        assert abs(fl["blnk_m"] - 1.76) < 0.001  # 176 cm → 1.76 m
        assert abs(fl["dist_m"] - 10.24) < 0.001  # 1024 cm → 10.24 m

    def test_variable_leader_heading(self):
        from ladcp.ingestion._pd0 import parse_pd0

        data = _make_minimal_ensemble(nbin=3)
        vl = parse_pd0(data)[0]["variable_leader"]
        assert abs(vl["heading_deg"] - 180.0) < 0.01  # 18000 * 0.01
        assert abs(vl["pitch_deg"] - 1.0) < 0.01  # 100 * 0.01
        assert abs(vl["roll_deg"] - 2.0) < 0.01  # 200 * 0.01
        assert abs(vl["temp_c"] - 2.0) < 0.01  # 200 * 0.01
        assert abs(vl["sound_vel_ms"] - 1500) < 0.1

    def test_velocity_scaling_and_bad_values(self):
        from ladcp.ingestion._pd0 import parse_pd0

        data = _make_minimal_ensemble(nbin=3)
        vel = parse_pd0(data)[0]["velocity"]  # shape (nbin, 4)
        assert vel.shape == (3, 4)
        assert abs(vel[0, 0] - 0.100) < 1e-4  # 100 * 0.001
        assert abs(vel[0, 1] - (-0.200)) < 1e-4
        assert np.isnan(vel[2, 3])  # -32768 → NaN

    def test_checksum_mismatch_drops_ensemble(self):
        from ladcp.ingestion._pd0 import parse_pd0

        data = bytearray(_make_minimal_ensemble(nbin=3))
        data[-1] ^= 0xFF  # corrupt checksum
        ensembles = parse_pd0(bytes(data))
        assert len(ensembles) == 0

    def test_two_ensembles(self):
        from ladcp.ingestion._pd0 import parse_pd0

        ens = _make_minimal_ensemble(nbin=3)
        ensembles = parse_pd0(ens + ens)
        assert len(ensembles) == 2
