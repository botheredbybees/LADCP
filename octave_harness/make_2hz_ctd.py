"""Generate LDEO_IX's expected `003.2Hz` ASCII CTD/nav file from the raw
P16N cast003 binary .cnv, per CONTINUATION_PLAN.md's "gift" section.

We do NOT have the original 003.2Hz LDEO used -- the recorded p-struct
(octave_harness/recorded_p_struct_attrs.txt, pulled from test_data/2015_P16N/
003.nc global attributes) tells us its exact field layout:

    ctd_fields_per_line: 11, ctd_header_lines: 0
    field 1 = time, field 2 = pressure, field 3 = temperature, field 4 = salinity
    field 10 = lat, field 11 = lon
    bad value flag: -999

    ctd_time_base/nav_time_base=0 means "elapsed time in SECONDS since
    p.time_start" (see loadctd.m/loadnav.m: `case 0: timctd =
    timctd/24/3600 + julian(p.time_start)`) -- NOT absolute Julian day as
    an earlier reading of this comment assumed. Confirmed by the recorded
    p.dt_profile = 11156.999994814396 (seconds, ~3.1 h -- matches this
    cast's real duration) and p.time_start = [2015 4 11 18 0 44]. Field 1
    below is therefore elapsed seconds since this file's own first scan,
    with the matching p.time_start printed for set_cast_params.m.

`003_01.cnv` is Sea-Bird binary format, ~24 Hz (interval = 0.0416667 s),
columns [pressure, t090C, cond mS/cm, oxygen ml/l, timeJ, oxygen mg/l, flag]
(confirmed against the file's own `# span` header lines). loadctd.m/loadnav.m
only read fields 1-4 and 10-11 (see their `switch` statements) -- fields
5-9 are `%*g`-skipped, so any filler value is fine there.

Caveats (documented per CONTINUATION_PLAN.md M2 instructions):
  - No per-scan lat/lon in the cnv (nquan=7, no lat/lon column) -- we use
    the recorded constant p.lat/p.lon for every scan (plan explicitly
    allows this for a first pass).
  - Our own Python pipeline's CTD loader (`ladcp.ingestion.ctd.load_ctd`)
    currently returns salinity=NaN for this file (column-name matching
    finds no `^sal` column) -- so there is no existing Python-computed
    salinity to mirror. We compute practical salinity via the PSS-78
    algorithm (Fofonoff & Millard 1983 / Seabird AN#14), reimplemented
    here from docs/legacy/pycurrents/src/pycurrents/data/seawater.py
    (read there, not imported or modified) since the Octave pipeline's
    downstream sound-speed correction needs a physically valid S.
  - Temperature is converted ITS-90 -> IPTS-68 (T68 = T90 * 1.00024)
    because the PSS-78 formula (and this-era LDEO code) assumes IPTS-68;
    the same T68 value is written to field 3.
"""
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
CNV_PATH = REPO / "test_data" / "2015_P16N" / "003_01.cnv"
OUT_PATH = REPO / "octave_harness" / "work" / "data" / "CTD" / "2Hz" / "003.2Hz"

# Recorded p-struct values (test_data/2015_P16N/003.nc global attrs).
RECORDED_LAT = -15.498335
RECORDED_LON = -150.19699
RECORDED_CTD_STARTTIME = 2457124.7505048458
RECORDED_CTD_ENDTIME = 2457124.879635498

# Julian-day base for this custom (midnight-epoch) julian.m convention,
# for calendar year 2015: full_julian = JULIAN_BASE_2015 + timeJ, where
# timeJ is Sea-Bird's "Julian Days" (day-of-year, Jan 1 00:00 = 1.0).
# Derived from octave_harness/ldeo_ix/julian.m: julian(2015,1,1,0) = 2457024.0,
# so day-of-year 1.0 maps to 2457024.0 => base = 2457024.0 - 1.0 = 2457023.0.
JULIAN_BASE_2015 = 2457023.0

# julian.m's custom (midnight-epoch) reference point, for converting a
# full_julian value back to a calendar datetime (needed for p.time_start).
JULIAN_REF_DATE = datetime(2015, 4, 11)
JULIAN_REF_VALUE = 2457124.0  # julian(2015,4,11,0) per julian.m


def _julian_to_datetime(full_julian: float) -> datetime:
    return JULIAN_REF_DATE + timedelta(days=full_julian - JULIAN_REF_VALUE)

DECIMATE = 12  # 24 Hz (interval 0.0416667 s) -> 2 Hz


def _parse_cnv(path: Path) -> np.ndarray:
    data = path.read_bytes()
    marker = b"*END*"
    idx = data.find(marker)
    if idx < 0:
        raise ValueError(f"*END* marker not found in {path}")
    j = idx + len(marker)
    while data[j : j + 1] in (b"\r", b"\n"):
        j += 1
    header = data[:idx].decode("ascii", errors="replace")
    nquan = int(next(l for l in header.splitlines() if "nquan" in l).split("=")[1])
    nvalues = int(next(l for l in header.splitlines() if "nvalues" in l).split("=")[1])
    remaining = len(data) - j
    expected = nquan * nvalues * 4
    if remaining < expected:
        raise ValueError(f"binary data section too short: have {remaining}, need {expected}")
    arr = np.frombuffer(data[j : j + expected], dtype="<f4").reshape(nvalues, nquan)
    return arr.astype(np.float64)


def _pss78_salinity(cond_sm: np.ndarray, temp68: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """PSS-78 practical salinity. C in S/m, T in IPTS-68 deg C, P in dbar.

    Reimplementation of the algorithm in
    docs/legacy/pycurrents/src/pycurrents/data/seawater.py:salinity()
    (IEEE J. Oceanic Eng. OE-5 No.1 Jan 1980 p.14 / Seabird AN#14).
    """
    A1, A2, A3 = 2.070e-5, -6.370e-10, 3.989e-15
    B1, B2, B3, B4 = 3.426e-2, 4.464e-4, 4.215e-1, -3.107e-3
    c0, c1, c2, c3, c4 = 6.766097e-1, 2.00564e-2, 1.104259e-4, -6.9698e-7, 1.0031e-9
    a0, a1, a2, a3, a4, a5 = 0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081
    b0, b1, b2, b3, b4, b5 = 0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144
    k = 0.0162

    R = cond_sm / 4.2914
    Rp_num = np.polyval([A3, A2, A1], pressure) * pressure
    Rp_den = np.polyval([B2, B1, 1], temp68) + B3 * R + B4 * R * temp68
    Rp = 1 + Rp_num / Rp_den
    rT = np.polyval([c4, c3, c2, c1, c0], temp68)
    RT = R / (Rp * rT)
    sq_RT = np.sqrt(RT)
    Sa = np.polyval([a5, a4, a3, a2, a1, a0], sq_RT)
    Sb = np.polyval([b5, b4, b3, b2, b1, b0], sq_RT)
    Tfac = (temp68 - 15) / (1 + k * (temp68 - 15))
    return Sa + Tfac * Sb


def main() -> None:
    arr = _parse_cnv(CNV_PATH)
    n = arr.shape[0]
    n_groups = n // DECIMATE
    trimmed = arr[: n_groups * DECIMATE].reshape(n_groups, DECIMATE, arr.shape[1])
    avg = trimmed.mean(axis=1)

    pressure = avg[:, 0]
    temp90 = avg[:, 1]
    cond_mScm = avg[:, 2]
    timeJ = avg[:, 4]

    temp68 = temp90 * 1.00024
    cond_sm = cond_mScm / 10.0
    salinity = _pss78_salinity(cond_sm, temp68, pressure)
    full_julian = JULIAN_BASE_2015 + timeJ
    elapsed_sec = (full_julian - full_julian[0]) * 86400.0

    time_start_dt = _julian_to_datetime(full_julian[0])
    time_end_dt = _julian_to_datetime(full_julian[-1])

    print(f"decimated {n} scans -> {n_groups} groups (dropped {n - n_groups * DECIMATE})")
    print(f"generated time span (full julian): {full_julian[0]!r} .. {full_julian[-1]!r}")
    print(f"recorded ctd_starttime/endtime:     {RECORDED_CTD_STARTTIME!r} .. {RECORDED_CTD_ENDTIME!r}")
    print(f"time_start (calendar): {time_start_dt.isoformat()}")
    print(f"time_end   (calendar): {time_end_dt.isoformat()}")
    print(
        "p.time_start = [{} {} {} {} {} {}];".format(
            time_start_dt.year, time_start_dt.month, time_start_dt.day,
            time_start_dt.hour, time_start_dt.minute,
            time_start_dt.second + time_start_dt.microsecond / 1e6,
        )
    )
    print(
        "p.time_end   = [{} {} {} {} {} {}];".format(
            time_end_dt.year, time_end_dt.month, time_end_dt.day,
            time_end_dt.hour, time_end_dt.minute,
            time_end_dt.second + time_end_dt.microsecond / 1e6,
        )
    )
    print(f"elapsed_sec range: {elapsed_sec[0]:.3f} .. {elapsed_sec[-1]:.3f} s")
    print(f"pressure range: {pressure.min():.1f} .. {pressure.max():.1f} dbar")
    print(f"salinity range: {salinity.min():.3f} .. {salinity.max():.3f} PSU")

    filler = np.zeros(n_groups)
    lat = np.full(n_groups, RECORDED_LAT)
    lon = np.full(n_groups, RECORDED_LON)

    cols = [elapsed_sec, pressure, temp68, salinity, filler, filler, filler, filler, filler, lat, lon]
    out = np.column_stack(cols)
    out = np.where(np.isfinite(out), out, -999.0)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(OUT_PATH, out, fmt="%.8f")
    print(f"wrote {OUT_PATH} ({n_groups} lines, {out.shape[1]} fields)")


if __name__ == "__main__":
    main()
