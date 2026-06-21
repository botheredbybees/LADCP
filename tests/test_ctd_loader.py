import numpy as np
import pytest
from pathlib import Path
from ladcp.ingestion.ctd import CTDTimeSeries, _parse_sbe_header, _map_column, load_ctd
from ladcp.ingestion._pd0 import _to_julian

BINARY_HEADER_LINES = [
    "* Sea-Bird SBE 9 Data File:\n",
    "# nquan = 4\n",
    "# name 0 = prDM: Pressure, Digiquartz [db]\n",
    "# name 1 = t090C: Temperature [ITS-90, deg C]\n",
    "# name 2 = sal00: Salinity, Practical [PSU]\n",
    "# name 3 = timeJ: Julian Days\n",
    "# bad_flag = -9.990e-29\n",
    "# file_type = binary\n",
    "# start_time = Jan 05 2015 10:30:00\n",
    "*END*\n",
]

def test_sbe_binary_header_parse():
    info = _parse_sbe_header(BINARY_HEADER_LINES)
    assert info['nquan'] == 4
    assert info['columns'][0] == 'prDM'
    assert info['columns'][1] == 't090C'
    assert info['columns'][2] == 'sal00'
    assert info['columns'][3] == 'timeJ'
    assert info['bad_flag'] == pytest.approx(-9.990e-29, rel=1e-3)
    assert info['file_type'] == 'binary'
    assert info['start_time_julian'] is not None

def test_column_name_mapping_pressure():
    assert _map_column('prDM') == 'pressure'
    assert _map_column('prSM') == 'pressure'
    assert _map_column('prE') == 'pressure'

def test_column_name_mapping_time_julian():
    assert _map_column('timeJ') == 'time_julian'

def test_column_name_mapping_temperature():
    assert _map_column('t090C') == 'temp'
    assert _map_column('t190C') == 'temp'


def _make_binary_cnv(tmp_path: Path, data: np.ndarray, bad_flag: float = -9.990e-29) -> Path:
    nquan = data.shape[1]
    header = (
        f"# nquan = {nquan}\n"
        "# name 0 = prDM: Pressure [db]\n"
        "# name 1 = t090C: Temperature\n"
        "# name 2 = sal00: Salinity\n"
        "# name 3 = timeJ: Julian Days\n"
        f"# bad_flag = {bad_flag}\n"
        "# file_type = binary\n"
        "*END*\n"
    ).encode()
    p = tmp_path / "test.cnv"
    p.write_bytes(header + data.astype('<f4').tobytes())
    return p

def test_bad_flag_masked_to_nan(tmp_path):
    BAD = -9.990e-29
    data = np.array([[100.0, BAD, 35.0, 2451545.0]], dtype=np.float64)
    p = _make_binary_cnv(tmp_path, data, bad_flag=BAD)
    result = load_ctd(p)
    assert np.isnan(result.temp_c[0])
    assert not np.isnan(result.pressure_dbar[0])

def test_format_dispatch_binary_flag(tmp_path):
    data = np.array([[100.0, 15.0, 35.0, 2451545.0],
                     [200.0, 10.0, 35.5, 2451545.1]], dtype=np.float64)
    p = _make_binary_cnv(tmp_path, data)
    result = load_ctd(p)
    assert isinstance(result, CTDTimeSeries)
    assert len(result.pressure_dbar) == 2
    assert abs(result.pressure_dbar[0] - 100.0) < 1.0
    assert abs(result.pressure_dbar[1] - 200.0) < 1.0

def test_format_dispatch_sbe_ascii(tmp_path):
    content = (
        "# nquan = 4\n"
        "# name 0 = prDM: Pressure [db]\n"
        "# name 1 = t090C: Temperature\n"
        "# name 2 = sal00: Salinity\n"
        "# name 3 = timeJ: Julian Days\n"
        "# bad_flag = -9.990e-29\n"
        "*END*\n"
        "100.000  15.000  35.000  2451545.000\n"
        "200.000  10.000  35.500  2451545.100\n"
    )
    p = tmp_path / "ascii_sbe.cnv"
    p.write_text(content)
    result = load_ctd(p)
    assert isinstance(result, CTDTimeSeries)
    assert len(result.pressure_dbar) == 2
    assert abs(result.pressure_dbar[1] - 200.0) < 1.0

def test_generic_ascii_elapsed_time(tmp_path):
    # I7N .1Hz style: elapsed seconds, fixed columns, 3-line header
    content = (
        "# station 002\n"
        "# instrument SBE9\n"
        "# fields time_s pressure_dbar temp_c salinity\n"
        "0.0  5.0  20.0  35.0\n"
        "1.0  10.0  19.5  35.1\n"
        "2.0  15.0  19.0  35.2\n"
    )
    p = tmp_path / "ctd.1hz"
    p.write_text(content)
    # Julian day 2451545.0 = 2000-01-01 00:00:00 UTC (J2000 reference)
    start_julian = 2451545.0
    result = load_ctd(
        p,
        skip_rows=3,
        col_time=0,
        col_pressure=1,
        col_temp=2,
        col_salinity=3,
        time_base="elapsed_s",
        time_start_julian=start_julian,
    )
    assert isinstance(result, CTDTimeSeries)
    assert len(result.pressure_dbar) == 3
    # time_julian = elapsed_s / 86400 + start_julian
    assert abs(result.time_julian[1] - (1.0 / 86400.0 + start_julian)) < 1e-9
    assert abs(result.pressure_dbar[2] - 15.0) < 1.0
    assert abs(result.salinity[0] - 35.0) < 0.01


from ladcp.ingestion._types import RDIData
from ladcp.ingestion.ctd import assign_bin_depths


def _make_rdi(nbin: int = 10, nens: int = 5, dist_m: float = 8.0, blen_m: float = 8.0) -> RDIData:
    zeros2d = np.zeros((nbin, nens))
    zeros1d = np.zeros(nens)
    return RDIData(
        u=zeros2d, v=zeros2d, w=zeros2d, e=zeros2d,
        heading=zeros1d, pitch=zeros1d, roll=zeros1d,
        time_julian=np.linspace(2451545.0, 2451545.5, nens),
        temp_c=zeros1d, sound_vel_ms=np.full(nens, 1500.0),
        echo=np.zeros((nbin, nens, 4), dtype=np.uint8),
        corr=np.zeros((nbin, nens, 4), dtype=np.uint8),
        pg=np.zeros((nbin, nens, 4), dtype=np.uint8),
        btrack_range_m=np.zeros((4, nens)),
        btrack_vel_ms=np.zeros((4, nens)),
        nbin=nbin, nens=nens, blen_m=blen_m, blnk_m=4.0,
        dist_m=dist_m, npng=4, coord_transform=0, serial=[0] * 8,
    )


def _make_ctd(nctd: int = 100) -> CTDTimeSeries:
    return CTDTimeSeries(
        time_julian=np.linspace(2451545.0, 2451545.5, nctd),
        pressure_dbar=np.linspace(0.0, 4000.0, nctd),
        temp_c=np.full(nctd, 15.0),
        salinity=np.full(nctd, 35.0),
    )


def test_assign_bin_depths_shape():
    rdi = _make_rdi(nbin=10, nens=5)
    ctd = _make_ctd()
    z_m, izm = assign_bin_depths(rdi, ctd)
    assert z_m.shape == (5,)
    assert izm.shape == (10, 5)


def test_assign_bin_depths_down_deeper():
    rdi = _make_rdi(nbin=5, nens=3, dist_m=8.0, blen_m=8.0)
    ctd = _make_ctd()
    z_m, izm = assign_bin_depths(rdi, ctd, looker="down")
    # Each successive bin is deeper than the instrument
    assert np.all(izm[0, :] > z_m)
    assert np.all(izm[1, :] > izm[0, :])


def test_assign_bin_depths_up_shallower():
    rdi = _make_rdi(nbin=5, nens=3, dist_m=8.0, blen_m=8.0)
    ctd = _make_ctd()
    z_m, izm = assign_bin_depths(rdi, ctd, looker="up")
    # Each successive bin is shallower than the instrument
    assert np.all(izm[0, :] < z_m)
    assert np.all(izm[1, :] < izm[0, :])


def test_assign_bin_depths_nan_propagation():
    # Use exact time-matching so np.interp returns NaN at ensemble index 1
    # (np.interp returns fp[i] exactly when x matches xp[i] — reliable NaN injection)
    nbin, nens = 4, 3
    rdi = _make_rdi(nbin=nbin, nens=nens)
    ctd_nan = CTDTimeSeries(
        time_julian=rdi.time_julian.copy(),           # exact same times
        pressure_dbar=np.array([100.0, np.nan, 300.0]),
        temp_c=np.full(nens, 15.0),
        salinity=np.full(nens, 35.0),
    )
    z_m, izm = assign_bin_depths(rdi, ctd_nan)
    assert np.isnan(z_m[1])
    assert np.all(np.isnan(izm[:, 1]))


def _make_binary_cnv_with_start_time(
    tmp_path: Path,
    data: np.ndarray,
    start_time: str,
    bad_flag: float = -9.990e-29,
) -> Path:
    nquan = data.shape[1]
    header = (
        f"# nquan = {nquan}\n"
        "# name 0 = prDM: Pressure [db]\n"
        "# name 1 = t090C: Temperature\n"
        "# name 2 = sal00: Salinity\n"
        "# name 3 = timeJ: Julian Days\n"
        f"# bad_flag = {bad_flag}\n"
        "# file_type = binary\n"
        f"# start_time = {start_time}\n"
        "*END*\n"
    ).encode()
    p = tmp_path / "test_doy.cnv"
    p.write_bytes(header + data.astype('<f4').tobytes())
    return p


def test_timej_day_of_year_converted_to_absolute_julian(tmp_path: Path):
    """SBE timeJ is day-of-year (1=Jan 1); load_ctd() must convert to absolute Julian day."""
    # Feb 5 2015 = day 36 of 2015; day 36.5 = noon on Feb 5
    data = np.array([
        [100.0, 15.0, 35.0, 5.0],   # Jan 5 00:00
        [200.0, 14.0, 35.1, 5.5],   # Jan 5 12:00
        [300.0, 13.0, 35.2, 6.0],   # Jan 6 00:00
    ], dtype=np.float64)
    p = _make_binary_cnv_with_start_time(tmp_path, data, "Jan 05 2015 10:30:00")
    result = load_ctd(p)

    jan1_2015 = _to_julian(2015, 1, 1, 0.0)
    expected_t0 = jan1_2015 + (5.0 - 1.0)   # day 5 → offset 4 days from Jan 1
    expected_t1 = jan1_2015 + (5.5 - 1.0)
    expected_t2 = jan1_2015 + (6.0 - 1.0)

    assert abs(result.time_julian[0] - expected_t0) < 1e-9
    assert abs(result.time_julian[1] - expected_t1) < 1e-9
    assert abs(result.time_julian[2] - expected_t2) < 1e-9
    # Values must be absolute Julian days (~2.4 million), not day-of-year (1–366)
    assert result.time_julian[0] > 1000


def _make_binary_cnv_latlon(
    tmp_path: Path,
    lat_vals: list[float],
    lon_vals: list[float],
    bad_lon_idx: int | None = None,
) -> Path:
    """Binary CNV with 'Store Lat/Lon Data' flag and nquan=6."""
    BAD = -9.990e-29
    rows = []
    for i, (la, lo) in enumerate(zip(lat_vals, lon_vals)):
        lo_val = BAD if (bad_lon_idx is not None and i == bad_lon_idx) else lo
        rows.append([100.0 + i * 100, 15.0 - i, 35.0 + i * 0.1,
                     2451545.0 + i * 0.1, la, lo_val])
    data = np.array(rows, dtype=np.float64)
    nquan = 6
    header = (
        f"# nquan = {nquan}\n"
        "# name 0 = prDM: Pressure [db]\n"
        "# name 1 = t090C: Temperature\n"
        "# name 2 = sal00: Salinity\n"
        "# name 3 = timeJ: Julian Days\n"
        "# name 4 = latitude: Latitude\n"
        "# name 5 = longitude: Longitude\n"
        f"# bad_flag = {BAD}\n"
        "# file_type = binary\n"
        "* Store Lat/Lon Data = Append to Every Scan\n"
        "*END*\n"
    ).encode()
    p = tmp_path / "test_latlon.cnv"
    p.write_bytes(header + data.astype("<f4").tobytes())
    return p


def test_ctd_loads_lat_lon(tmp_path: Path):
    """CNV with 'Store Lat/Lon Data' flag populates lat/lon of correct length."""
    p = _make_binary_cnv_latlon(tmp_path, [-30.0, -30.01], [-140.0, -139.99])
    result = load_ctd(p)
    assert result.lat is not None
    assert result.lon is not None
    assert len(result.lat) == 2
    assert len(result.lon) == 2
    assert abs(result.lat[0] - (-30.0)) < 0.01
    assert abs(result.lon[0] - (-140.0)) < 0.01
    # Standard columns must still be populated correctly
    assert abs(result.pressure_dbar[0] - 100.0) < 1.0
    assert abs(result.pressure_dbar[1] - 200.0) < 1.0


def test_ctd_no_lat_lon_returns_none(tmp_path: Path):
    """CNV without 'Store Lat/Lon Data' flag returns lat=None, lon=None."""
    data = np.array([[100.0, 15.0, 35.0, 2451545.0]], dtype=np.float64)
    p = _make_binary_cnv(tmp_path, data)
    result = load_ctd(p)
    assert result.lat is None
    assert result.lon is None


def test_ctd_bad_flag_lat_lon_becomes_nan(tmp_path: Path):
    """SBE bad-flag sentinel in a lat/lon column becomes NaN."""
    p = _make_binary_cnv_latlon(
        tmp_path,
        [-30.0, -30.01],
        [-140.0, -139.99],
        bad_lon_idx=0,  # row 0 lon is bad
    )
    result = load_ctd(p)
    assert result.lon is not None
    assert np.isnan(result.lon[0])
    assert not np.isnan(result.lon[1])
    assert not np.isnan(result.lat[0])
