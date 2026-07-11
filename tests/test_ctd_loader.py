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


from ladcp.ingestion.ctd import compute_ship_velocity


def test_compute_ship_velocity_linear_track():
    """Straight eastward track at 1.0 m/s returns (u≈1.0, v≈0.0)."""
    import math
    # At lat=0 deg, 1 deg longitude = 111320 m. Moving east at 1 m/s:
    # Δlon_deg/s = 1.0 / 111320
    lat0 = 0.0
    speed_mps = 1.0
    dt_s = np.array([0.0, 100.0, 200.0, 300.0])
    east_m = speed_mps * dt_s
    lon = lon0 = -140.0
    lon_arr = lon0 + east_m / (math.cos(math.radians(lat0)) * 111320.0)
    lat_arr = np.full(4, lat0)
    t0 = 2457100.0  # arbitrary Julian day
    time_jul = t0 + dt_s / 86400.0
    u_ship, v_ship = compute_ship_velocity(lat_arr, lon_arr, time_jul)
    assert abs(u_ship - 1.0) < 0.01, f"u_ship={u_ship:.4f} expected ≈1.0"
    assert abs(v_ship) < 0.01, f"v_ship={v_ship:.4f} expected ≈0.0"


def test_compute_ship_velocity_insufficient_data():
    """All-NaN lat/lon returns (0.0, 0.0)."""
    lat = np.full(5, np.nan)
    lon = np.full(5, np.nan)
    t0 = 2457100.0
    time_jul = t0 + np.arange(5) / 86400.0
    u_ship, v_ship = compute_ship_velocity(lat, lon, time_jul)
    assert u_ship == 0.0
    assert v_ship == 0.0


# --- pressure_to_depth (Saunders & Fofonoff p2z, loadctd.m parity) ---

from ladcp.ingestion.ctd import estimate_ctd_adcp_lag, pressure_to_depth


def test_pressure_to_depth_saunders_check_value():
    # Documented check value in loadctd.m::p2z: depth = 9712.654 m
    # for p = 1000 bars (= 10000 dbar), lat = 30 deg.
    z = pressure_to_depth(np.array([10000.0]), lat_deg=30.0)
    assert abs(z[0] - 9712.654) < 0.01


def test_pressure_to_depth_shallower_than_fallback_at_depth():
    # At 4400 dbar the Saunders depth is ~89 m shallower than p*1.00445
    # (the izm registration offset root cause, octave_harness REPORT.md P2).
    p = np.array([4400.0])
    z = pressure_to_depth(p, lat_deg=-15.0)
    assert 80.0 < (p[0] * 1.00445 - z[0]) < 100.0


def test_assign_bin_depths_lat_uses_saunders():
    rdi = _make_rdi(nbin=4, nens=3)
    ctd = _make_ctd()
    z_m, _ = assign_bin_depths(rdi, ctd, looker="down", lat_deg=-15.0)
    p_interp = np.interp(rdi.time_julian, ctd.time_julian, ctd.pressure_dbar)
    np.testing.assert_allclose(z_m, pressure_to_depth(p_interp, lat_deg=-15.0))


def test_assign_bin_depths_time_offset_shifts_pressure_sampling():
    rdi = _make_rdi(nbin=4, nens=3)
    ctd = _make_ctd()
    tau_days = 30.0 / 86400.0
    z_shift, _ = assign_bin_depths(rdi, ctd, time_offset_days=tau_days)
    p_expected = np.interp(
        rdi.time_julian + tau_days, ctd.time_julian, ctd.pressure_dbar
    )
    np.testing.assert_allclose(z_shift, p_expected * 1.00445)


# --- estimate_ctd_adcp_lag (loadctd.m besttlag equivalent) ---


def _synthetic_cast(tau_s: float, dtctd_s: float = 0.5, nscans: int = 8000):
    """CTD pressure triangle cast + ADCP whose clock lags CTD by tau_s.

    ADCP timestamps t_a correspond to physical events at t_a + tau_s, so the
    aligned CTD pressure for ADCP sample t_a is p_ctd(t_a + tau_s):
    estimate_ctd_adcp_lag() should recover lagdt_days = +tau_s / 86400.
    """
    rng = np.random.default_rng(42)
    t0 = 2457000.0
    t_phys = t0 + np.arange(nscans) * dtctd_s / 86400.0
    half = nscans // 2
    # descent then ascent at ~1 m/s with some variability
    w_pkg = np.where(np.arange(nscans) < half, 1.0, -1.0)
    w_pkg = w_pkg * (1.0 + 0.3 * np.sin(np.arange(nscans) / 97.0))
    depth = np.cumsum(w_pkg) * dtctd_s
    depth -= depth.min() - 10.0
    ctd = CTDTimeSeries(
        time_julian=t_phys,
        pressure_dbar=depth / 1.00445,
        temp_c=np.full(nscans, 5.0),
        salinity=np.full(nscans, 35.0),
    )
    # ADCP samples every ~1.4 s, clock offset -tau (labels earlier than physical)
    idx = np.arange(0, nscans - 200, 3)
    t_adcp = t_phys[idx] - tau_s / 86400.0
    w_adcp = w_pkg[idx] + rng.normal(0.0, 0.05, idx.size)
    return t_adcp, w_adcp, ctd


def test_estimate_ctd_adcp_lag_recovers_synthetic_offset():
    tau_s = 12.0
    t_adcp, w_adcp, ctd = _synthetic_cast(tau_s)
    lag_scans, lagdt_days, corr = estimate_ctd_adcp_lag(
        t_adcp, w_adcp, ctd, lat_deg=-15.0
    )
    assert corr > 0.9
    assert abs(lagdt_days * 86400.0 - tau_s) <= 1.0  # within 2 CTD scans
    assert abs(lag_scans - round(tau_s / 0.5)) <= 1  # scan-quantized search


def test_estimate_ctd_adcp_lag_zero_offset():
    t_adcp, w_adcp, ctd = _synthetic_cast(0.0)
    lag_scans, lagdt_days, corr = estimate_ctd_adcp_lag(
        t_adcp, w_adcp, ctd, lat_deg=-15.0
    )
    assert corr > 0.9
    assert abs(lagdt_days * 86400.0) <= 1.0


# --- interval-based time fallback (I7N-style cnv with no time column) ---


def _make_ascii_cnv_no_time(tmp_path: Path) -> Path:
    """I7N 2018 raw-cnv variant: 24 Hz ASCII, no timeJ/timeS column; time is
    implicit via '# interval = seconds' + '# start_time' (NMEA)."""
    header = (
        "* Sea-Bird SBE 9 Data File:\n"
        "# nquan = 4\n"
        "# nvalues = 5\n"
        "# name 0 = prDM: Pressure, Digiquartz [db]\n"
        "# name 1 = t090C: Temperature [ITS-90, deg C]\n"
        "# name 2 = c0S/m: Conductivity [S/m]\n"
        "# name 3 = flag:  0.000e+00\n"
        "# interval = seconds: 0.0416667\n"
        "# bad_flag = -9.990e-29\n"
        "# start_time = Apr 28 2018 13:30:51 [NMEA time, header]\n"
        "*END*\n"
    )
    rows = "\n".join(
        f"{10.0 + i:.3f} {5.0:.4f} {3.2:.5f} 0.000e+00" for i in range(5)
    )
    p = tmp_path / "i7n.cnv"
    p.write_text(header + rows + "\n")
    return p


def test_interval_fallback_builds_time_axis(tmp_path):
    p = _make_ascii_cnv_no_time(tmp_path)
    ctd = load_ctd(p)
    assert len(ctd.time_julian) == 5
    # spacing == interval (per-step tolerance is float64 quantization at
    # absolute-Julian magnitude, ~5e-5 s; the mean is much tighter)
    dt_s = np.diff(ctd.time_julian) * 86400.0
    np.testing.assert_allclose(dt_s, 0.0416667, atol=5e-5)
    assert abs(dt_s.mean() - 0.0416667) < 1e-5
    # first scan at start_time (Apr 28 2018 13:30:51)
    expected0 = _to_julian(2018, 4, 28, 13 + 30 / 60.0 + 51 / 3600.0)
    assert abs(ctd.time_julian[0] - expected0) < 1e-9
    np.testing.assert_allclose(ctd.pressure_dbar, 10.0 + np.arange(5))


def test_no_time_and_no_interval_still_raises(tmp_path):
    header = (
        "# nquan = 2\n"
        "# name 0 = prDM: Pressure [db]\n"
        "# name 1 = t090C: Temperature\n"
        "*END*\n"
    )
    p = tmp_path / "noint.cnv"
    p.write_text(header + "10.0 5.0\n11.0 5.0\n")
    with pytest.raises(ValueError, match="time"):
        load_ctd(p)
