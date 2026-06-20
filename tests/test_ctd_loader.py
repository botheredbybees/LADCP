import numpy as np
import pytest
from ladcp.ingestion.ctd import CTDTimeSeries, _parse_sbe_header, _map_column

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


import numpy as np
from pathlib import Path
from ladcp.ingestion.ctd import load_ctd

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
