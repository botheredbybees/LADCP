"""Unit tests for SBE hex decoder."""
from __future__ import annotations
from pathlib import Path
import textwrap
import pytest

from ladcp.ingestion.sbe_hex import XmlconCoeffs, load_xmlcon, HexHeader, parse_hex_header


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
