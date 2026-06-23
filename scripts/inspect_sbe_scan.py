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
