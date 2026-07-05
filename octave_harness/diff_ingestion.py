"""M1 -- compare LDEO_IX's loadrdi.m ingestion (run under Octave, dumped to
octave_harness/dumps/m1_loadrdi.mat) against our own ladcp.ingestion.rdi.load_rdi()
on the same two raw files (P16N cast 003 DL/UL PD0).

Two kinds of fields are compared:

1. Static per-instrument config (fixed-leader-derived): exact match expected,
   no alignment needed.
2. Per-ensemble time series (time, heading, pitch, roll, temperature, sound
   velocity): loadrdi.m's updown() merges DL+UL onto a single joint ensemble
   axis via an integer lag shift (see the recorded log: "shift ADCP timeseries
   by lag: 1"). Rather than assume that shift, we find it empirically per
   instrument by nearest-time matching against our own per-file time_julian,
   then diff the matched pairs.

Velocities (d.ru/rv/rw/re) are NOT compared here: loadrdi.m applies the
beam->earth rotation inline during ingestion (this PD0 is in BEAM
coordinates -- see the log's "DETECTED BEAM coordinates: rotating to EARTH
coordinates"), whereas our architecture keeps ingestion (raw beam frame) and
transforms (src/ladcp/transforms) as separate layers. A fair velocity
comparison requires our own beam2earth() output and belongs in M3, where the
full pipeline stage context (bin-mapping, tilt correction) is available.
"""
from pathlib import Path

import numpy as np
import scipy.io as sio

from ladcp.ingestion.rdi import load_rdi

REPO = Path(__file__).resolve().parent.parent
DUMP_PATH = REPO / "octave_harness" / "dumps" / "m1_loadrdi.mat"
DL_PATH = REPO / "test_data" / "2015_P16N" / "003DL000.000"
UL_PATH = REPO / "test_data" / "2015_P16N" / "003UL000.000"


def _load_octave_dump():
    m = sio.loadmat(DUMP_PATH, struct_as_record=False, squeeze_me=True)
    return m["d"], m["p"]


def _config_table(d, p, dl, ul) -> list[tuple]:
    rows = []
    checks = [
        ("nbin", p.nbin_d, dl.nbin, p.nbin_u, ul.nbin),
        ("blen_m", p.blen_d, dl.blen_m, p.blen_u, ul.blen_m),
        ("blnk_m", p.blnk_d, dl.blnk_m, p.blnk_u, ul.blnk_m),
        ("dist_m", p.dist_d, dl.dist_m, p.dist_u, ul.dist_m),
        ("beam_angle_deg", p.beamangle, dl.beam_angle_deg, p.beamangle, ul.beam_angle_deg),
    ]
    for name, oct_dl, py_dl, oct_ul, py_ul in checks:
        rows.append((name, "DL", float(oct_dl), float(py_dl), float(oct_dl) - float(py_dl)))
        rows.append((name, "UL", float(oct_ul), float(py_ul), float(oct_ul) - float(py_ul)))
    return rows


def _nearest_time_match_and_diff(oct_time, oct_field, py_time, py_field, angular=False):
    """Match each Octave ensemble to the nearest-in-time Python ensemble
    (per-sample nearest-neighbor, not a fixed index shift) -- necessary
    because ping intervals are non-constant/staggered in this cast (see the
    recorded log's "non-constant ping rate ... staggered pinging?").
    """
    order = np.argsort(py_time)
    py_time_sorted = py_time[order]
    py_field_sorted = py_field[order]

    idx = np.searchsorted(py_time_sorted, oct_time)
    idx = np.clip(idx, 1, len(py_time_sorted) - 1)
    left = idx - 1
    use_right = np.abs(py_time_sorted[idx] - oct_time) < np.abs(py_time_sorted[left] - oct_time)
    nearest = np.where(use_right, idx, left)

    matched_time = py_time_sorted[nearest]
    matched_field = py_field_sorted[nearest]
    time_err = np.abs(matched_time - oct_time)

    diff = oct_field - matched_field
    if angular:
        diff = (diff + 180.0) % 360.0 - 180.0

    return {
        "n_matched": len(oct_time),
        "mean_time_err_days": float(np.mean(time_err)),
        "max_time_err_days": float(np.max(time_err)),
        "max_abs_diff": float(np.nanmax(np.abs(diff))),
        "rms_diff": float(np.sqrt(np.nanmean(diff**2))),
    }


def main() -> None:
    d, p = _load_octave_dump()
    dl = load_rdi(DL_PATH)
    ul = load_rdi(UL_PATH)

    print("=" * 78)
    print("M1 ingestion diff: LDEO_IX loadrdi.m (Octave) vs ladcp.ingestion.rdi.load_rdi()")
    print("=" * 78)

    print("\n--- Static config fields (exact match expected) ---")
    print(f"{'field':<16}{'inst':<5}{'octave':>14}{'python':>14}{'diff':>12}")
    for name, inst, ov, pv, diff in _config_table(d, p, dl, ul):
        print(f"{name:<16}{inst:<5}{ov:>14.6g}{pv:>14.6g}{diff:>12.3g}")

    oct_time = np.asarray(d.time_jul, dtype=float)
    oct_hdg = np.asarray(d.hdg, dtype=float)
    oct_pit = np.asarray(d.pit, dtype=float)
    oct_rol = np.asarray(d.rol, dtype=float)
    oct_temp = np.asarray(d.temp, dtype=float)
    oct_sv = np.asarray(d.sv, dtype=float)

    print("\n--- Per-ensemble time series (nearest-time match, per-sample) ---")
    print(f"{'field':<14}{'inst':<5}{'n':>7}{'mean_dt(d)':>12}{'max_dt(d)':>11}{'max|diff|':>12}{'rms diff':>12}")
    fields = [
        ("time_julian", oct_time, oct_time, dl.time_julian, ul.time_julian, False),
        ("heading_deg", oct_hdg[0], oct_hdg[1], dl.heading, ul.heading, True),
        ("pitch_deg", oct_pit[0], oct_pit[1], dl.pitch, ul.pitch, False),
        ("roll_deg", oct_rol[0], oct_rol[1], dl.roll, ul.roll, False),
        ("temp_c", oct_temp[0], oct_temp[1], dl.temp_c, ul.temp_c, False),
        ("sound_vel_ms", oct_sv[0], oct_sv[1], dl.sound_vel_ms, ul.sound_vel_ms, False),
    ]
    for name, oct_dl_field, oct_ul_field, py_dl_field, py_ul_field, angular in fields:
        for inst, oct_field, py_field, py_time in [
            ("DL", oct_dl_field, py_dl_field, dl.time_julian),
            ("UL", oct_ul_field, py_ul_field, ul.time_julian),
        ]:
            r = _nearest_time_match_and_diff(oct_time, oct_field, py_time, py_field, angular=angular)
            print(
                f"{name:<14}{inst:<5}{r['n_matched']:>7}"
                f"{r['mean_time_err_days']:>12.2e}{r['max_time_err_days']:>11.2e}"
                f"{r['max_abs_diff']:>12.4g}{r['rms_diff']:>12.4g}"
            )

    print(
        "\nVelocities (d.ru/rv/rw/re) intentionally NOT compared here -- Octave's "
        "loadrdi.m applies the beam->earth rotation inline (this PD0 is BEAM-coordinate), "
        "our load_rdi() keeps raw beam-frame data (ingestion/transforms are separate "
        "layers). See M3 for the fair post-transform comparison."
    )


if __name__ == "__main__":
    main()
