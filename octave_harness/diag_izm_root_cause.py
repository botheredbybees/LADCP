"""Root-cause test for the depth-varying izm registration offset (REPORT.md P2).

Hypothesis (two mechanisms, both in loadctd.m, both absent from Python's
assign_bin_depths()):

  H1  pressure->depth formula: loadctd.m::p2z uses Saunders & Fofonoff (1976)
      with the cast latitude (p.poss(1) = -15); Python's assign_bin_depths()
      has the same formula in its lat_deg branch but NO CALLER PASSES lat_deg,
      so every pipeline run uses the shallow-water fallback z = p * 1.00445.
      Predicted signature: 0 at surface growing to ~ +90 m (Octave shallower)
      at 4400 dbar -- matches P2's depth-correlated component.

  H2  CTD-ADCP clock offset: loadctd.m detects lag = 47 CTD scans (~ -23.5 s,
      besttlag on w_ctd vs ADCP w; see recorded_p_struct_attrs.txt log),
      re-interpolates pressure at d.time_jul + lagdt AND shifts d.time_jul
      itself by lagdt. Python does neither. Predicted signature: depth error
      = descent rate x 23.5 s (~20 m early cast, ~0 at the bottom where w=0)
      -- matches P2's -20.8 m near-start component.

Method: diff Octave step09 d.z (instrument depth, negative-down, the
row-uniform component of izm) against Python z variants computed from the
same CTD file, nearest-time matched exactly like diff_stages.py:

  A  baseline        z = -1.00445 * p(t)                (current pipeline)
  B  A + p2z         z = -p2z(p(t), lat)                (H1 only)
  C  B + lag interp  z = -p2z(p(t + lagdt), lat)        (H1 + H2 pressure part)
  D  C + time shift  as C, and py time axis shifted by lagdt for matching
                     (H1 + H2, full loadctd.m emulation)

Pass criterion: A reproduces P2's stats (mean ~ +34 m, corr(diff, depth)
~ -0.80); D collapses to ~0 (residual << 10 m rms, no depth correlation).
"""
from pathlib import Path

import numpy as np
import scipy.io as sio

import sys

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from ladcp.ingestion.ctd import load_ctd  # noqa: E402
from ladcp.ingestion.rdi import load_rdi  # noqa: E402

DUMPS = REPO / "octave_harness" / "work" / "dumps"
DATA_DIR = REPO / "test_data" / "2015_P16N"

LAT_DEG = -15.0          # loadctd.m: lat = p.poss(1) = -15 (integer degrees part)
LAG_SCANS = 47           # recorded log: "best lag W: 47 CTD scans ~ -24 seconds"
DT_CTD_S = 0.5           # 2 Hz CTD file
LAGDT_DAYS = -(LAG_SCANS * DT_CTD_S) / 86400.0  # loadctd.m: lagdt = -lag*dtctd


def p2z(p_dbar: np.ndarray, lat: float) -> np.ndarray:
    """loadctd.m::p2z -- Saunders & Fofonoff (1976), EOS-80 refit. p in dbar."""
    p = p_dbar / 10.0  # to bars, as in the MATLAB source
    x = np.sin(np.radians(lat)) ** 2
    gr = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-5 * p
    depth = (((-1.82e-11 * p + 2.279e-7) * p - 2.2512e-3) * p + 97.2659) * p
    return depth / gr


def nearest_match(key_a, key_b):
    order = np.argsort(key_b)
    kb = key_b[order]
    idx = np.clip(np.searchsorted(kb, key_a), 1, len(kb) - 1)
    left = idx - 1
    use_right = np.abs(kb[idx] - key_a) < np.abs(kb[left] - key_a)
    nearest = np.where(use_right, idx, left)
    return order[nearest], np.abs(kb[nearest] - key_a)


def report(label, z_oct, z_py, t_err_s):
    diff = z_oct - z_py  # both negative-down; positive diff = Octave shallower
    finite = np.isfinite(diff)
    d = diff[finite]
    depth = z_oct[finite]
    corr = float(np.corrcoef(d, depth)[0, 1]) if d.size > 2 else float("nan")
    print(
        f"{label:<44} mean={d.mean():>7.2f}  median={np.median(d):>7.2f}  "
        f"rms={np.sqrt((d**2).mean()):>7.2f}  min={d.min():>7.2f}  "
        f"max={d.max():>7.2f}  corr(diff,depth)={corr:>6.3f}  "
        f"[match dt: mean {t_err_s.mean()*86400:.1f} s]"
    )


def main() -> None:
    step9 = sio.loadmat(DUMPS / "step09.mat", struct_as_record=False, squeeze_me=True)
    d9 = step9["d"]
    oct_time = np.asarray(d9.time_jul, dtype=float)
    z_oct = np.asarray(d9.z, dtype=float)  # negative-down instrument depth

    rdi = load_rdi(DATA_DIR / "003DL000.000")
    ctd = load_ctd(DATA_DIR / "003_01.cnv")
    t_adcp = rdi.time_julian

    p_at = lambda t: np.interp(t, ctd.time_julian, ctd.pressure_dbar)  # noqa: E731

    variants = {
        "A baseline (p*1.00445, no lag)": (-1.00445 * p_at(t_adcp), t_adcp),
        "B p2z(lat), no lag": (-p2z(p_at(t_adcp), LAT_DEG), t_adcp),
        "C p2z(lat) + lag-shifted pressure": (
            -p2z(p_at(t_adcp + LAGDT_DAYS), LAT_DEG), t_adcp),
        "D = C + lag-shifted time axis (full emu)": (
            -p2z(p_at(t_adcp + LAGDT_DAYS), LAT_DEG), t_adcp + LAGDT_DAYS),
    }

    print(f"Octave step09 d.z: n={z_oct.size}, "
          f"range [{z_oct.min():.1f}, {z_oct.max():.1f}] m")
    print(f"lagdt = {LAGDT_DAYS*86400:+.1f} s, lat = {LAT_DEG} deg\n")
    print("diff = z_oct - z_py (negative-down: positive diff = Octave SHALLOWER)")
    for label, (z_py, t_py) in variants.items():
        idx, t_err = nearest_match(oct_time, t_py)
        report(label, z_oct, z_py[idx], t_err)

    # --- residual structure: is what's left of variant B a pure time offset? ---
    # If diff_B = w * tau (descent rate x clock offset), fitting tau and
    # re-interpolating pressure at t + tau should collapse the residual.
    idx, t_err = nearest_match(oct_time, t_adcp)
    z_b = -p2z(p_at(t_adcp), LAT_DEG)[idx]
    diff_b = z_oct - z_b
    w_oct = np.gradient(z_oct, oct_time * 86400.0)  # m/s, package velocity
    ok = np.isfinite(diff_b) & np.isfinite(w_oct)
    tau = float(np.sum(diff_b[ok] * w_oct[ok]) / np.sum(w_oct[ok] ** 2))
    r_wd = float(np.corrcoef(diff_b[ok], w_oct[ok])[0, 1])
    print(f"\nresidual-B vs descent rate: corr={r_wd:.3f}, "
          f"least-squares tau = {tau:+.1f} s (diff ~ w*tau)")

    z_e = -p2z(p_at(t_adcp + tau / 86400.0), LAT_DEG)
    idx_e, t_err_e = nearest_match(oct_time, t_adcp)
    report(f"E p2z(lat) + fitted tau={tau:+.1f}s pressure", z_oct, z_e[idx_e], t_err_e)


if __name__ == "__main__":
    main()
