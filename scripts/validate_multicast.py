"""Multi-cast pipeline validation vs archived LDEO_IX outputs.

For every reference NetCDF in --ref-dir (LDEO_IX ladcp2cdf output, named
NNN.nc), find the matching raw DL/UL PD0 and CTD time-series files in
--raw-dir, run ladcp.pipeline.process_cast with parameters read from the
reference file's own attributes, and report per-cast u/v RMSE (and
correlation) against the archived profiles on well-observed bins
(nvel >= 3), plus a summary.

Anti-overfitting context: the pipeline was tuned against exactly one cast
(P16N 2015 003). This script exists to measure it on casts and cruises it
has never seen. When the archived run used a SADCP constraint, ladcp2cdf
embedded the exact profile in the reference NC (z_sadcp/u_sadcp/...) and
CastParams.from_ldeo_nc reads it back, so the rerun applies the same
constraint set as the archive; uship/vship barotropic constraints are
likewise read from the reference attrs.

Usage:
  uv run python scripts/validate_multicast.py \
      --raw-dir test_data/2018_I7N/raw --ref-dir test_data/2018_I7N/processed_uv \
      [--casts 001,002,...] [--dl-glob "{cast}DL*.000"] [--ul-glob "{cast}UL*.000"] \
      [--ctd-glob "{cast}*.cnv"] [--out results.csv]

The globs are .format()-ed with {cast} = the 3-digit cast id; adjust to
the cruise's naming convention.
"""
from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ladcp.pipeline import CastParams, process_cast  # noqa: E402


def rmse_vs_ref(res, ref_path: Path) -> dict:
    ds = netCDF4.Dataset(str(ref_path))
    rz = np.array(ds.variables["z"][:])
    ru = np.array(ds.variables["u"][:])
    rv = np.array(ds.variables["v"][:])
    rn = np.array(ds.variables["nvel"][:]) if "nvel" in ds.variables \
        else np.full(rz.shape, 99)
    ds.close()
    iu = np.interp(rz, res.z, res.u, left=np.nan, right=np.nan)
    iv = np.interp(rz, res.z, res.v, left=np.nan, right=np.nan)
    ok = np.isfinite(ru) & np.isfinite(iu) & (rn >= 3)
    n = int(ok.sum())
    if n < 10:
        return dict(n=n, u_rmse=np.nan, v_rmse=np.nan, r_u=np.nan)
    return dict(
        n=n,
        u_rmse=float(np.sqrt(np.mean((iu[ok] - ru[ok]) ** 2))),
        v_rmse=float(np.sqrt(np.mean((iv[ok] - rv[ok]) ** 2))),
        r_u=float(np.corrcoef(iu[ok], ru[ok])[0, 1]),
        max_depth=float(rz[ok].max()),
    )


def find_one(raw_dir: Path, pattern: str, cast: str) -> Path | None:
    hits = sorted(raw_dir.glob(pattern.format(cast=cast)))
    return hits[0] if hits else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", type=Path, required=True)
    ap.add_argument("--ref-dir", type=Path, required=True)
    ap.add_argument("--casts", default=None,
                    help="comma-separated cast ids; default = all refs")
    ap.add_argument("--dl-glob", default="{cast}DL*.000")
    ap.add_argument("--ul-glob", default="{cast}UL*.000")
    ap.add_argument("--ctd-glob", default="{cast}*.cnv")
    ap.add_argument("--out", type=Path, default=None, help="CSV output path")
    args = ap.parse_args()

    refs = sorted(args.ref_dir.glob("[0-9][0-9][0-9].nc"))
    if args.casts:
        wanted = {c.strip().zfill(3) for c in args.casts.split(",")}
        refs = [r for r in refs if r.stem in wanted]

    rows = []
    print(f"{'cast':>5} {'n':>5} {'u_rmse':>8} {'v_rmse':>8} {'r(u)':>7} "
          f"{'maxz':>7}  status")
    for ref in refs:
        cast = ref.stem
        dl = find_one(args.raw_dir, args.dl_glob, cast)
        ul = find_one(args.raw_dir, args.ul_glob, cast)
        ctd = find_one(args.raw_dir, args.ctd_glob, cast)
        if not (dl and ul and ctd):
            missing = [n for n, p in (("DL", dl), ("UL", ul), ("CTD", ctd))
                       if p is None]
            print(f"{cast:>5} {'-':>5} {'-':>8} {'-':>8} {'-':>7} {'-':>7}  "
                  f"SKIP (missing {','.join(missing)})")
            continue
        try:
            params = CastParams.from_ldeo_nc(ref)
            res = process_cast(dl, ul, ctd, params)
            st = rmse_vs_ref(res, ref)
            rows.append(dict(cast=cast, **st))
            print(f"{cast:>5} {st['n']:>5} {st['u_rmse']:>8.4f} "
                  f"{st['v_rmse']:>8.4f} {st['r_u']:>7.3f} "
                  f"{st.get('max_depth', float('nan')):>7.0f}  ok")
        except Exception as e:  # keep going; report at the end
            rows.append(dict(cast=cast, n=0, u_rmse=np.nan, v_rmse=np.nan,
                             r_u=np.nan, error=repr(e)))
            print(f"{cast:>5} {'-':>5} {'-':>8} {'-':>8} {'-':>7} {'-':>7}  "
                  f"ERROR {type(e).__name__}: {e}")
            traceback.print_exc(limit=3)

    good = [r for r in rows if np.isfinite(r.get("u_rmse", np.nan))]
    print(f"\ncasts attempted: {len(rows)}   succeeded: {len(good)}")
    if good:
        u = np.array([r["u_rmse"] for r in good])
        v = np.array([r["v_rmse"] for r in good])
        print(f"u RMSE: median {np.median(u):.4f}  mean {u.mean():.4f}  "
              f"p90 {np.percentile(u, 90):.4f}  worst {u.max():.4f}")
        print(f"v RMSE: median {np.median(v):.4f}  mean {v.mean():.4f}  "
              f"p90 {np.percentile(v, 90):.4f}  worst {v.max():.4f}")
        print(f"casts under 0.05/0.05: "
              f"{sum(1 for r in good if r['u_rmse'] < .05 and r['v_rmse'] < .05)}"
              f"/{len(good)}")

    if args.out:
        import csv
        keys = ["cast", "n", "u_rmse", "v_rmse", "r_u", "max_depth", "error"]
        with open(args.out, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            wr.writeheader()
            wr.writerows(rows)
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
