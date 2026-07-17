# Bulk validation report — I7N 2018 (final) & A16N 2013 (final)

_Last updated 2026-07-17 local. Produced per BULK_VALIDATION_BRIEF.md:
raw casts from NCEI, validated with `scripts/validate_multicast.py` against
archived LDEO_IX outputs. Per-batch CSVs/logs in `test_data/<cruise>/bulk/`;
concatenated results in each cruise's `validate_all.csv`._

**Provenance note:** initial Haiku-agent runs suffered process-duplication
and path bugs; all I7N results below come from a single serial pipeline
run over Content-Length-verified downloads (batches 01–03 from the
verified agent pipeline of 2026-07-11, batches 04–13 from the
deterministic rerun of 2026-07-12). A16N agent-era batch CSVs were
quarantined to `test_data/2013_A16N/bulk/agent_era_csvs/` (produced while
concurrent downloads held file locks — untrusted); all 95 A16N casts were
then re-run from verified files in a single serial pipeline (batches
01–10, completed 2026-07-12 14:05). The re-run's per-cast RMSE values are
numerically identical to the quarantined agent-era run to the digits
reported, so the agent-era numbers are corroborated rather than
contradicted — the quarantine was the correct caution, but no corruption
was actually present in this cruise's results.

## Cruise 1: I7N 2018 — COMPLETE, 124/124 casts

**Counts:** 124 attempted, 124 succeeded, 0 harness errors, 0 missing.

| tier | casts |
|---|---|
| pass both (u<0.05 and v<0.05) | **53** |
| pass u only | 20 |
| marginal (u 0.05–0.2) | 41 |
| bad (u 0.2–1) | 0 |
| **EXPLODED (u RMSE > 1, up to ~10⁹)** | **10** |

**Distributions (all 124):** u median 0.0427, p90 0.1172; v median
0.0439, p90 0.2562. **Excluding the 10 exploded casts:** u median
0.0412; deep casts (>4000 m, n=81) median u 0.0469 vs shallow 0.0312.

**Exploded casts** (numerical blowup, RMSE ~10⁶–10¹⁰ — solver failure,
not disagreement): 018, 042, 060, 062, 086, 099, 102, 103, 118, 119.
Max depths 3440–4560 m (mid-pack, not the deepest). The distribution is
bimodal — no cast lands between u 0.2 and u 1 — indicating an
ill-conditioning/degeneracy failure mode in the inverse for specific
casts rather than gradual degradation. **This is the headline new
finding of the bulk run** (never seen in the 2-cast spot checks) and a
debuggable target: likely near-empty superensembles or a degenerate
constraint row. Investigate one exploded cast (e.g. 060) end-to-end.

**Interpretation vs known leads:** the 41 marginal casts are dominated
by v-misses (v is the systematic limiter, consistent with the
lanarrow/ps.shear porting queue). The u median at 0.0427 on a
124-cast unseen cruise — with SADCP constraints reconstructed from the
reference files — is strong evidence the core pipeline generalizes.

## Cruise 2: A16N 2013 — COMPLETE, 95/95 casts (001–095; 096+ excluded, no uplooker)

**Counts:** 95 attempted, 95 succeeded, 0 harness errors, 0 missing, 0
download-failed.

| tier | casts |
|---|---|
| pass both (u<0.05 and v<0.05) | **15** |
| pass u only (v ≥ 0.05) | 4 |
| pass u regardless of v | 19 |

**Distributions (all 95 succeeded):** u median 0.3353, p90 1.1292, worst
3 940 694.16 (cast 085). v median 0.3308, p90 0.9778, worst 786 403.12
(cast 085).

**Errored casts:** none.

**Top 10 worst casts by u_rmse:**

| Cast | n | u_rmse | v_rmse | r_u | max_depth (m) |
|------|---|--------|--------|-----|---------------|
| 085 | 670 | 3940694.1628 | 786403.1229 | 0.0821 | 5491.5 |
| 078 | 638 | 1.9132 | 1.0552 | -0.2476 | 5223.6 |
| 072 | 638 | 1.7521 | 2.3886 | 0.3812 | 5220.4 |
| 089 | 650 | 1.3365 | 0.7011 | 0.3084 | 5331.0 |
| 071 | 643 | 1.3224 | 1.6695 | -0.6550 | 5258.3 |
| 069 | 638 | 1.3078 | 0.5838 | -0.1098 | 5217.9 |
| 077 | 633 | 1.2638 | 0.7364 | 0.5312 | 5180.4 |
| 062 | 629 | 1.1941 | 0.7144 | 0.2719 | 5140.1 |
| 064 | 649 | 1.1308 | 0.4740 | -0.0126 | 5304.2 |
| 082 | 657 | 1.1292 | 1.4601 | 0.3843 | 5381.2 |

**Observations:**

- Deep casts (>4000 m): 59 succeeded; all 59 have u_rmse ≥ 0.05. Per the
  brief's calibration notes, large RMSE on deep A16N casts is the
  expected/known result (pending the ps.shear port), not a setup error.
- Shallow/mid casts (≤4000 m): 36 succeeded; 19 pass u_rmse < 0.05.
- 22 casts have negative u correlation (r_u < 0): 013, 033, 035, 038,
  049, 050, 051, 053, 064, 065, 066, 069, 071, 074, 078, 079, 081, 086,
  087, 092, 093, 094 — all among the deep, failing casts, consistent
  with a systematic (not just noisy) breakdown at depth rather than
  scatter around a good fit.
- Cast 085 is a single outlier by ~6 orders of magnitude (u_rmse ≈
  3.9×10⁶ vs. the next-worst cast at u_rmse ≈ 1.9) while its correlation
  (r_u = 0.08) is unremarkable — the same qualitative signature as the
  10 "EXPLODED" I7N casts (numerical blowup, not gradual disagreement).
  This is one cast, not the bimodal cluster seen in I7N, but it suggests
  the same ill-conditioning failure mode is cruise-independent and worth
  folding into that investigation rather than treating as an A16N-only
  quirk.
- RuntimeWarnings and CTD-lag warnings observed in harness logs are
  normal per the brief; no `ERROR` rows occurred in any batch log
  (confirmed by grep over all 10 logs).
