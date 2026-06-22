# RMSE Closure Investigation — P16N Cast 003

**Status:** Blocked on missing raw data. Resume when MATLAB intermediate arrays are available.

**Goal:** Close u_RMSE from 0.0718 m/s to < 0.050 m/s (the `test_inverse_u_rmse` / `test_inverse_v_rmse` integration test targets). The tests are currently `xfail`.

---

## Current State (2026-06-22)

| Metric | Python | MATLAB ref |
|---|---|---|
| u RMSE (full profile) | 0.0718 m/s | — |
| corr 0–500 m | +0.91 | — |
| corr 500–1000 m | +0.89 | — |
| corr 1000–1500 m | **−0.39** | — |
| corr 1500–2000 m | **−0.41** | — |
| n_se | 524 | 827 |

The 0–1000 m range matches well. A systematic anti-correlation appears at 1000–2000 m and degrades the overall RMSE.

---

## Root Cause: Super-Ensemble Count Discrepancy

MATLAB produces 827 super-ensembles; our Python produces 524. MATLAB's `getinv.m` line 42 computes `dz` dynamically:

```matlab
ps = setdefv(ps, 'dz', medianan(abs(diff(di.izm(:,1)))));
```

For an 8 m cell-size instrument this gives `dz ≈ 8 m`, not our hardcoded `dz = 16 m`. The smaller window means finer temporal resolution of `u_ctd`. Testing `prepare_superensembles(ens, dz=None)` (which already auto-computes from bin spacing) is untested — and RMSE didn't clearly improve from the dz=16 run.

---

## What Has Been Ruled Out

| Hypothesis | Test | Result |
|---|---|---|
| GPS constraint over-weighting | `barofac = 0.0` | Anti-corr persists (same corrs) |
| Bottom-track constraint | `botfac = 0.0` | Marginal change; not causal |
| `smoofac` difference | Both use 0 | Not a factor |
| Adaptive `velerr` (MATLAB: ~2.5 m/s) | Tested | RMSE 0.0718 → 0.0712 (negligible) |
| `zero_mask` excluding valid single-ping windows | Minor effect | RMSE not significantly changed |

The GPS sensitivity test (`barofac = 0.0`) is the key negative result: the anti-correlation is intrinsic to the observations/matrix structure, not to constraint weighting.

---

## What the Data Shows

Raw super-ensemble `ruav` (reference-bin averaged relative velocity) at z=959 m (downcast, SE 59) ≈ −0.052 m/s. The reference `u_ocean` at the reference bin water depth ≈ +0.08 m/s. This implies:

```
u_ctd_true = u_ocean_ref − ruav ≈ +0.08 − (−0.052) ≈ +0.13 m/s
```

MATLAB's saved `uctd` at equivalent depths is +0.05–0.07 m/s. The `ru` values are internally consistent (the observations themselves are not clearly wrong), but the inverse fails to recover the large positive `u_ctd` at 1000–1600 m depth.

---

## What Is Needed to Proceed

To distinguish between candidate root causes, we need MATLAB's intermediate arrays for cast 003:

1. **`di.ru` / `di.rv`** — raw super-ensemble relative velocities (before inversion). Compare with our `se.ru` bin-by-bin at the same ensemble indices. A discrepancy here means the super-ensemble formation (`prepare_superensembles`) is wrong.

2. **`di.izm`** — MATLAB's super-ensemble bin depths. Verify whether MATLAB's 827 SEs use 8 m windows or something else.

3. **MATLAB `uctd` vs `zctd`** — the full instrument velocity time series. The anti-correlation at 1000–2000 m implies MATLAB recovers large positive `u_ctd` (±0.10–0.15 m/s) at those depths while our solver cannot.

4. **The full A matrix structure** — specifically, which depth bins `j` each MATLAB observation maps to. If MATLAB's `dz` for the inverse is 8 m (not 10 m as in our `InverseParams`), the column assignments differ.

### How to Extract from MATLAB

In MATLAB, after running `getinv`:
```matlab
% di is the prepinv output
save('di_cast003.mat', 'di');
% dr is the getinv output  
save('dr_cast003.mat', 'dr');
```

Then `di.ru`, `di.izm`, `dr.uctd`, `dr.zctd` in the `.mat` files can be compared directly.

---

## Pending Hypothesis (Not Yet Tested)

**`dz = None` in `prepare_superensembles`** (auto-compute from bin spacing). This would give `dz ≈ 8 m` matching MATLAB's 827 SEs. The integration test currently hard-codes `dz=16.0` on line 192 of `tests/integration/test_inverse_p16n_cast003.py`. Testing the `dz=None` path may or may not help — but cannot be conclusively evaluated without MATLAB intermediates to compare against.

---

## Implementation Tasks (when data is available)

- [ ] Extract `di.ru`, `di.izm`, `dr.uctd`, `dr.zctd` from MATLAB for P16N cast 003
- [ ] Compare `di.ru[:, se_idx]` with `se.ru[:, se_idx]` bin-by-bin at depths 900–1200 m
- [ ] If `di.ru ≈ se.ru`: problem is in `compute_inverse`, not `prepare_superensembles`. Test different `params.dz` for the depth bins (8 m vs 10 m) and examine recovered `u_ctd` profile.
- [ ] If `di.ru ≠ se.ru`: problem is in `prepare_superensembles`. Likely causes: reference bin ordering, `izr` selection, `_medianan` vs MATLAB `medianan`, or UL bin reversal.
- [ ] Once root cause found: fix, run integration test, promote `xfail` tests if RMSE < 0.05 m/s.
