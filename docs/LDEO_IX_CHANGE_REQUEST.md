# Change request: LDEO_IX LADCP processing software

**Prepared by:** Australian Antarctic Division / RSV *Nuyina* science data
systems, 2026-07-11.
**Against:** LDEO_IX `Version IX_14beta` as vendored in this repository
(`docs/legacy/`). The reference dataset used for verification (GO-SHIP
P16N 2015, cast 003) was originally processed with `IX_13beta`; all
findings below were confirmed present in the `IX_14beta` sources.
**How found:** during a line-by-line reimplementation of the LDEO_IX
horizontal-velocity pipeline in Python, validated by a differential
harness that runs the unmodified MATLAB code under GNU Octave 9.2 and
diffs every processing stage against the Python port (see
`octave_harness/REPORT.md`, sections M1-M4 and P1-P6). The port now
reproduces LDEO_IX to machine precision through super-ensemble formation,
which required understanding — and in the cases below, questioning — the
original code's exact behavior.

**Audience:** intended for submission upstream (A.M. Thurnherr / LDEO)
and as a permanent record for anyone processing data with LDEO_IX.

---

## Summary

| # | File / lines | Severity | Type | One-line description |
|---|---|---|---|---|
| 1 | `loadrdi.m` 475-478, 497-500 | Medium | Correctness | Non-pinging check tests the v-gradient twice, never the u-gradient (misnamed variables + copy-paste) |
| 2 | `loadnav.m` 72 vs 137/196 | Medium | Correctness | `nav_time_base` default set on `p` but consumed from `f` — the default is silently dead |
| 3 | `prepinv.m` 613 + `outlier.m` 47 | Low-Medium | Probable dead code | Bottom-track outlier editing of super-ensembles can never execute (orientation check always false) |
| 4 | `sounds.m` header | Low | Documentation | Documented check value (1731.995) does not hold for the code as written (actual 1732.139394) |
| 5 | `getinv.m` 352/722, `plotraw.m` 223 | Low | Portability | `do` used as identifier — reserved keyword in Octave |
| 6 | `end_processing_step.m` 25 vs `begin_processing_step.m` 25 | Low | Portability | Checkpoint saved without `.mat`, loaded with `.mat` — resume broken under Octave |
| 7 | `plotraw.m` (dependency) | Low | Packaging | `makebars.m` referenced but absent from the distribution |
| 8 | `medianan.m` usage in `prepinv.m` 530-551 | Info | Documentation | `medianan(x, round(n/2))` degenerates to a plain NaN-mean — name misleads readers/porters |

None of these findings invalidate LDEO_IX's scientific results on the
data we tested: items 1-3 affect edge-case editing paths whose practical
effect on our validation cast was negligible, and items 4-8 are
documentation/portability. They are reported because (a) the *intent* of
the code is clearly different from its behavior in items 1-3, and (b)
items 5-7 block running the software under Octave, which is increasingly
how groups without MATLAB licenses will run it.

---

## 1. `loadrdi.m`: non-pinging detection never tests the u-gradient

**Location:** lines 475-478 (downlooker) and 497-500 (uplooker).

```matlab
drw=medianan(abs(diff(d.rw(d.izd,:))));
dru=medianan(abs(diff(d.rv(d.izd,:))));     % <- named dru, computed from rv
drv=medianan(abs(diff(d.ru(d.izd,:))));     % <- named drv, computed from ru
nbad=find(abs(drw)<0.005 & abs(dru)<0.005 & abs(dru)<0.005);
%                                            ^^^^^^^^^ dru tested twice
```

Two compounding defects:

1. The variable names are swapped relative to their contents (`dru` holds
   the **v**-gradient, `drv` holds the **u**-gradient).
2. The condition tests `dru` twice and `drv` never — so the u-velocity
   bin-to-bin gradient plays no part in the "dead instrument / not
   pinging" detection. Only the w- and v-gradients are effective.

The identical pattern appears in the uplooker block (lines 497-500).

**Impact:** an ensemble with flat w and v but structured u would be
flagged as non-pinging and have its weights NaN'd despite carrying real
data; conversely the intended three-component test is weaker by one
component. On P16N 003 this path flagged 17 uplooker ensembles; testing
all three components does not change that count for this cast, so the
practical impact there is nil — but the code plainly does not implement
its evident intent.

**Suggested fix:**

```matlab
drw=medianan(abs(diff(d.rw(d.izd,:))));
drv=medianan(abs(diff(d.rv(d.izd,:))));
dru=medianan(abs(diff(d.ru(d.izd,:))));
nbad=find(abs(drw)<0.005 & abs(dru)<0.005 & abs(drv)<0.005);
```

(and the same in the uplooker block).

**Note for reproducibility:** our Python port deliberately replicates the
buggy behavior (only w- and v-gradients tested) to remain bit-compatible
with existing LDEO_IX output; see `ladcp.qa.editing.build_ldeo_weights`.
If this fix is adopted upstream, ports should follow.

## 2. `loadnav.m`: `nav_time_base` default is set on the wrong struct

**Location:** line 72 sets the default; lines 137 and 196 consume it.

```matlab
p = setdefv(p,'nav_time_base',0);    % line 72: default onto p
...
switch f.nav_time_base               % line 137: read from f
...
if f.nav_time_base ~= 0              % line 196: read from f
```

The default is placed on the parameter struct `p`, but every consumer
reads `f.nav_time_base`. Unless the user explicitly sets
`f.nav_time_base` in their cast-parameters file, the `switch` errors (or
in permissive interpreters, misbehaves).

**Precedent:** `loadctd.m`'s changelog records the *identical* bug for
`ctd_time_base`, fixed in 2014 ("moved ctd_time_base from p. to f."). The
fix was never ported to `loadnav.m`.

**Suggested fix:** `f = setdefv(f,'nav_time_base',0);` on line 72,
mirroring the 2014 `loadctd.m` fix.

## 3. `prepinv.m`: bottom-track outlier editing of super-ensembles is dead code

**Location:** `prepinv.m` line 613 calls `[di,p]=outlier(di,p)`;
`outlier.m` line 47 gates its bottom-track branch on:

```matlab
if size(dummyb,2)==4, ibvel=1; else, ibvel=0; end
```

When `outlier()` is called from `loadrdi.m` (its original context),
`d.bvel` is `(n_ens, 4)` and the branch runs. When called from
`prepinv.m`, `di.bvel` is `(4, n_se)` — components along dim 1 — so
`size(...,2)` equals the number of super-ensembles, the check is false,
and **the bottom-track outlier pass on super-ensembles never executes**,
silently, on every cast.

**Impact:** if the check is an intentional "only edit BT when a 4-column
bvel is present" guard, the code works by coincidence of orientation and
deserves a comment; if BT super-ensemble editing was intended (the
symmetric treatment of water-column and BT data elsewhere in `outlier.m`
suggests it was), it has never run. Either way the current form is a trap
for maintainers and porters — we spent nontrivial harness time proving
which behavior the reference outputs actually contain (they contain the
skip: our validation only matched Octave's step-10 dumps to machine
precision once our port *disabled* BT editing in this call).

**Suggested fix:** make the intent explicit — either
`[di,p]=outlier(di,p);  % NB: BT branch intentionally inactive here` or
transpose-aware handling inside `outlier.m`.

## 4. `sounds.m`: documented check value does not hold for the code as written

**Location:** header comment, lines 15-16:

```matlab
% Checkvalues:
%  SVEL=1731.995 :Salinity=40.0, Temp.=40.0, Pres.=10000.0
```

Running the unmodified function under Octave 9.2:

```
>> sounds(10000.0, 40.0, 40.0)
ans = 1732.139394
```

i.e. 0.144 m/s off the documented UNESCO/Chen & Millero check value. We
did not isolate which coefficient drifted in the FORTRAN→MATLAB
translation; the ~1×10⁻⁴ relative error is immaterial where only the
ratio of sound speeds enters (the `getdpthi.m` corrections), but the
function's absolute output fails its own documented acceptance test,
which matters if it is ever reused for absolute sound-speed work.

**Suggested fix:** either correct the coefficients against UNESCO 44 /
Chen & Millero (1977) so the check value passes, or amend the comment to
record the actual output of the implementation. Ports that must match
existing LDEO_IX output byte-for-byte need the *current* arithmetic (our
port pins 1732.139394 in its tests for exactly this reason).

## 5. `do` used as an identifier — Octave reserved keyword

**Location:** `getinv.m` line 352 (`do=d;`, later consumed at line 722,
`de.do=do;`) and `plotraw.m` line 223 (`function checkbeam(t,ax,do)`).

`do` is a reserved keyword in GNU Octave (do-until loops). Both files
fail to *parse* under Octave, which blocks the entire pipeline (not just
plotting): `getinv.m` is the inverse solver.

**Suggested fix:** rename (`d_orig` / `is_bottom` are the names our
harness used; any non-keyword works). No behavioral change in MATLAB.

## 6. Checkpoint save/load filename mismatch — resume broken under Octave

**Location:** `end_processing_step.m` line 25 saves without an extension:

```matlab
eval(sprintf('save %s_%d',f.checkpoints,pcs.cur_step));
```

while `begin_processing_step.m` line 25 (with its Jul 2016 changelog
entry "added .mat to checkpoint filename") loads with one:

```matlab
load(sprintf('%s_%d.mat',f.checkpoints,pcs.cur_step-1));
```

MATLAB's `save` auto-appends `.mat`, so the pair works there. Octave's
`save` does not, so every checkpoint load fails and `p.checkpoints`-based
resume silently never works under Octave.

**Suggested fix:** add `.mat` explicitly in `end_processing_step.m`
(matching the 2016 change to the load side).

## 7. `makebars.m` missing from the distribution

`plotraw.m` calls `makebars()` (diagnostic bar overlay), but the function
is absent from the IX_14beta distribution we vendored (checked the full
tree including the original tarball). Under MATLAB with no toolbox
providing it, `plotraw` errors. Cosmetic (plotting only), but it makes a
default full run crash.

## 8. (Informational) `medianan(x, na)` with `na = round(n/2)` is a plain mean

`prepinv.m` lines 530-551 compute per-window averages via
`medianan(..., iav)` with `iav = round(length(ur)/200*p.avpercent)` and
`avpercent = 100`, i.e. `na = round(n_win/2)`. `medianan.m` treats `na`
as a *half*-window — it averages the `2*na+1` central sorted values,
clipped to the available range — and `2*round(n/2)+1 ≥ n+1` always covers
every finite sample. The expression is therefore exactly a NaN-mean for
the default configuration; the sort and median machinery contribute
nothing.

Not a bug (the generality presumably serves `avpercent < 100`
configurations), but the naming cost us a real defect in our port: a
reasonable reading of "medianan with an averaging count" as a *trimmed*
mean is wrong on every cell. A one-line comment at the call sites
("avpercent=100 ⇒ this is meannan") would spare future porters and
reviewers the same trap.

---

## Verification notes

All items were confirmed against the vendored sources by inspection, and
items 1, 3, 4 additionally by execution under GNU Octave 9.2.0 via the
differential harness in `octave_harness/` (whose stage dumps also serve
as the reference for our Python port's parity tests: velocities and
weights agree with the unmodified LDEO_IX code to ≤5×10⁻⁸ m/s rms through
super-ensemble formation on GO-SHIP P16N 2015 cast 003). The two `do`
renames (item 5) and the checkpoint-extension fix (item 6) are the only
modifications our harness had to make to the M-code to run it at all;
they are recorded with diffs in `octave_harness/PATCHES.md`.
