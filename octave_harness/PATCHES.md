# Patches to copied LDEO_IX M-code

Per CONTINUATION_PLAN.md's hard rules: `docs/legacy/` is untouched (read-only).
Everything here applies to the copy in `octave_harness/ldeo_ix/` only. Every
edit is listed below; nothing else was changed.

## Version note

The in-tree code (`octave_harness/ldeo_ix/default.m`) reports
`Version IX_14beta`. The recorded p-struct in `test_data/2015_P16N/003.nc`
(`GEN_Software_orig`) shows LDEO actually processed this cast with
`Version IX_13beta`. This mismatch is **not** fixed here per the plan's
explicit instruction -- noted for the report only.

## Genuine bugs in the in-tree code (fixed to get the pipeline running)

1. **`getinv.m`**: `do` used as a plain variable name (`do=d;` / `de.do=do;`).
   `do` is a reserved keyword in Octave (`do ... until`), a valid identifier
   in MATLAB. Renamed the *variable* to `d_orig` (the struct field name
   `de.do` is untouched -- field names aren't keyword-restricted). Two-line
   change, no behavior difference.

2. **`plotraw.m`**: same class of bug -- `function checkbeam(t,ax,do)` uses
   `do` as a function parameter. Renamed to `is_bottom` (call sites use
   positional args only, so this is a safe, local rename).

3. **`loadnav.m`**: a real, pre-existing bug, not introduced by us. Its own
   `setdefv()` call (line 72) defaults `nav_time_base` onto the `p` struct,
   but the `switch` that consumes it (line 137) reads `f.nav_time_base` --
   with only `p.nav_time_base` set, this dies with "structure has no member
   'nav_time_base'". `loadctd.m`'s changelog (line 85) shows the *identical*
   `ctd_time_base` bug was fixed on 2014-03-21 (default moved from `p` to
   `f`); the fix was apparently never ported to `loadnav.m`. **Not
   patched** -- worked around instead by setting `f.nav_time_base` directly
   in `set_cast_params.m` (see below), since editing call/default
   consistency inside `loadnav.m` felt like more surgery than the harness
   warranted for a single missing default.

## `end_processing_step.m` -- harness instrumentation

- Added the M2 per-step dump (`save('-v6', 'dumps/stepNN.mat', 'd','p')`,
  plus `di`/`dr` when they exist, wrapped in try/catch).
- Fixed a real save-vs-load mismatch while adding the above: the existing
  checkpoint save (`eval(sprintf('save %s_%d', f.checkpoints, pcs.cur_step))`)
  has no `.mat` extension. MATLAB's `save name` auto-appends `.mat`; Octave's
  does not. `begin_processing_step.m`'s `load(sprintf('%s_%d.mat', ...))`
  expects the extension, so checkpoint resume silently failed (checkpoint
  file "003_1" written, "003_1.mat" never found) until this was fixed to
  `save %s_%d.mat`. This was essential for iterating on the plotting-stub
  fixes below without re-running `loadrdi.m` (~4 min) every time.

## Missing/incompatible functions -- fixed via `octave_harness/stubs/`

`docs/legacy/` and the LDEO_IX_Software.tar do not contain `makebars.m`
anywhere -- genuinely missing (not an Octave/MATLAB API difference), used
by `plotraw.m` for a diagnostic bar overlay. Stubbed with harmless
placeholder outputs (see stub file).

`interp1q` is a MATLAB builtin not implemented in Octave (confirmed against
gnuoctave/octave:9.2.0) -- stubbed as a thin wrapper over `interp1(...,'linear')`.

All other stubs are pure no-op plotting functions (CONTINUATION_PLAN.md
explicitly anticipated this fallback): `figure`, `plot`, `subplot`, `hold`,
`axis`, `title`, `xlabel`, `ylabel`, `text`, `legend`, `colorbar`, `pcolor`,
`contourf`, `streamer`, `orient`, `print`, `pause`, `clf`, `gca`, `grid`,
`colormap`, `imagesc`, `shading`, `fill`, `caxis`, `bar`, `axes`. Each was
added lazily, one crash at a time, per the plan's guidance ("don't
pre-write fifty"). `set.m` is the one exception with real logic: it
no-ops only when called on a numeric (graphics-handle) first argument,
falling through to `builtin('set', ...)` otherwise, per the plan's
explicit caution not to blanket-stub `set`/`get`.

None of these touch `src/ladcp/` or change any numerical result -- they
only silence a diagnostic-plotting subsystem we don't need for the
stage-by-stage numeric comparison.

## `octave_harness/ldeo_ix/set_cast_params.m` -- new file (harness-specific)

Not a patch to an LDEO file -- `set_cast_params.m` is a per-cruise/per-cast
script every LDEO_IX deployment must supply (`process_cast.m` calls it and
errors if it's missing). Built entirely from the recorded p-struct in
`test_data/2015_P16N/003.nc` global attributes; see the file's own comments
for which fields are genuine inputs vs. computed outputs we deliberately
left unset.
