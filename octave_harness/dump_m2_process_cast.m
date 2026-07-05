% M2 -- run process_cast(3) through step 14, dumping d/p/di/dr after each
% end_processing_step call. See CONTINUATION_PLAN.md milestone M2.
addpath('/work/octave_harness/ldeo_ix');
addpath('/work/octave_harness/stubs');
cd('/work/octave_harness/work');

mkdir('dumps');  % per-step dumps land here (octave_harness/work/dumps/, gitignored)

stn = 3;
% stop=0 ("don't stop"): stop=2 ("stop after all steps") looked like the
% right choice from the docstring, but end_processing_step.m only resets
% pcs.stop to 0 when it equals 1 -- with stop=2 it never resets and
% `keyboard` (Octave's interactive debug prompt) fires after every single
% step, hanging non-interactively. stop=0 runs straight through.
%
% BEGIN_STEP is overridable via the OCTAVE_HARNESS_BEGIN_STEP env var so
% stub-fixing iterations can resume from the last good checkpoint instead
% of re-running loadrdi.m (~4 min) every time. Defaults to 1 (full run).
begin_step_str = getenv('OCTAVE_HARNESS_BEGIN_STEP');
if isempty(begin_step_str)
  begin_step = 1;
else
  begin_step = str2num(begin_step_str);
end
process_cast(stn, begin_step, 0);

disp('M2 process_cast run complete');
