% M1 -- run loadrdi.m alone on P16N cast 003 DL/UL PD0 files, dump d/p structs.
% See CONTINUATION_PLAN.md milestone M1.
addpath('/work/octave_harness/ldeo_ix');

f = struct();
f.ladcpdo = '/work/test_data/2015_P16N/003DL000.000';
f.ladcpup = '/work/test_data/2015_P16N/003UL000.000';

p = struct();
p.name = '003';
% Recorded magnetic deviation from 003.nc GEN_Magnetic_deviation_deg -- avoids
% depending on magdev.m (lat/lon/date) for this ingestion-only milestone.
p.drot = 12.318441;
% default.m normally sets this before process_cast calls loadrdi; loadrdi.m's
% own setdefv() calls don't cover it (outlier.m reads p.debug directly).
p.debug = 0;

[d, p, de] = loadrdi(f, p);

mkdir('/work/octave_harness/dumps');
save('-v6', '/work/octave_harness/dumps/m1_loadrdi.mat', 'd', 'p');
disp('M1 dump complete');
