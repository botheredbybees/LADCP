function d = edit_data(d,p)
% function d = edit_data(d,p)
%
% perform data editing (e.g. sidelobes, previous-ping interference, &c)
%
% HISTORY:
%  Jul  3, 2004: - implemented side-lobe contamination editing
%  Jul 15, 2004: - implemented PPI contamination editing
%  Jul 16, 2004: - implemented cross-contamination editing
%  Jul 18, 2004: - added parameter defaults
%  Jul 21, 2004: - moved bin masking from [loadrdi.m]
%  Jul 23, 2004: - implemented ensemble skipping
%  May 18, 2006: - incorporated simple fix for multi-ping ensembles
%		   provided by Mattew Alford
%		 - disabled edit_spike_filter by default
%		 - changed edit_spike_filter_max_curv default value
%  Aug 13, 2013: - BUG: edit_sidelobes did not work for UL only data
%  Jul 13, 2014: - automatically edit bins 1 when blanking distance is zero
%  Jan 23, 2015: - BUG: automatic zero blanking editing had a typo
%		 - BUG: automatic zero blanking editing did not work with DL-only data
%		 - added p.edit_{dn,up}_bad_ensembles
%  Jul 21, 2015: - made bin-masking more permissive (allow indices > #bins)

%======================================================================
%                    E D I T _ D A T A . M 
%                    doc: Sat Jul  3 17:13:05 2004
%                    dlm: Tue Jul 21 11:39:18 2015
%                    (c) 2004 A.M. Thurnherr
%                    uE-Info: 126 23 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

%----------------------------------------------------------------------
% Initialize & Set Default Parameters
%----------------------------------------------------------------------

d.ts_edited = d.ts; % save for plotting/inspection

% Set to list of bins to remove from data.
p=setdefv(p,'edit_mask_dn_bins',[]);
p=setdefv(p,'edit_mask_up_bins',[]);

% Set to 1 to remove side-lobe contaminated data near seabed and
% surface.
p=setdefv(p,'edit_sidelobes',1);

% Set to 1 to implement time-domain spike filter on the data; 
% this removes interference from other acoustic instruments but,
% more importantly, can get rid of PPI when staggered pings
% are used. USING THIS FILTER WITH A CURVATURE CRITERION THAT IS
% TOO STRICT CAN SERIOUSLY DEGRADE THE VELOCITY PROFILES. THEREFORE,
% IT HAS BEEN DISABLED BY DEFAULT IN VERSION IX_3.
p=setdefv(p,'edit_spike_filter',0);

% Spike filtering is done using 2nd-difference
% peak detection in time. This parameter gives the maximum
% target-strength 2nd derivative that's allowed. Set to larger
% values to weaken the filtering. (Check figure 14 to see if
% filter is too strong or too weak.) Setting the value of this
% parameter too low will seriously degrade the velocity profiles,
% without any apparent sign of trouble. The optimal value of this
% parameter depends on the instrument type.
p=setdefv(p,'edit_spike_filter_max_curv',5);

% Set to 1 to remove data contaminated by previous-ping interference.
% NOTES: - Using the spike filter seems to work more robustly, as long
%     	   as staggered pings are used. Great care must be taken,
%	   however, to chose a good value for edit_spike_filter_max_curv.
p=setdefv(p,'edit_PPI',0);

% PPI layer thickness in meters; the value is taken directly from Eric
% Firing's default (2*clip_margin = 180m).
p=setdefv(p,'edit_PPI_layer_thickness',180);

% max distance from seabed at which PPI should be removed. This is
% an observed parameter and depends on the clarity of the water.
% Check Figure 14 to see whether this should be changed.
p=setdefv(p,'edit_PPI_max_hab',1000);

% set this vector to enable skipping of ensembles; skipping vector
% is wrapped around, i.e. [1 0] skips all odd ensembles, [0 1 0] skips
% ensembles 2 5 8 11.... This filter is useful to process the casts
% with only half the data to see whether the two halves agree, which
% probably means that the cast can be trusted. Note that if staggered
% ping setup is used to avoid PPI, the skipping vector should leave
% adjacent ensembles intact, i.e. use something like [1 1 0 0] and
% [0 0 1 1].
p=setdefv(p,'edit_skip_ensembles',[]);

% the following vectors can be used to edit out blocks of bad ensembles,
% caused, for example, by intermittent hardware failures. NOTE: it is
% assumed that the ensembles are numbered consecutively, starting with 1,
% i.e. the ensemble numbers in the data files are ignored. This may 
% not work as intended if the data files are trimmed during the DL/UL
% merge in [loardrdi.m]

p=setdefv(p,'edit_dn_bad_ensembles',[]);
p=setdefv(p,'edit_up_bad_ensembles',[]);

%----------------------------------------------------------------------
% Bin Masking
%----------------------------------------------------------------------

if length(d.zu)>0 && p.blnk_u==0
  disp(sprintf(' bin masking               : masking uplooker bin 1 because of zero blanking distance'));
  p.edit_mask_up_bins = [1 p.edit_mask_up_bins];
end
if p.blnk_d==0
  disp(sprintf(' bin masking               : masking downlooker bin 1 because of zero blanking distance'));
  p.edit_mask_dn_bins = [1 p.edit_mask_dn_bins];
end

if ~isempty(p.edit_mask_dn_bins) | ~isempty(p.edit_mask_up_bins)

  nbad = 0;
  if length(d.zu) > 0
    for bi=1:length(p.edit_mask_up_bins)
      if p.edit_mask_up_bins(bi)<=p.nbin_u
         bn = length(d.zu)+1 - p.edit_mask_up_bins(bi);
	 nbad = nbad + length(find(isfinite(d.weight(bn,:))));
	 d.weight(bn,:) = NaN; d.ts_edited(bn,:) = NaN;
      end
    end
  end
  if length(d.zd) > 0
    for bi=1:length(p.edit_mask_dn_bins)
      if p.edit_mask_dn_bins(bi)<=p.nbin_d
        bn = length(d.zu) + p.edit_mask_dn_bins(bi);
	nbad = nbad + length(find(isfinite(d.weight(bn,:))));
	d.weight(bn,:) = NaN; d.ts_edited(bn,:) = NaN;
      end
    end
  end

  disp(sprintf(' bin masking               : set %d weights to NaN',nbad));

end % if bin masking enabled

%----------------------------------------------------------------------
% Side-Lobe Contamination
%----------------------------------------------------------------------

if p.edit_sidelobes

  nbad = 0;
  
  % first, the uplooker: d.z is -ve distance of ADCP from surface;
  % Cell_length is in cm, i.e. 0.015*Cell_length is 1.5 x bin size
  % in m --- the same value used by Eric Firing's software
  
  if length(d.zu)==0 && d.zd(1)<0			% UL only (in DL structures)

    for b=1:length(d.zd)
      zlim(b,:) = (1 - cos(pi*d.down.Beam_angle/180)) * d.z ...
		- 0.015*d.down.Cell_length;
    end
    ibad = find(d.izm > zlim);
    nbad = nbad + length(find(isfinite(d.weight(ibad))));
    d.weight(ibad) = NaN; d.ts_edited(ibad) = NaN;
  
  elseif length(d.zu > 0)				% DL/UL combo
  
    for b=1:length(d.zu)+length(d.zd)
      zlim(b,:) = (1 - cos(pi*d.up.Beam_angle/180)) * d.z ...
		- 0.015*d.up.Cell_length;
    end
    ibad = find(d.izm > zlim);
    nbad = nbad + length(find(isfinite(d.weight(ibad))));
    d.weight(ibad) = NaN; d.ts_edited(ibad) = NaN;
  
  end
  
  % now, the downlooker: p.zbottom is the +ve depth of the sea bed; therefore,
  % -d.z - p.zbottom is the -ve distance from the sea bed
  
  for b=1:length(d.zu)+length(d.zd)
    zlim(b,:) = -p.zbottom ...
	      + (1 - cos(pi*d.down.Beam_angle/180)) * (d.z+p.zbottom) ...
	      + 0.015*d.down.Cell_length;
  end
  ibad = find(d.izm < zlim);
  nbad = nbad + length(find(isfinite(d.weight(ibad))));
  d.weight(ibad) = NaN; d.ts_edited(ibad) = NaN;
  
  disp(sprintf(' side-lobe contamination   : set %d weights to NaN',nbad));

end %if p.edit_sidelobes

%----------------------------------------------------------------------
% Time-Domain Spike Filter
%----------------------------------------------------------------------

if p.edit_spike_filter

  nbad = 0;
  for b=1:length(d.zd)+length(d.zu)
    ibad = find(diff(diff(d.ts(b,:))) < -1*p.edit_spike_filter_max_curv) + 1;
    nbad = nbad + length(find(isfinite(d.weight(b,ibad))));
    d.weight(b,ibad) = NaN; d.ts_edited(b,ibad) = NaN;
  end
  disp(sprintf(' spike filter              : set %d weights to NaN',nbad));

end %if p.edit_spike_filter

%----------------------------------------------------------------------
% Previous-Ping Interference
%----------------------------------------------------------------------

if p.edit_PPI

  % NB: at present, PPI filtering is only implemented for the downlooker
  
  nbad = 0;
  
  % calc ping-intervals; dt(1) contains the difference (in seconds)
  % between the first two ensembles (t(2) - t(1)).
  
  dt = diff(d.time_jul)*86400;

  % adjust for multi-ping ensembles; THIS FIX, PROVIDED BY MATTHEW ALFORD,
  % DOES NOT WORK FOR IRREGULAR PINGING RATES (e.g. 10s ensembles of 3
  % pings each with 1s between pings)

  dt = dt/d.down.Pings_per_Ensemble;
  
  % use the mean sound speed below the approximate expected PPI depth;
  % this is anal but not very expensive; using 1500m/s would be nearly
  % as good.
  
  if isfield(d,'ctdprof_z') & isfield(d,'ctdprof_ss')
    guess_z = -p.zbottom + 1500 * meannan(dt) / 2;
    SS = meannan(d.ctdprof_ss(find(d.ctdprof_z > -guess_z)));
    if ~isfinite(SS), SS = 1500; end
  else
    SS = 1500;
  end
  
  % calculate the depth limits to remove the PPI for all (but the first)
  % ensembles; the beam-angle limits were found to be too conservative
  % and were replaced by a nominal layer_tickness.
  
  %PPI_min_beam_angle = 0.0 * d.down.Beam_angle;
  %PPI_max_beam_angle = 1.2 * d.down.Beam_angle;
  %PPI_max_z = -p.zbottom + SS * dt/2 * cos(pi/180 * PPI_min_beam_angle);
  %PPI_min_z = -p.zbottom + SS * dt/2 * cos(pi/180 * PPI_max_beam_angle);
  
  PPI_hab = SS * dt/2 * cos(pi/180 * d.down.Beam_angle);
  PPI_hab(find(PPI_hab > p.edit_PPI_max_hab)) = inf;
  
  PPI_min_z = -p.zbottom + PPI_hab - p.edit_PPI_layer_thickness / 2;
  PPI_max_z = PPI_min_z + p.edit_PPI_layer_thickness;
  
  % remove the contaminated data from the downlooker bins
  
  for b=length(d.zu)+1:length(d.zu)+length(d.zd)
    ibad = find(d.izm(b,2:end) > PPI_min_z & d.izm(b,2:end) < PPI_max_z) + 1;
    nbad = nbad + length(find(isfinite(d.weight(b,ibad))));
    d.weight(b,ibad) = NaN; d.ts_edited(b,ibad) = NaN;
  end
  
  disp(sprintf(' previous-ping interference: set %d weights to NaN',nbad));

end %if p.edit_PPI

%----------------------------------------------------------------------
% Ensemble Skipping
%----------------------------------------------------------------------

if ~isempty(p.edit_skip_ensembles)

  nskipped = 0;
  iskip = [];
  for i=1:length(p.edit_skip_ensembles)
    if p.edit_skip_ensembles(i)
      iskip = [iskip i:length(p.edit_skip_ensembles):length(d.time_jul)];
    end
  end

  for b=1:length(d.zd)+length(d.zu)
    nskipped = nskipped + length(find(isfinite(d.weight(b,iskip))));
    d.weight(b,iskip) = NaN; d.ts_edited(b,iskip) = NaN;
  end

  disp(sprintf(' ensemble skipping         : set %d weights to NaN',nskipped));

end % if p.edit_skip_ensembles enabled

%----------------------------------------------------------------------
% remove blocks of bad ensembles from UL and/or DL data
%----------------------------------------------------------------------

if ~isempty(p.edit_dn_bad_ensembles) || ~isempty(p.edit_up_bad_ensembles)
  if ~isempty(p.edit_dn_bad_ensembles)
    dn_bad = 0; ibad = p.edit_dn_bad_ensembles;
    for b=length(d.zu)+1:length(d.zu)+length(d.zd)
      dn_bad = dn_bad + length(find(isfinite(d.weight(b,ibad))));
      d.weight(b,ibad) = NaN; d.ts_edited(b,ibad) = NaN;
    end
  end
  if ~isempty(p.edit_up_bad_ensembles)
    up_bad = 0; ibad = p.edit_up_bad_ensembles;
    for b=length(d.zu):-1:1
      up_bad = up_bad + length(find(isfinite(d.weight(b,ibad))));
      d.weight(b,ibad) = NaN; d.ts_edited(b,ibad) = NaN;
    end
  end
  disp(sprintf(' DL/UL ensemble editing    : set %d/%d weights to NaN',dn_bad,up_bad));
end

%----------------------------------------------------------------------
% Plot Results of Editing
%----------------------------------------------------------------------

bin_no = [0];
if length(d.zu) > 0, bin_no = [-length(d.zu):1 bin_no]; end
if length(d.zd) > 0, bin_no = [bin_no 1:length(d.zd)]; end

figure(14);
clf;
orient landscape;
colormap([[1 1 1]; jet(128)]);

subplot(2,1,1);
imagesc(1:size(d.ts,2),bin_no,...
	[d.ts(1:length(d.zu),:); ...
	 ones(1,size(d.ts,2))*NaN; ...
	 d.ts(size(d.ts,1)-length(d.zd)+1:end,:)...
        ]);
csc = caxis;
xlabel('Ensemble #');
ylabel('Bin #');
title('Before Data Editing');

subplot(2,1,2);
imagesc(1:size(d.ts,2),bin_no,...
	[d.ts_edited(1:length(d.zu),:); ...
	 ones(1,size(d.ts,2))*NaN; ...
	 d.ts_edited(size(d.ts,1)-length(d.zd)+1:end,:)...
        ]);
csc = caxis;
xlabel('Ensemble #');
ylabel('Bin #');
title('After Data Editing');

streamer([p.name,'  Figure 14']);

