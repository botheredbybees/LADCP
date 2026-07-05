%======================================================================
%                    L O A D C T D . M 
%                    doc: Sat Jun 26 15:56:43 2004
%                    dlm: Fri Sep 14 18:15:35 2018
%                    (c) 2004 M. Visbeck & A. Thurnherr
%                    uE-Info: 545 36 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

function [d,p]=loadctd(f,d,p)
% function [d,p]=loadctd(f,d,p)

% This routine works for generic ASCII files, containing CTD time series,
% with fields pressure, temperature, salinity, and time.

%====================
% TWEAKABLES
%====================

% PRE-AVERAGING OF GPS TIME (IN DAYS)
p=setdefv(p,'navtime_av',2/60/24);

% set ctdmaxlag.=100 to the maximum pings that the ADCP data can be shifted to
% best match W calculated from CTD pressure time series (loadctd)
% If you have good times set it to 10... if your time base is questionable
% you can set it 100 or more
p=setdefv(p,'ctdmaxlag',150);

% NUMBER OF DATA CHUNKS TO USE FOR TIME LAGGING
p=setdefv(p,'ctdmaxlagnp',600);

% restrict time range to profile and disregard data close to surface
% p.cut = 0 dont restrict
% p.cut > 0 restrict time to adcp depth below a depth of p.cut
p=setdefv(p,'cut',10);

% INTERPOLATE IRREGULAR CTD TIME SERIES
p=setdefv(p,'interp_ctd_times',1);

% FILE LAYOUT
f = setdefv(f,'ctd_header_lines',0);
f = setdefv(f,'ctd_fields_per_line',4);
f = setdefv(f,'ctd_time_field',1);
f = setdefv(f,'ctd_pressure_field',2);
f = setdefv(f,'ctd_temperature_field',3);
f = setdefv(f,'ctd_salinity_field',4);
f = setdefv(f,'ctd_badvals',-9e99);

% TIME BASE
% 	0 for elapsed time in seconds
% 	1 for year-day (1.0 = Jan 1, 00:00)
% 	2 for Visbeck's Gregorian (see gregoria.m)
f = setdefv(f,'ctd_time_base',0);

%======================================================================

% MODIFICATIONS BY ANT:
%   Jun 26, 2004: - totally re-wrote the file-reading, which reduces
%		    run time by factor 6 or so for deep casts
%   Jun 29, 2004: - estimate time offset at 90% max depth during downcast
%		    instead of at max depth (which is not well defined
%		    if there is a bottle stop)
%   Jun 30, 2004: - added handling for different time bases
%   		  - added optional interpolation code
%   Jul  1, 2004: - removed ipos argument
%   Dec 18, 2004: - BUG: header_lines > 0 did not work
%   Jul 19, 2006: - removed NaN before call to interp1 to avoid Matlab 7.2 warning
%   Jan  5, 2007: - merged from 2 different old versions
%		  - added CTD file layout into p structure
%   Jan 26, 2007: - BUG: file layout default was in p (rather than f) structure
%   May 23, 2008: - added error message on imaginary soundspeed values
%		  - added f.ctd_badvals
%   Jun 26, 2008: - change default of p.interp_ctd_times to 1, because this
%		    improves time matching of irregularly spaced CTD time
%		    series without affecting regularly spaced ones
%   Jun 30, 2008: - BUG: typo in log output
%   Jul  2, 2008: - removed "not recomended" (sic) comment from log file when
%		    profile is created from time series
%   Jul 17, 2008: - moved some code from [loadnav.m] to fix bug associated with
%		    adjusting of start/end positions in case of significant
%		    ADCP vs GPS/CTD clock offset (unclear whether this happened
%		    only when elapsed time was used)
%   Sep 18, 2008: - BUG: moved code had assumed that nav data were loaded
%   Jan  7, 2009: - tightened use of exist()
%   Jun 16, 2009: - BUG: patching short nav time series did not work correctly
%   Mar 21, 2014: - BUG: f.ctd_time_base used p.ctd_time_base set as default
%   May 27, 2015: - removed confusing diagnostic message regarding adjusting NAV time
%   May 28, 2015: - added error message when there are no valid vertical velocities
%   Apr 18, 2018: - BUG: ADCP-time shift warning was meaningless with elapwsed time_base
%   Sep 14, 2018: - BUG: code move in 2008 broke working with single GPS file for entire cruise

% read SEABIRD ctd timeseries file
disp(['LOADCTD: load CTD time series ',f.ctd])
if ~exist(f.ctd,'file')
 warn=([' LOADCTD can not find ',f.ctd]);
 p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
 disp(warn)
 return
end

% construct input format
cur_field = 1; input_format = '';
for i=1:f.ctd_fields_per_line
  switch i,
    case f.ctd_pressure_field,
      i_press = cur_field; cur_field = cur_field + 1;
      input_format = [input_format ' %g'];
    case f.ctd_temperature_field
      i_temp = cur_field; cur_field = cur_field + 1;
      input_format = [input_format ' %g'];
    case f.ctd_salinity_field
      i_salin = cur_field; cur_field = cur_field + 1;
      input_format = [input_format ' %g'];
    case f.ctd_time_field
      i_time = cur_field; cur_field = cur_field + 1;
      input_format = [input_format ' %g'];
    otherwise
      input_format = [input_format ' %*g'];
  end
end
if cur_field ~= 5
  error('File format definition error');
end

% open input & skip header
fp=fopen(f.ctd);
header_lines = f.ctd_header_lines;
while header_lines > 0
  header_lines = header_lines - 1;
  fgets(fp);
end

% read time series
[A,nread] = fscanf(fp,input_format,[4,inf]);
fclose(fp);

% remove bad values
ibad = find(A == f.ctd_badvals);
if (length(ibad) > 0)
  A(ibad) = nan;
  disp(sprintf(' Warning: %d bad values (%g) removed from CTD data',length(ibad),f.ctd_badvals));
end

% CTD time
timctd=A(i_time,:)';
switch f.ctd_time_base
  case 0 % elapsed time in seconds
    timctd = timctd/24/3600 + julian(p.time_start);
  case 1 % year-day
    timctd = timctd + julian([p.time_start(1) 1 0 0 0 0]);
end

disp(sprintf(' read %d CTD scans; median delta_t = %.2f seconds',...
		length(timctd),median(diff(timctd))*24*3600));

% pressure, temperature, salinity
data=A([i_press i_temp i_salin],:)';

% interpolate to regular time series
if p.interp_ctd_times
  min_t = min(timctd);
  max_t = max(timctd);
  delta_t = median(diff(timctd));
  data = interp1q(timctd,A([i_press i_temp i_salin],:)',[min_t:delta_t:max_t]');
  timctd = [min_t:delta_t:max_t]';
  disp(sprintf(' interpolated to %d CTD scans; delta_t = %.2f seconds',...
  		length(timctd),median(diff(timctd))*24*3600));
end

% calc LADCP depth
w=meannan(d.rw);
if sum(isfinite(w)) == 0
    error('No valid vertical velocities --- aborting');
end
ii=find(~isfinite(w));
w(ii)=0;
dt=diff(d.time_jul)*24*3600;
dt=mean([0,dt;dt,0]);
z=cumsum(w.*dt);

% guess time offset; in Martin's version, this was done based on bottom
% time, which does not work well when there is a bottle stop at the bottom, as
% the LADCP depth can drift slowly. Now, the guess is done based on the time
% when the package reaches 90% of the maximum depth/pressure.
[zmax,izmax]=max(z);
[pmax,ipmax]=max(data(:,1));
i_ctd_near_bottom   = min(find(data(:,1) > 0.9*pmax));
i_ladcp_near_bottom = min(find(z         > 0.9*zmax));
delta_t = (d.time_jul(i_ladcp_near_bottom) - timctd(i_ctd_near_bottom))*24*60;

disp([' 90% LADCP depth  ',int2str(0.9*zmax),' at ',datestrj(d.time_jul(i_ladcp_near_bottom))])
disp([' 90% CTD pressure ',int2str(0.9*pmax),' at ',datestrj(timctd(i_ctd_near_bottom))])

% in case elapsed time is given, estimate absolute time
if f.ctd_time_base == 0
  p.ctdtimoff = delta_t/60/24;
  timctd = timctd + p.ctdtimoff;
  delta_t = 0;
end

% now, handle warnings; NB: cannot happen if elapsed time is given
if abs(delta_t)>5
  disp('WARNING WARNING ');
  warn=(sprintf(' estimated time offset between CTD/LADCP is %.1f minutes',delta_t));
  disp(warn); p.warn(size(p.warn,1)+1,1:length(warn))=warn;
end
if abs(delta_t)>60
  p.ctdtimoff = round(delta_t/60)/24;
  p.ctdtimoff = delta_t/60/24;
  disp(sprintf(' maybe the hour was set wrong: adjusting CTD time by %d days',...
 	         p.ctdtimoff));
  timctd = timctd + p.ctdtimoff;
end

if p.navdata && f.nav_time_base == 0
  if ~strcmp(f.ctd,f.nav)
    disp(sprintf(' adjusting GPS time to CTD time (%+d seconds)',floor(p.ctdtimoff*24*3600)))
  end
  d.navtime_jul = d.navtime_jul + p.ctdtimoff;
end

% plot profile
figure(4);
clf;
orient tall;

subplot(325)
tim0=timctd(1);
plot((timctd-tim0)*24,-data(:,1))
axis tight
title('Cut CTD profile')
xlabel('Time in hours')
ylabel('CTD pressure [dBar]')
hold on
ax=axis;
%   start time
if p.cut~=0
 i10=find(abs(data(1:ipmax,1)-p.cut/2)<p.cut/2);
 if length(i10)>0
  istart=i10(end);
  p=setdefv(p,'ctdtime',1);
  plot(24*(timctd(i10)-tim0),-data(i10,1),'.g')
 else
  istart=1;
  if abs(data(1,1)-p.cut)<10
   p.ctdtime=1;
  else
   disp([' WARNING first CTD depth is: ',int2str(data(1,1)),' [dB]'])
   disp([' WARNING ignore pressure time series '])
   p.ctddepth=0;
   p.ctdtime=0;
   return
  end
 end
else
 istart=1;
 p.ctdtime=0;
end

p.ctd_starttime=timctd(istart);
plot((timctd(istart)-tim0)*24,-data(istart,1),'or')
plot([1 1]*24*(timctd(istart)-tim0),ax(3:4),'r--')

%   end time
if p.cut~=0
 i10=find(abs(data(ipmax:end,1)-p.cut/2)<p.cut/2)+ipmax-1;
 if length(i10)>0
  iend=i10(1);
  plot((timctd(i10)-tim0)*24,-data(i10,1),'.g')
 else
  iend=length(timctd);
  p.ctdtime=0;
 end
else
 iend=length(timctd);
end
p.ctd_endtime=timctd(iend);
plot((timctd(iend)-tim0)*24,-data(iend,1),'or')
plot([1 1]*24*(timctd(iend)-tim0),ax(3:4),'r--')
pause(0.01)

if p.ctdtime==1
 disp([' Changed Start Time : ',datestrj(d.time_jul(1)),...
      '  to ',datestrj(p.ctd_starttime)])
 p.time_start=gregoria(p.ctd_starttime);
 disp([' Changed End   Time : ',datestrj(d.time_jul(end)),...
      '  to ',datestrj(p.ctd_endtime)])
 p.time_end=gregoria(p.ctd_endtime);
end

% ===================================================================
% at this point you NEED to have read the raw CTD data
%
% timectd to be in the same time units as the ADCP data
% data a matrix starting with pressure (not depth) in dbar
%
% if you have synchronized CTD ADCP data jump below
% =================================================================

% fix time for interpolation
if min(timctd)>max(d.time_jul) | max(timctd)<min(d.time_jul)
  disp('CTD timeseries does not overlap WRONG STATION????')
  disp(['  CTD-data : ',datestrj(timctd(1)),'  to ',datestrj(timctd(end))])
  disp([' ADCP-data : ',datestrj(d.time_jul(1)),'  to ',datestrj(d.time_jul(end))])
  disp(' will try to use W-data to get depth')
  p.ctddepth=0;
  return
end

p.ctddepth=1;

% check for latitude
if existf(p,'poss')
  lat = p.poss(1);
elseif p.navdata
  lat = medianan(d.slat);
else
  disp('WARNING WARNING ');
  warn=(' do not have required latitude for accurate pressure->depth conversion; using 0 instead');
  disp(warn); p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  lat = 0;
end

ii=find(diff(timctd)==0);
data(ii,:)=[];
timctd(ii)=[];
[ln,lv]=size(data);

% check for spikes in pressure data
%
dtadcp=medianan(diff(d.time_jul));
dtctd=medianan(diff(timctd),length(timctd)/5);
disp([' median CTD time difference ',num2str(dtctd*24*3600),' s'])
% fix jitter in CTD time base
timctd2=cumsum(round((diff(timctd))/dtctd)*dtctd);
timctd2=[timctd(1); timctd(1)+timctd2];


ibad=1:10;
sbad=0;
nn=0;
while length(ibad)>9 & nn<10
 z=-p2z(data(:,1),lat);
 wctd=[-diff(z)./(diff(timctd2)*24*3600)];
 ibad=find(abs(wctd)>(3*std(wctd)));
 ibad=[ibad; ibad+1];
 data(ibad,:)=[];
 timctd(ibad,:)=[];
 timctd2(ibad,:)=[];
 sbad=sbad+length(ibad);
 nn=nn+1;
end
disp([' removed ',int2str(sbad),' pressure spikes'])
z=-p2z(data(:,1),lat);

% compute wctd similar to what ADCP would do
nshift=max([1,round(dtadcp/dtctd)])+1;
nshift2=fix(nshift/2);
disp([' use ',int2str(nshift2),' time base for W_ctd'])
wctd=z+nan;
i2=(nshift):(length(z)-nshift2);
i1=1:length(i2);

% get W from CTD 
wctd(i1+nshift2)=-(z(i2)-z(i1))./((timctd2(i2)-timctd2(i1))*24*3600);

if sbad>10
  warn=[' removed ',int2str(sbad),' pressure spikes during: ',int2str(nn),' scans'];
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
end

% interpolate on ADCP time series
timctd(1)=0;
timctd(end)=1e32;
data_int=interp1(timctd,data,d.time_jul','nearest');
igood = find(isfinite(wctd));
d.wctd=interp1(timctd(igood),wctd(igood),d.time_jul','nearest')';
d.z=-p2z(data_int(:,1),lat)';
disp([' CTD max depth : ',int2str(-min(d.z))])

if p.ctdmaxlag>0
 % prepare for time lag check
 dt=mean(diff(d.time_jul));
 w=medianan(d.rw,2);

 iw=sum(isfinite(d.rw));
 ii=find(iw<4);
 w(ii)=nan;

 % check for timelag
 dtctd=medianan(diff(timctd));

 % make up array to check for lag
 ctdtw=[timctd, wctd];
 adcptw=[d.time_jul', w'];
 
 [lag,co]=besttlag(ctdtw,adcptw,p.ctdmaxlag,p.ctdmaxlagnp);
 lagdt=-lag*dtctd;
 disp([' best lag W: ',int2str(lag),' CTD scans ~ ',...
             int2str(lagdt*24*3600),' seconds  corr:',num2str(co)]);
 p=setdefv(p,'ctdmincorr',max(0.9-100/length(wctd),0.6));
 if co<p.ctdmincorr
 % try whole cast
  [lag,i1,i2,co]=bestlag(d.wctd,w,p.ctdmaxlag);
  lagdt=-lag*dt;
  disp([' best lag W: ',int2str(lag),' CTD scans ~ ',...
             int2str(lagdt*24*3600),' seconds  corr:',num2str(co)]);
 end
 % recompute W-ctd
 data_int=interp1(timctd,data,d.time_jul'+lagdt,'linear');
 igood = find(isfinite(wctd));
 d.wctd=interp1(timctd(igood),wctd(igood),d.time_jul'+lagdt,'nearest')';
 d.z=-p2z(data_int(:,1),lat)';
 i1=1:length(w);

% plot two small sections for visual control
 ii=[-12:12]+fix(length(i1)*0.7);
 ii=[ii, [-12:12]+fix(length(i1)*0.3)];
 ii=ii(find(ii>0 & ii<length(i1)));

 subplot(326)
 plot(w(i1(ii)),'-b')
 hold on
 plot(d.wctd(i1(ii)),'-r')
 axis tight
 ax=axis;
 plot([25 25],ax(3:4),'-k')
 text(5,mean(w(i1(ii(25:end)))),'down cast')
 text(30,mean(w(i1(ii(1:25)))),'up cast')
 title(['best lag W: ',int2str(lag),' scans ~ ',...
               int2str(lagdt*24*3600),' sec.  c:',num2str(co)]);
 ylabel('W used for lag correlation')
 xlabel('sample scans (b)ADCP (r)CTD')
 grid
 pause(0.1)

 if p.debug==9, keyboard, elseif p.debug>1, pause, end


 % decide if the correlation is good enought to shift ADCP time
 if abs(lag)<p.ctdmaxlag & co>p.ctdmincorr
  disp(' adjust ADCP time to CTD time and shift depth record ')
  d.time_jul=d.time_jul+lagdt;
  if f.ctd_time_base~=0 && lagdt*24*3600>10
   disp('WARNING WARNING WARNING')
   warn=[' shifted ADCP timeseries by ',int2str(lagdt*24*3600),' seconds '];
   disp(warn)
   p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  end
 else
  error('Cannot determine time offset between CTD and LADCP time series --- aborting');

  % The following old code has been disabled for IX_10 because it is known not to work correctly
  % for combined NAV/CTD data, because the NAV time series *must* be matched to the LADCP data
  % for the GPS constraint to be applied.

  text(1,0,'LOW CORRELATION not used','color','r')
  disp('WARNING WARNING WARNING')
  warn=(' lag too big or correlation too weak: ignore CTD time series');
  disp(warn)
  p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  p.ctddepth=0;
  if abs(lag)<p.ctdmaxlag 
   disp(' will still use max ctd depth to constrain int W')
   p.zpar(2)=-min(d.z);
  end
  rmfield(d,'z');
  rmfield(d,'wctd');
  return
 end
end

% save other CTD data 
d.ctd_data=data_int(:,2:end)';

% pressure data
d.ctd_p=data_int(:,1)';

% if you have a temperatue time series provide here
d.ctd_temp=data_int(:,2)';

% if you have sound speed time series provide here
% d.ctd_ss=data_int(:,iss-1)';

% have salinity time series provide here
d.ctd_s=data_int(:,3)';

[pmax,imax]=max(d.ctd_p);

pp=[max(d.ctd_p(1),0):pmax]';
ipp=1:imax;
[dum,is]=sort(d.ctd_p(ipp));
ipps=ipp(is);

if existf(d,'ctdprof_p')==0
 % make ctd-profile from data;
 d.ctdprof_p=pp;
 d.ctdprof_lat=lat;
 d.ctdprof_z=interp1q(d.ctd_p(ipps)',-d.z(ipps)',pp);
 d.ctdprof_t=interp1q(d.ctd_p(ipps)',d.ctd_temp(ipps)',pp);
 d.ctdprof_s=interp1q(d.ctd_p(ipps)',d.ctd_s(ipps)',pp);

% get N^2
 if exist('sw_bfrq','file')
  d.ctdprof_N2=sw_bfrq(d.ctdprof_s,d.ctdprof_t,d.ctdprof_p,d.ctdprof_lat);
  d.ctdprof_N2(end+1)=d.ctdprof_N2(end);
  d.ctdprof_ss=sw_svel(d.ctdprof_s,d.ctdprof_t,d.ctdprof_p);
  if ~isreal(d.ctdprof_ss)
    error('sound-speed profile has imaginary values --- check CTD data');
  end
 else
  disp(' download SW routines to get more accurate soundspeed and N^2 ')
 end
 disp(' made CTD profile from time series data. ')
end

if p.ctdtime==1 		% clear all velocity data outside time range
 ii=(d.time_jul>=p.ctd_starttime) & (d.time_jul<=p.ctd_endtime);
 if length(ii)~=sum(ii)
  d=cutstruct(d,ii);
%%% CODE DISABLED BECAUSE IT IS UNCLEAR WHETHER IT IS CONSISTENT WITH
%%% REMAINDER OF CHANGES FOR VERSION IX_6
%%%  if ~existf(d,'slon')
%%%   % don't have time variable navigation data, adjust end position
%%%   % to account for shorter cast time
%%%   slat=p.poss(1)+p.poss(2)/60;
%%%   slon=p.poss(3)+p.poss(4)/60;
%%%   elat=p.pose(1)+p.pose(2)/60;
%%%   elon=p.pose(3)+p.pose(4)/60;
%%%   elon=slon+(elon-slon)*sum(ii)/length(ii);
%%%   elat=slat+(elat-slat)*sum(ii)/length(ii);
%%%   p.pose=[fix(elat), (elat-fix(elat))*60, fix(elon), (elon-fix(elon))*60];
%%%   disp(' adjusted end position to account for short CTD time series')
%%%  end
 end
end

%----------------------------------------------------------------------
% Code moved from [loadnav.m]:
%	- set start/end positions
%	- p.navtime_jul should be consistent with both timctd and d.time_jul
%	  at this stage
%----------------------------------------------------------------------

if p.navdata && f.nav_time_base == 0
  if min(d.navtime_jul)>max(d.time_jul) | max(d.navtime_jul)<min(d.time_jul)
    disp('NAV timeseries does not overlap WRONG STATION????')
    disp(['  NAV-data : ',datestrj(d.navtime_jul(1)),'	to ',datestrj(d.navtime_jul(end))])
    disp([' ADCP-data : ',datestrj(d.time_jul(1)),'  to ',datestrj(d.time_jul(end))])
    disp(' will ignore nav data')
    p.navdata = 0;
    d.navtime_jul = d.time_jul; % make sure same length
    d.slat = d.navtime_jul*NaN;
    d.slon = d.slat;
  else
    if d.navtime_jul(1) > d.time_jul(1)
      disp('NAV timeseries starts after ADCP timeseries: used first NAV value to patch ')
      disp(['  NAV-data : ',datestrj(d.navtime_jul(1)),'  to ',datestrj(d.navtime_jul(end))])
      disp([' ADCP-data : ',datestrj(d.time_jul(1)),'  to ',datestrj(d.time_jul(end))])
      d.navtime_jul(1) = d.time_jul(1);
    end
    
    if d.navtime_jul(end) < d.time_jul(end)
      disp('NAV timeseries ends before ADCP timeseries: used last NAV value to patch ')
      disp(['  NAV-data : ',datestrj(d.navtime_jul(1)),'  to ',datestrj(d.navtime_jul(end))])
      disp([' ADCP-data : ',datestrj(d.time_jul(1)),'  to ',datestrj(d.time_jul(end))])
      d.navtime_jul(end) = d.time_jul(end);
    end
  
    % find valid
    ii=find(diff(d.navtime_jul)>0);
    ii=[ii;length(d.navtime_jul)];
  
    % average over p.navtime_av days
    dt2m=[-p.navtime_av:(1/3600/24):p.navtime_av]';
    slon=medianan(interp1q(d.navtime_jul(ii),d.slon(ii),julian(p.time_start)+dt2m));
    slat=medianan(interp1q(d.navtime_jul(ii),d.slat(ii),julian(p.time_start)+dt2m));
    p.nav_start=[fix(slat), (slat-fix(slat))*60, fix(slon), (slon-fix(slon))*60];
    elon=medianan(interp1q(d.navtime_jul(ii),d.slon(ii),julian(p.time_end)+dt2m));
    elat=medianan(interp1q(d.navtime_jul(ii),d.slat(ii),julian(p.time_end)+dt2m));
    p.nav_end=[fix(elat), (elat-fix(elat))*60, fix(elon), (elon-fix(elon))*60];
  
    % interpolate on RDI data
    % this also shortens vectors to length of d.time_jul, which may be the
    % only thing that is really needed
    d.slon=interp1q(d.navtime_jul(ii),d.slon(ii),d.time_jul')';
    d.slat=interp1q(d.navtime_jul(ii),d.slat(ii),d.time_jul')';
  end

  p=setdefv(p,'poss',[NaN NaN NaN NaN]);
  [slat,slon] = pos2str(p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60);
  disp([' update start pos  from:',slat,'  ',slon])
  p.poss=p.nav_start;
  [slat,slon] = pos2str(p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60);
  disp(['			to:',slat,'  ',slon])
  
  p=setdefv(p,'pose',[NaN NaN NaN NaN]);
  [slat,slon] = pos2str(p.pose(1)+p.pose(2)/60,p.pose(3)+p.pose(4)/60);
  disp([' update end pos    from:',slat,'  ',slon])
  p.pose=p.nav_end;
  [slat,slon] = pos2str(p.pose(1)+p.pose(2)/60,p.pose(3)+p.pose(4)/60);
  disp(['			to:',slat,'  ',slon])

end % if p.navdata

%----------------------------------------------------------------------
% End of code moved from [loadnav.m]
%----------------------------------------------------------------------

p.zpar=[max(0,-d.z(1)), -min(d.z), max(0,-d.z(end))];

% ==================================================================
%
% This is what you NEED from the ctd data to make a difference
%
% d.z = depth (not pressure) from the CTD data interpolated to the
%        ADCP ping time 
%       NOTE that the depth is measure from the surface, i.e. only
%       negative values are to be expected
%
%
%====================================================================


function depth=p2z(p,lat)
% !!!!!! USES Z=0 AT P=0  (I.E. NOT 1ATM AT SEA SURFACE)
%	pressure to depth conversion using
%	saunders&fofonoff's method (deep sea res.,
%	1976,23,109-111)
%	formula refitted for alpha(p,t,s) = eos80
%	units:
%		depth         z        meter
%		pressure      p        dbars  (original in bars, but below
%                                              division by 10 is included)
%		latitude      lat      deg
%	checkvalue:
%		depth =       9712.654  m
%	for
%		p     =         1000.     bars
%		lat   =           30.     deg
%	real lat,p
        if nargin < 2, lat=54; end
        p=p/10.;
	x=sin(lat/57.29578);
	x=x*x;
	gr=9.780318*(1.0+(5.2788e-3+2.36e-5*x)*x)+1.092e-5*p;
	depth=(((-1.82e-11*p+2.279e-7).*p-2.2512e-3).*p+97.2659).*p;
	depth=depth./gr;
%
function a=datestrj(b)
% 
% print julian date string
%
a=datestr(b-julian([0 1 0 0 0 0]));

%
%==============================================================
function [hours]=hms2h(h,m,s);
%HMS2H converts hours, minutes, and seconds to hours
%
%  Usage:  [hours]=hms2h(h,m,s);   or [hours]=hms2h(hhmmss);
%
if nargin== 1,
   hms=h;
   h=floor(hms/10000);
   ms=hms-h*10000;
   m=floor(ms/100);
   s=ms-m*100;
   hours=h+m/60+s/3600;
else
   hours=h+(m+s/60)/60;
end

%==============================================================
function a=cutstruct(a,ii)
% reduce array size in structure
lz=length(ii);
iok=find(ii==1);
if isstruct(a)
  fnames = fieldnames(a);
  for n=1:size(fnames,1)
   dummy = getfield(a,fnames{n});
   [ly,lx]=size(dummy);
   if ly==lz
    a=setfield(a,fnames{n},dummy(iok,:));
   elseif lx==lz
    a=setfield(a,fnames{n},dummy(:,iok));
   end
  end
end

