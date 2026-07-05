%======================================================================
%                    L O A D S A D C P . M 
%                    doc: Sun Jun 27 23:42:04 2004
%                    dlm: Fri Jun 12 11:35:40 2015
%                    (c) 2004 ladcp@
%                    uE-Info: 97 33 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%  Jun 30, 2004: - added Figure 9 title
%  Jul  9, 2004: - made lack of SADCP data in time window into real warning
%  Jan  7, 2009: - tightened use of exist()
%  Mar  4, 2010: - changed default of p.sadcp_dtok from 5min to zero
%  May 22, 2012: - removed code that took GPS information from SADCP data, 
%		   because these data are unlikely to be accurate enough
%		   for the ship-drift constraint; if they are, the user
%		   should verify and make a GPS file during pre-processing
%  Jun 12, 2015: - added code to set p.lat and p.lon when only SADCP 
%		   nav data are available

function   [di,p]=loadsadcp(f,di,p)
% function   [di,p]=loadsadcp(f,di,p)
%
%
%
%  need to make array  di.svel=[depth,U,V,Vel_error];
%       will use array di.sadcp_lon
%                      di.sadcp_lat
%
% 

% find SADCP profiles within time of station +/- slack (in fractional days)
p = setdefv(p,'sadcp_dtok',0);

if existf(f,'sadcp')==1
 if exist(f.sadcp,'file')

  disp(['LOADSADCP: load SADCP data file ',f.sadcp])
  load(f.sadcp);
  
% By Now you should have:
%           tim_sadcp(t) - Julian Days
%           lat_sadcp(t) - Degrees N
%           lon_sadcp(t) - Degrees E 
%           u_sadcp(z,t) - m/s
%           v_sadcp(z,t) - m/s
%           z_sadcp(z,1) - Meter (Positive Depth) 

  ii=find(tim_sadcp>(julian(p.time_start)-p.sadcp_dtok) & tim_sadcp<(julian(p.time_end)+p.sadcp_dtok));
  if length(ii)<1 
   warn = ' no SADCP data found in time window';
   p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
   disp(warn)
   return; 
  end

% check if position matches
  if abs(p.lon)+abs(p.lat)~=0 & (abs(mean(lon_sadcp(ii))-p.lon)>0.1/cos(p.lat*pi/360) |...
      abs(mean(lat_sadcp(ii))-p.lat)>0.1)
   disp(' position of SADCP data is more than 0.1 degree away from LADCP')
   disp(' no SADCP data used')
   
   figure(9)
   clf
   orient tall
   subplot(211)
   text(0,0.5,' SADCP data is more than 0.1 degree away from LADCP',...
            'color','r','fontsize',16,'fontweight','bold')
   axis off
   subplot(212)
   plot(lon_sadcp(ii),lat_sadcp(ii),'rp')
   hold on
   if existf(di,'slon')
    plot(di.slon,di.slat,'.g')
   end
   if abs(p.lon)+abs(p.lat)~=0
    plot(p.lon,p.lat,'pr')
    plot(p.poss(3)+p.poss(4)/60,p.poss(1)+p.poss(2)/60,'bp','markersize',15)
    plot(p.pose(3)+p.pose(4)/60,p.pose(1)+p.pose(2)/60,'kp','markersize',15)
   end
   title('ship nav (g.) start (bp) end (kp) SADCP (rp)')
   xlabel('longitude')
   ylabel('latitude')
   streamer([p.name,' Figure 9']);
   pause(0.001)
   return; 
  end % if no SADCP data in LADCP region
  
  % interpolate SADCP navigation to LADCP time series
  di.sadcp_lon=interp1(tim_sadcp,lon_sadcp(:),di.time_jul);
  di.sadcp_lat=interp1(tim_sadcp,lat_sadcp(:),di.time_jul);

  % set position from SADCP nav
  % NB: this is only used to set lat lon in the output files; magdev is NOT calculated from these
  if abs(p.lon)+abs(p.lat)==0
    p.lon = medianan(di.sadcp_lon);
    p.lat = medianan(di.sadcp_lat);
  end

%  if abs(p.lon)+abs(p.lat)==0
%   error('as of version IX_9, using GPS info from SADCP data stream is no longer supported');
%   slat=di.sadcp_lat(1);
%   slon=di.sadcp_lon(1);
%   elat=di.sadcp_lat(end);
%   elon=di.sadcp_lon(end);
%   p.poss=[fix(slat), (slat-fix(slat))*60, fix(slon), (slon-fix(slon))*60];
%   p.pose=[fix(elat), (slat-fix(elat))*60, fix(elon), (slon-fix(elon))*60];
%  end
 
%% if no other ship navigation exists, use SADCP navigation
%  if existf(di,'slon')==0
%   error('as of version IX_9, using GPS info from SADCP data stream is no longer supported');
%   di.slon=di.sadcp_lon;
%   di.slat=di.sadcp_lat;
%  else
%   if sum(isfinite(di.slon+di.slat))==0
%    error('as of version IX_9, using GPS info from SADCP data stream is no longer supported');
%    di.slon=di.sadcp_lon;
%    di.slat=di.sadcp_lat;
%   end
%  end

  if length(find(isfinite(u_sadcp(:,ii))))<1, 
   disp(' no finite SADCP data found ')
   return; 
  end

  disp([' found ',int2str(length(ii)),' SADCP profiles ']) 
  u=squeeze(u_sadcp(:,ii));
  v=squeeze(v_sadcp(:,ii));

% compute velocity standard deviation
  if numel(u) == length(u)
     v_err = u*0+0.1;
  else
     v_err=(nstd(u')+nstd(v'))';
     u=meannan(u')';
     v=meannan(v')';
  end
  ij=find(v_err==0);
  v_err(ij)=0.1;
  nvel=sum(isfinite(u+v)')'+v_err*0;
  v_err=v_err*max(nvel)./nvel;
  z=z_sadcp;

% make output array
  izok=find(isfinite(u+v) & z<(p.maxdepth));
  if isfinite(p.zbottom),
   izok=izok(find(z(izok)<(p.zbottom-30)));
  end
  di.svel=[z(izok),u(izok),v(izok),v_err(izok)];
 else
  disp([' can not find SADCP data file:',f.sadcp])
  return
 end
else
  disp(' no SADCP data file given')
  return
end

% plot some results

figure(9)
clf
subplot(211)
plot(squeeze(u_sadcp(:,ii)),-z_sadcp,'-r')
hold on
plot(squeeze(v_sadcp(:,ii)),-z_sadcp,'--g')
ax=axis;
axis(ax)
grid
streamer(p.name);
ylabel('depth [m]')
xlabel('velocity [m/s]')
title('SADCP U(-r) and V(--g) profiles')

subplot(212)
plot(di.slon,di.slat,'.g')
hold on
plot(lon_sadcp(ii),lat_sadcp(ii),'rp')
plot(p.poss(3)+p.poss(4)/60,p.poss(1)+p.poss(2)/60,'bp','markersize',15)
plot(p.pose(3)+p.pose(4)/60,p.pose(1)+p.pose(2)/60,'kp','markersize',15)
title('ship nav (g.) start (bp) end (kp) SADCP (rp)')
xlabel('longitude')
ylabel('latitude')
pause(0.001)


% ------------------------------------------------------------
function y = nstd(x,flag,dim)
%NSTD   Standard deviation, ignoring NaN.
%   Same as STD, but NaN's are ignored.
%

%   Copyright (c) 1997 by Toby Driscoll.
%   Adapted from STD.M, written by The MathWorks, Inc.
%   added backward compatibility	G.Krahmann, LODYC Paris


  if nargin<2, flag = 0; end
  if nargin<3, 
    dim = min(find(size(x)~=1));
    if isempty(dim), dim = 1; end
  end

% Avoid divide by zero.
  if size(x,dim)==1, y = zeros(size(x)); return, end

  tile = ones(1,max(ndims(x),dim));
  tile(dim) = size(x,dim);

  xc = x - repmat(meannan(x),tile);  % Remove mean
  mask = isnan(xc);
  xc(mask) = 0;
  s = sum(~mask,dim);
  s(s==0) = NaN;
  if flag,
    y = sqrt(sum(conj(xc).*xc,dim)./s);
  else
    z = (s==1);
    s(z) = Inf; 
    y = sqrt(sum(conj(xc).*xc,dim)./(s-1));
  end

