function [d,p]=getbtrack(d,p);
% function [d,p]=getbtrack(d,p);
%
% create own bottom track in addition to the one used before
%

% force own distance to RDI bottom track
p=setdefv(p,'bottomdist',0);
% "manual" bottom track for down looker
p=setdefv(p,'btrk_mode',3);
% p.btrk_ts is in dB to detect bottom above bin1 level (for own btm track)
p=setdefv(p,'btrk_ts',10);
% p.btrk_below gives bin used below target strength maximum
p=setdefv(p,'btrk_below',0.5);
% p.btrk_range gives range of allowed bottom track ranges
p=setdefv(p,'btrk_range',[300 50]);
% maximum allowed difference between reference layer W and W bottom track
p=setdefv(p,'btrk_wlim',0.05);

p=setdefv(p,'btrk_used',0);

disp('GETBTRACK creates own bottom track in addition to RDI')

% convert echo amplitude to relative target strength
% at=0.039; % attenuation dB/m for 150 kHz
% at=0.062; % attenuation dB/m for 300 kHz
if d.down.Frequency==300,
  p=setdefv(p,'ts_att_dn',0.06);
else
  p=setdefv(p,'ts_att_dn',0.039);
end
d.tg(d.izd,:)=targ(d.ts(d.izd,:),d.zd,p.ts_att_dn);

if length(d.izu)>0
 if d.up.Frequency==300,
  p=setdefv(p,'ts_att_up',0.06);
 else
  p=setdefv(p,'ts_att_up',0.039);
 end
 d.tg(d.izu,:)=targ(d.ts(d.izu,:),d.zu,p.ts_att_up);
end

% save bin number where down looker starts
ib1=d.izd(1)-1;
nbin=length(d.izd);

disp(['  in: p.btrk_mode ',int2str(p.btrk_mode),' and p.btrk_used ',int2str(p.btrk_used)])

if p.btrk_mode>=1
    if p.btrk_ts>0

      disp(' using increased bottom echo amplitudes to create bottom track')
      % don't use first bin
      fitb1=1;
      zd=d.zd(fitb1:end);
      ead=d.tg(d.izd(fitb1:end),:);
      % fit parabola to locate bottom
      [zmead,mead,imead]=localmax2(zd',ead);
      imead=imead+fitb1-1;
      dts=mead-ead(1,:);

      % decide which bin to use for bottom velocities
      dz=abs(diff(d.zd(1:2)));
      % imeadbv=round(imead+p.btrk_below);
      imeadbv=round((zmead-d.zd(1))/dz+1+p.btrk_below);

      if p.btrk_used==1
       if p.bottomdist==0
        % check RDI bottom track only if non zero
        ii=find(d.hbot<min(p.btrk_range));
        if length(ii)>0
         disp([' found ',int2str(length(ii)),' bottom depth below btrk_range ',...
            int2str(min(p.btrk_range))])
         d.bvel(ii,:)=nan;
         d.hbot(ii)=nan;
        end
        ii=find(d.hbot>max(p.btrk_range));
        if length(ii)>0
         disp([' found ',int2str(length(ii)),' bottom depth above btrk_range ',...
            int2str(max(p.btrk_range))])
         d.bvel(ii,:)=nan;
         d.hbot(ii)=nan;
        end
       end
       % save RDI bottom track
       d.bvel_rdi=d.bvel;
       d.hbot_rdi=d.hbot;
       ii=find(isfinite(d.bvel(:,1)+d.bvel(:,2)));
       if length(ii)<10, 
        disp(' found less than 10 RDI bottom track values, try own')
        p.btrk_used=0; 
       end
      end

      disp([' use ',num2str(p.btrk_below),...
         ' bins below maximum target strength for own bottom track velocity'])
      % locate acceptable bottom tracks (don't accept first/last bin)
      ii=find(dts>p.btrk_ts & ...
              zmead>min(p.btrk_range) & zmead<max(p.btrk_range) & ...
               imeadbv<(nbin-1) & imeadbv>fitb1);
 
      if length(ii)>0

      % save bottom distance data
       d.hbot_own=d.hbot+NaN;
       d.hbot_own(ii)=zmead(ii);

      % force bottom distance if RDI mode fails to report distance
       if p.bottomdist | p.btrk_mode==2 | p.btrk_used~=1
        d.hbot=d.hbot_own;
        disp([' created ',int2str(length(ii)),...
		' bottom distances'])
        if p.bottomdist, p.btrk_used = 12; end
       else
        disp([' created ',int2str(length(ii)),...
		' bottom distances keeping original'])
       end

       % make bottom velocity data
       bv=d.bvel+NaN;

       for j=1:length(ii)
          ji=ii(j);
          bv(ji,1)=medianan(d.ru(ib1+imeadbv(ji)+[-1,0,0,1],ji));
          bv(ji,2)=medianan(d.rv(ib1+imeadbv(ji)+[-1,0,0,1],ji));
          bv(ji,3)=medianan(d.rw(ib1+imeadbv(ji)+[-1,0,0,1],ji));
          bv(ji,4)=medianan(d.re(ib1+imeadbv(ji)+[-1,0,0,1],ji));
       end

       % check for W-bot 
       wref=medianan(d.rw(d.izd,:),2);
       ii=find(abs(wref'-bv(:,3))>p.btrk_wlim);
       disp([' removed ',int2str(length(ii)),...
           ' bottom track profiles W_btrk - W_ref difference > ',num2str(p.btrk_wlim)])
       bv(ii,:)=nan;

       % check for outlier
       bv=boutlier(bv,d.hbot_own,p);
       d.bvel_own=bv;
       ii=find(isfinite(bv(:,1)+bv(:,2)));

       if (p.btrk_used~=1 & p.btrk_used~=12) | p.btrk_mode==2
        p.btrk_used = 2;
        d.bvel=d.bvel_own;
        disp([' created ',int2str(length(ii)),...
		' bottom track data from normal velocities'])
       else
        disp([' created ',int2str(length(ii)),...
		' bottom track velocities keeping original'])
       end

      else
        if (p.btrk_used~=1 & p.btrk_used~=12)
         p.btrk_used = -1;
        end
        disp(' did not find any bottom echos to create own bottom track ')
      end
    else
      disp(' no valid own bottom track. Increase target strength difference ? ')
    end
else
 disp('force no bottom track data ')
 p.btrk_used = -1;
 d.bvel=d.bvel+nan;
 d.hbot=d.hbot+nan;
end
% summary output
disp([' out: p.btrk_mode ',int2str(p.btrk_mode),' and p.btrk_used ',int2str(p.btrk_used)])

%=================================
function [bvel,p] = boutlier(bvel,hbot,p);
% function [bvel,p] = boutlier(bvel,hbot,p);
%
%
% input  : bvel bottom track velocity
%           n       factor for outlier rejection
%                   n(1)*rms(data) first sweep
%                   n(2)*rms(data) second sweep etc.
%
% output :	d		changed LADCP analysis data structure
%
% version 0.1	last change 27.6.2000

n=p.outlier;
nblock=p.outlier_n;

dummyb = bvel*0;
if size(dummyb,2)~=4, disp(' no data '), return, end

si = size(dummyb);
sn = ceil(si(1)/nblock);


lob=length(find(isnan(dummyb)));

for i=1:length(n)
  for m=1:sn
    ind = (m-1)*nblock+[1:nblock];
    ii = find( ind<=si(1) );
    ind = ind(ii);
    % bottom track
     bvelt(ind,1)=bvel(ind,1)-medianan(bvel(ind,1));
     bu = find(abs(bvelt(ind,1))>n(i)*rms(bvelt(ind,1)));
     dummyb(ind(bu),:)=nan;
     bvelt(ind,2)=bvel(ind,2)-medianan(bvel(ind,2));
     bv = find(abs(bvelt(ind,2))>n(i)*rms(bvelt(ind,2)));
     dummyb(ind(bv),:)=nan;
     bw = find(abs(bvel(ind,3))>n(i)*rms(bvel(ind,3)));
     dummyb(ind(bw),:)=nan;
     hbot(ind)=hbot(ind)-medianan(hbot(ind));
     bh = find(abs(hbot(ind))>n(i)*rms(hbot(ind)));
     dummyb(ind(bh),:)=nan;
  end
  bvel = bvel + dummyb; 
  hbot = hbot + dummyb(:,1)'; 
end  
disp([' boutlier removed ',int2str(length(find(isnan(dummyb)))-lob),...
      ' bottom track velocities '])
return

%
function y = rms(x,dim)
%RMS    root-mean square. For vectors, RMS(x) returns the standard
%       deviation.  For matrices, RMS(X) is a row vector containing
%       the root-mean-square of each column. The difference to STD is
%       that here the mean is NOT removed.  
%       RMS(X,DIM) returns the root-mean-square of dimension DIM
%
%	See also STD,COV.
%

%       Uwe Send, IfM Kiel, Apr 1992
% added NaN handling   Gerd Krahmann, IfM Kiel, Oct 1993, Jun 1994
% removed bug in NaN handling   G.Krahmann, Aug 1994
% added compatibility to MATLAB 5	G.Krahmann, LODYC Paris, Jul 1997


  if nargin<2
    dim=min(find(size(x)>1));
  end

  if all(isnan(x))
    y=nan;
    return
  end    

  x = shiftdim(x,dim-1);
  s = size(x);
  so = s(1);
  s(1) = 1;

  for n = 1:prod(s)
    good = find(~isnan(x(:,n)));
    if ~isempty(good)
      y(1,n) = norm( x(good,n) ) / sqrt(length(good));
    else
      y(1,n) = NaN;
    end
  end
  y = reshape(y,s);
  y = shiftdim( y, ndims(x)-dim+1 );


%================================================
function [ts,bcs]=targ(ea,dis,at,bl,eas,ap)
% function [ts,bcs]=targ(ea,dis,at,bl,eas,ap)
% Target strength of EA for volume scatterer
% ea = echoamp in  dB
% dis = distance in  m
% at = attenuation dB/m
% bl = pulse/bin legth in  m
% eas = source level
% ap = aperature in degree
% M. Visbeck 2004

% make distance matrix if needed
[lr,lc]=size(ea);

if size(dis,2)==1 | size(dis,1)==1
 dis=dis(:);
 dis=repmat(dis,[1,lc]);
end

if nargin<3
 at=0.039; % attenuation dB/m for 150 kHz
 at=0.062; % attenuation dB/m for 300 kHz
end


%binlength
if nargin<4 , bl=median(abs(diff(dis(:,1)))); end

% source level in dB
if nargin <5, eas=100; end

% beam aperature in DEGREE convert to radian
if nargin <6, ap=2; end
al=ap*pi/180; 

% radius of top and bottom of each bin
r1=tan(al)*(dis-bl/2);
r2=tan(al)*(dis+bl/2);

% ensonified volume 
v=pi*bl/3*(r1.^2+r2.^2+r1.*r2);

% transmission loss
tl=20*log10(dis)+at*dis;

% target strength
ts=ea-eas+2*tl-10*log10(v);

if nargout>1
 %backscatter cross section
 bcs=exp(ts./10);
end
