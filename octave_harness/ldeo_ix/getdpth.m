function [d,p]=getdpth(d,p)
% function [d,p]=getdpth(d,p)
% LADCP-2 processing software v 5.0
%
%  -- make depth from raw data 
% Matin Visbeck and Gerd Krahmann April-2000
% changed December 2002

disp('GETDEPTH: Integrating depth from vertical velocity')

% set default start, deepst  and end depth
p=setdefv(p,'zpar',[0 NaN 0]);
p=setdefv(p,'cut',15);
p=setdefv(p,'dzbelow',[2 -1]*medianan(abs(diff(d.zd))));
p.zpar=abs(p.zpar);

% get time difference for w-values
dt=diff(d.time_jul)*24*3600;
dt=mean([0,dt;dt,0]);

d.izmflag = d.rw*0;

p=setdefv(p,'ctddepth',0);

figure(4)
orient tall


% two sweeps to properly remove bottom values

d=getmeanw(d,p);

for n=1:2

  % set depth below which to delete all data 
  dzbelow=p.dzbelow(n);

  disp([' starting run ',int2str(n),' to find bottom depth'])

  % integrate vertical velocity to obtain depth
  %  assume zero depth at beginning and end of cast
  %  first sweep to get sound speed profile


  if n==1
   zzd=cumsum(d.wm.*dt);
   [zzdmax,zzdimax]=max(zzd);
   disp([' maxdepth form down-trace is ',num2str(zzdmax)])
   zzu=fliplr(cumsum(fliplr(-d.wm.*dt)));
   [zzumax,zzuimax]=max(zzu);
   disp([' maxdepth from up-trace   is ',num2str(zzumax)])
   % find mean maximum depth
   zzimax=round((zzdimax+zzuimax)/2);
   zzmax=(zzd(zzimax)+zzu(zzimax))/2;
   % merge up and down trace depth by forcing them to agree
   zz=[zzd(1:zzimax)/zzd(zzimax)*zzmax,zzu(zzimax+1:end)/zzu(zzimax)*zzmax];
   zz0=zz;
   if existf(d,'ctdprof_ss')
    disp(' take soundspeed from CTD profile ')
    zctd=d.ctdprof_z;
    zctd(1)=-1e5;
    zctd(end)=1e5;
    ss=interp1q(zctd,d.ctdprof_ss,zz')';
   else
    disp(' make soundspeed based on pressure and ADCP temp')
    pp=press(abs(zz));
    ss=sounds(pp,d.temp(1,:),34.5);
   end
   % sound speed correction
   if d.soundc==0
    sc=ss./d.sv(1,:);
    d.wm=d.wm.*sc;
    if existf(d,'hsurf')
     d.hsurf=d.hsurf.*sc;
    end
    if existf(d,'hbot')
     d.hbot=d.hbot.*sc;
    end
   end
  end

  % to surface corrections
  if ~isfinite(p.zpar(1)), p.zpar(1)=10; end
  if ~isfinite(p.zpar(3)), p.zpar(3)=10; end
  zsc=zz(1)-p.zpar(1);
  zec=zz(end)-p.zpar(3);

  [dum,ibot]=max(zz);
  % for very shallow stations turn of surface detection
  if dum<100, p.surfdist=0; disp(' shallow station no surface detection '),end
  ii=1:ibot;
  iok=ii(find(zz(ii)<200 & zz(ii)>30 ));
  if length(iok)>2
   iok2=1:iok(end);
   subplot(321)
   plot(iok2,-zz(iok2),'.r',iok2,iok2*0,'-k'),
   hold on
   title('start depth (.red) (.blue) Surface distance')
   xlabel('time in ensembles')
   ylabel('depth in meter')
   ax=axis;
   ax(4)=max(ax(4),20);
   axis(ax);
  end

  ii=ibot:length(zz);
  iok=ii(find(zz(ii)<200 & zz(ii)>30 ));
  if length(iok)>2
   iok2=iok(1):length(zz);
   subplot(322)
   plot(iok2,-zz(iok2),'.r',iok2,iok2*0,'-k'),
   hold on
   title('end depth')
   xlabel('time in ensembles')
   ylabel('depth in meter')
   ax=axis;
   ax(4)=max(ax(4),20);
   axis(ax);
  end

  % surface distance to find start depth
  if existf(d,'hsurf') & p.surfdist
   if sum(isfinite(d.hsurf))>10
    [zmax,ibot]=max(zz);
   % start depth
    ii=find(isfinite(d.hsurf(1:ibot)));
    iok=ii(find(zz(ii)<200 & zz(ii)>40 & abs(d.hsurf(ii)-zz(ii))<30));
    if ~isempty(iok)
     % use surface return for start depth
      start_depth_corr=median(d.hsurf(iok)-zz(iok));
    % temporary plot of surface detection
      iok2=1:iok(end);
      subplot(321)
      if n==1
       plot(iok,d.hsurf(iok)-zz(iok),'.',iok2,-zz(iok2),'.r'),
       hold on, plot(iok2,iok2*0+start_depth_corr,'--b') 
      else
       plot(iok,d.hsurf(iok)-zz(iok),'.',iok2,-zz(iok2),'.r'),
      end
      if std(d.hsurf(iok)-zz(iok))<20
       zsc=-start_depth_corr;
       disp([' start depth correction from surface return: ',int2str(zsc)])
      end
      ax=axis;
      ax(4)=max(ax(4),20);
      axis(ax);
    end

   % surface distance to find end depth
    ii=find(isfinite(d.hsurf(ibot:end)))+ibot-1;
    iok=ii(find(zz(ii)<200 & zz(ii)>40 & abs(d.hsurf(ii)-zz(ii))<30));
    if ~isempty(iok)
     % use surface return for end depth
      end_depth_corr=median(d.hsurf(iok)-zz(iok));
    % temporary plot of surface detection
      subplot(322)
      iok2=iok(1):length(zz);
      if n==1
       plot(iok,d.hsurf(iok)-zz(iok),'.',iok2,-zz(iok2),'.r'),
       hold on, plot(iok2,iok2*0+end_depth_corr,'--b') 
      else
       plot(iok,d.hsurf(iok)-zz(iok),'.',iok2,-zz(iok2),'.r'),
      end
      if std(d.hsurf(iok)-zz(iok))<20
       zec=-end_depth_corr;
       disp([' end depth correction from surface return: ',int2str(zec)])
      end
      ax=axis;
      ax(4)=max(ax(4),20);
      axis(ax);
    end
   end 
  end 

  disp([' last depth from int W is :',num2str(zz0(end))])
  disp([' should be                :',num2str(p.zpar(3))])
  % correct for start and end depth
  d.z_uncorr = -zz0;
  zz1=linspace(zsc,zec,length(zz));
  if length(zz1)>1
    zz=zz-zz1;
  else
    zz=zz-(zz(1)-p.zpar(1));
  end

  disp([' maximum depth from int W is :',num2str(max(zz))])
  disp([' should be                   :',num2str(p.zpar(2))])
  % correct for maximum depth
  if isfinite(p.zpar(2))
    zz=zz/max(zz)*p.zpar(2);
  end

  if p.ctddepth==0 
   % save results only if CTD-depth was not availabel
   d.z=-zz;
   p.ladcpdepth=1;
   disp(' use LADCP depth from integrated W')
  else
   p.ladcpdepth=0;
  end

  if p.ctddepth==1 & n>1
   d.z_ladcp=-zz;
   zz=-d.z;
   dz=d.z_ladcp-d.z;
   ii=find(isfinite(dz));
   p.ladcpr_CTD_depth_std=[mean(dz), std(dz)];
   disp(' use CTD time series depth ')
   disp([' LADCP minus CTD depth mean: ',num2str(p.ladcpr_CTD_depth_std(1)),...
          '  std: ',num2str(p.ladcpr_CTD_depth_std(2))]);
  end
  [p.maxdepth,ibottom]=max(-d.z);

  % bottom depth
  p.zbottom=NaN;
  subplot(312)
  iok=find((max(zz)-zz)<200);
  iok2=iok(1):iok(end);
  plot(iok2,d.z(iok2),'.r'),
  if sum(isfinite(d.hbot))>10
    % look for bottom only close to deepest CTD depth
    iok=find((max(zz)-zz)<200 & d.hbot>0 & abs(d.wm)>0);
    if ~isempty(iok)
      % fit polynom to bottom depth time series
      c=polyfit(iok,d.hbot(iok)-d.z(iok),1);
      if n>1
      % use deepest point to set bottom depth
       zbottomerr= polyval(c,iok)-(d.hbot(iok)-d.z(iok)) ;
       iok=iok(find(abs(zbottomerr)<1.5*std(zbottomerr) | abs(zbottomerr)<50 ));
       c=polyfit(iok,d.hbot(iok)-d.z(iok),2);
       zbottomerr= polyval(c,iok)-(d.hbot(iok)-d.z(iok)) ;
       p.zbottom=polyval(c,ibottom);
      else
       p.zbottom=medianan(d.hbot(iok)-d.z(iok));
       zbottomerr= p.zbottom-(d.hbot(iok)-d.z(iok)) ;
      end
      p.zbottomerror = medianan(abs(zbottomerr)); 
    % temporary plot of bottom detection
     iok2=iok(1):iok(end);
      subplot(312)
      plot(iok,-d.hbot(iok)+d.z(iok),'.',iok2,d.z(iok2),'.r'),
      hold on, plot(iok2,iok2*0-p.zbottom,'--k') 
      hold on, plot(iok2,-polyval(c,iok2),'-b')
    % remove outlier
      ii=find(abs(zbottomerr)>2*std(zbottomerr) | abs(zbottomerr)>100 );
      if n==1
       hold on, plot(iok(ii),-d.hbot(iok(ii))+d.z(iok(ii)),'*g'), 
      end
      title('bottom (--k)  bottom distance (.b)')
      xlabel('time in ensembles')
      ylabel('depth in meter')

      d.hbot(iok(ii))=NaN;
      d.bvel(iok(ii),:)=NaN;
    else
      p.zbottomerror = nan;
    end
   % check if bottom is shallower that maxctd-depth an
    if ((p.zbottom-p.maxdepth<-20 & isfinite(p.zbottom)) |...
       	p.zbottomerror > 20 )
      disp('  no bottom found')
      disp(['   given maximum profile depth : ',int2str(p.maxdepth)])
      disp(['   extracted bottom depth      : ',int2str(p.zbottom)])
      disp(['        bottom depth error     : ',int2str(p.zbottomerror)])
      p.zbottom=NaN;
    else
      disp(['  bottom found at ',int2str(p.zbottom),' +/- ',...
                                 int2str(p.zbottomerror),' m'])
      if (p.zbottom<p.maxdepth)
        disp('  extracted bottom within 10m above given maximum profile depth')
      end
    end
  end
  pause(0.1)

  [izm1,izm]=meshgrid([fliplr(d.zu),-d.zd],d.z);
  d.izm=izm'+izm1';

  % flag all data below bottom as bad
  if ~isnan(p.zbottom)
    ii = find(d.izm<-p.zbottom-dzbelow);
    d.izmflag(ii)=NaN; 
  end

  if length(d.zu)>0
   % flag all data close to the surface as bad
   ii = find(d.izm>-d.zu(1));
   d.izmflag(ii)=NaN; 
  end
  

end 


% set velocities deeper than bottom to NaN
bad = find( isnan(d.izmflag) & isfinite(d.ru) );
if ~isempty(bad)
  disp([' removing ',int2str(length(bad)),...
	' values  below recognized bottom'])
end

%d.ru = d.ru+d.izmflag;
%d.rv = d.rv+d.izmflag;
%d.rw = d.rw+d.izmflag;
d.weight = d.weight+d.izmflag;

% compute pressure from depth
d.p=press(abs(d.z));

if d.z(1)<-50, 
 warn=[' first LADCP depth is ',int2str(d.z(1))];
 disp(warn)
 p.warn(size(p.warn,1)+1,1:length(warn))=warn;
end

if d.z(end)<-50, 
 warn=[' last LADCP depth is ',int2str(d.z(end))];
 disp(warn)
 p.warn(size(p.warn,1)+1,1:length(warn))=warn;
end

% get sound speed time series if not available already
if ~existf(d,'ss')
 if existf(d,'ctd_ss')
    disp(' take soundspeed from CTD time series')
  d.ss=d.ctd_ss;
 elseif existf(d,'ctdprof_ss')
    disp(' take soundspeed from CTD profile ')
    zctd=d.ctdprof_z;
    zctd(1)=-1e5;
    zctd(end)=1e5;
    d.ss=interp1q(zctd,d.ctdprof_ss,d.z')';
 else
   if existf(d,'ctd_temp')
    disp(' make soundspeed based on pressure and CTD temp')
    d.ss=sounds(d.p,d.ctd_temp,34.5);
   else
    disp(' make soundspeed based on pressure and ADCP temp')
    d.ss=sounds(d.p,d.temp(1,:),34.5);
   end
 end
end

% correct velocity for sound speed
if d.soundc==0
    disp(' correct velocities for sound speed ')
    sc=meshgrid(d.ss./d.sv(1,:),d.izd);
    d.ru(d.izd,:)=d.ru(d.izd,:).*sc;
    d.rv(d.izd,:)=d.rv(d.izd,:).*sc;
    d.rw(d.izd,:)=d.rw(d.izd,:).*sc;
    if length(d.zd)~=length(d.ru(:,1))
     sc=meshgrid(d.ss./d.sv(2,:),d.izu);
     d.ru(d.izu,:)=d.ru(d.izu,:).*sc;
     d.rv(d.izu,:)=d.rv(d.izu,:).*sc;
     d.rw(d.izu,:)=d.rw(d.izu,:).*sc;
    end
    d.soundc=1;
else
    disp(' will not correct for sound speed twice')
end

% cut raw data to only include profile
i1=find(d.z(1:ibottom)>-p.cut);
i2=find(d.z(ibottom:end)>-p.cut)+ibottom-1;
if length(i1)==0
 i1=1;
end
if length(i2)==0
 i2=length(d.z);
end
ii=d.z*0;
ic=i1(end):i2(1);
ii(ic)=1;
if (sum(ii)~=length(ii)) & p.cut>0 & existf(d,'cutindx')~=1
 disp(' remove data at begining and end of cast')
 d=cutstruct(d,ii);
 subplot(321), ax=axis; plot([1 1]*ic(1),ax(3:4),'--k')
 subplot(322), ax=axis; plot([1 1]*ic(end),ax(3:4),'--k')
 %p.time_start=gregoria(d.time_jul(1));
 %p.time_end=gregoria(d.time_jul(end));
 p.zpar(1)=max([0,-d.z(1)]);
 p.zpar(3)=max([0,-d.z(end)]);
end

streamer([p.name,'  Figure 4']);
pause(0.01)

%-----------------------------------------------------------------
function d=getmeanw(d,p)
% function d=getmeanw(d,p)
  [d.wm1,dum]=medianan(d.rw(p.wizr,:)+d.izmflag(p.wizr,:),1);
  ii=find(sum(isfinite(dum))<2);
  d.wm1(ii)=NaN;
  d.wm=d.wm1;
  ii=2:(length(d.wm1)-1);
  d.wm(ii)=meannan([d.wm1(ii-1);d.wm1(ii);d.wm1(ii);d.wm1(ii+1)]);
  % try to replace bad W-data by propagating them from the ends
  for dummy=1:3
   ii=find(~isfinite(d.wm));
   if length(ii)==0, break, end
   dat=[];
   for nn=1:dummy
    i1=ii-nn;
    i=find(i1<1);
    i1(i)=1;
    i2=ii+nn;
    i=find(i2>length(d.wm));
    i2(i)=length(d.wm);
    dat=[dat;d.wm(i1);d.wm(i2)];
   end
   d.wm(ii)=meannan(dat);
  end
  % replace missing w-profiles by zero
  ii=find(~isfinite(d.wm));
  d.wm(ii)=0;


%==============================================================
function a=cutstruct(a,ii)
% reduce array size in structure
lz=length(ii);
iok=find(ii==1);
a.cutindx=[iok(1) iok(end)];
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
