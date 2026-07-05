%======================================================================
%                    P R E P I N V . M 
%                    doc: Wed Jan  7 16:46:29 2009
%                    dlm: Fri Apr  1 08:52:25 2016
%                    (c) 2009 A.M. Thurnherr
%                    uE-Info: 17 29 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jan  7, 2009: - tightened use of exist()
%   Oct 27, 2009: - color-coded down/upcast in fig. 6
%   Mar  4, 2010: - BUG: special code for shallow casts (no superensembles)
%			 was broken => removed
%   		  - BUG: down/upcast separation bombed wen no CTD data were
%			 available
%   Jun  9, 2014: - improved messages
%   Apr  1, 2016: - cosmetics

function [di,p,d]=prepinv(d,p,dr)
% function [di,p,d]=prepinv(d,p,dr)
% LADCP-2 processing software version 3.0
% prepare for inverse slolve
%
% - average velocities within depth range
% 
%  Martin Visbeck, LDEO, 6/10/99

% average over all ensembles in the depth range avdz 
%  set default avdz to one bin length
p=setdefv(p,'avdz',medianan(abs(diff(d.izm(:,1)))));
p=setdefv(p,'oversample',1.0);
p=setdefv(p,'avpercent',100);
p=setdefv(p,'rotup2down',0);
p=setdefv(p,'offsetup2down',0);
p=setdefv(p,'tilt_weight',10);
p=setdefv(p,'avens',NaN);
p=setdefv(p,'soundcorr',1);
p=setdefv(p,'nav_error',30);
d=setdefv(d,'slon',d.time_jul*NaN);
d=setdefv(d,'slat',d.time_jul*NaN);
p=setdefv(p,'superens_std_min',d.down.Single_Ping_Err/sqrt(d.down.Pings_per_Ensemble));

disp('PREPINV: prepare data for inversion, form Super-Ensembles')
if isfinite(p.avens)
 disp([' average profiles over (p.avens) ',num2str(p.avens),' ensembles'])
else
 disp([' average profiles over (p.avdz) ',num2str(p.avdz),' meter'])
end

% remove bottom "bad" track values 
% This can be because of range 
if isfinite(p.zbottom)
 p=setdefv(p,'btrk_range',[300 50]);
 % backward compatible
 if length(p.btrk_range)<2, p.btrk_range(2)=0; end

 ii=find(abs(d.z+p.zbottom)>max(p.btrk_range) & isfinite(d.bvel(:,1)') );
 if length(ii)>0
  disp([' discarded ',int2str(length(ii)),...
    ' bottom tracks velocities because of height above bottom > ',int2str(max(p.btrk_range))])
  d.bvel(ii,:)=NaN;
  d.bvels(ii,:)=NaN;
 end

 ii=find(abs(d.z+p.zbottom)<min(p.btrk_range) & isfinite(d.bvel(:,1)') );
 if length(ii)>0
  disp([' discarded ',int2str(length(ii)),...
    ' bottom tracks velocities because of height above bottom < ',int2str(min(p.btrk_range))])
  d.bvel(ii,:)=NaN;
  d.bvels(ii,:)=NaN;
 end

 ii=find(abs(d.z+p.zbottom-d.hbot)> 100 );
 if length(ii)>0
  disp([' discarded ',int2str(length(ii)),...
    ' bottom distances because of depth difference > 100'])
  d.hbot(ii)=NaN;
 end
else
 d.bvel(:,:)=NaN;
 d.bvels(:,:)=NaN;
end


% reduce weight for large tilts
if p.tilt_weight>0 & existf(d,'tilt_weight')~=1
 fac=1.-tanh(d.tilt/p.tilt_weight)/2;
 d.weight=d.weight.*meshgrid(fac,d.weight(:,1));
 disp([' reduce weight for larger tilts 0.5 at ',...
        num2str(p.tilt_weight),' degree'])
 d.tilt_weight=p.tilt_weight;
end

% prepare for heading averaging
i=sqrt(-1);
u1d=exp(-i*(d.hdg(1,:))*pi/180); 
if length(d.izu)>1
 u1u=exp(-i*(d.hdg(2,:))*pi/180); 
 % get mean heading offset from COMPASS comparison
 hoff=compoff(u1d,u1u);
 % hoff=angle(u1d/u1u)*180/pi;
 p.up_dn_comp_off=hoff;
 disp([' mean heading offset from compasses = ',num2str(hoff),' deg'])
 
 % check tilt sensors
 % rotate by compass offset
  diary off
  hoff2=fminsearch('checktilt',0,[],[d.rol(2,:);d.pit(2,:);d.rol(1,:);d.pit(1,:)]);
  diary on
  disp([' mean heading offset from pitch/roll = ',num2str(hoff2),' deg'])
  p.up_dn_pit_rol_comp_off=hoff2;
  [d.rol(3,:),d.pit(3,:)]=uvrot(d.rol(2,:),d.pit(2,:),-hoff2);
  p.up_dn_rol_off=mean(d.rol(1,:)-d.rol(3,:));
  p.up_dn_pit_off=mean(d.pit(1,:)-d.pit(3,:));
 
% plot compass and tilt meter differences
  figure(6)
  clf
  orient tall

% find downcast/upcast separation
  [btm,btmi] = min(d.z);

% plot compass comparison between up and down instrument
  da=-angle(u1d)+angle(u1u*exp(i*hoff*pi/180));
  ii=find(da>pi);
  da(ii)=da(ii)-2*pi;
  ii=find(da<-pi);
  da(ii)=da(ii)+2*pi;
  subplot(311)
  plot(d.hdg(1,1:btmi),da(1:btmi)*180/pi,'r.')
  hold on
  plot(d.hdg(1,btmi:end),da(btmi:end)*180/pi,'b.')
  ylabel('Heading difference up-down')
  title([' Heading offset : ',num2str(hoff),' (r/b: down-/upcast)'])
  xlabel('Heading down')
  grid
  dhmax=sort(abs(da*180/pi));
  dhmax=dhmax(fix(end*0.95));
  if dhmax > 15 & ~existf(d,'hrot')
   warn=[' Large compass deviation: ',num2str(dhmax)];
   p.warn(size(p.warn,1)+1,1:length(warn))=warn;
  end
  dhmax=max([10 dhmax*1.3]); 
  axis([0 360 [-1 1]*dhmax])

% plot tilt meter difference
  tiltmax=p.tiltmax(1);
  subplot(312)
  pa=d.pit(1,:)-d.pit(3,:)-p.up_dn_pit_off;
  plot(d.pit(1,1:btmi),pa(1:btmi),'r.')
  hold on
  plot(d.pit(1,btmi:end),pa(btmi:end),'b.')
  ylabel('Pitch difference up-down')
  xlabel('Pitch of down instrument')
  title([' Pitch offset : ',num2str(p.up_dn_pit_off)])
  grid;
  axis([[-1 1]*tiltmax [-1 1]*5])
 
  subplot(313)
  ra=d.rol(1,:)-d.rol(3,:)-p.up_dn_rol_off;
  plot(d.rol(1,1:btmi),ra(1:btmi),'r.')
  hold on
  plot(d.rol(1,btmi:end),ra(btmi:end),'b.')
  ylabel('Roll difference up-down')
  xlabel('Roll of down instrument')
  title(['Roll offset : ',num2str(p.up_dn_rol_off)])
  grid;
  axis([[-1 1]*tiltmax [-1 1]*5])
 
  streamer([p.name,' Figure 6']);
  pause(0.01)
end


% offset upward looking ADCP to downward looking ADCP
if (p.offsetup2down~=0 & length(d.zd)~=length(d.ru(:,1)) & exist('dr','var'))
  if p.rotup2down==2
   % will not rotate to match velocities and correct offset
   p.rotup2down=1;
  end
  i=sqrt(-1);
  disp(' remove first guess ocean velocity from raw data')
  % OCEAN velocity
  [ib,it]=size(d.ru);
  z=-d.izm+d.ru*0;
  ii=find(z>=min(dr.z) & z<=max(dr.z));
  l.uoce(:,1)=interp1q(dr.z,dr.u,z(ii));
  l.uoce(:,2)=interp1q(dr.z,dr.v,z(ii));
  [prof,bin]=meshgrid(1:it,1:ib);
  l.ru=full(sparse(bin(ii),prof(ii),l.uoce(:,1)));
  l.rv=full(sparse(bin(ii),prof(ii),l.uoce(:,2)));
  l.ru(ib,it)=NaN;
  l.rv(ib,it)=NaN;
  ii=find(~(z>=min(dr.z) & z<=max(dr.z)));
  l.ru(ii)=NaN;
  l.rv(ii)=NaN;

  uu=medianan(d.ru(d.izu,:)+i*d.rv(d.izu,:)...
             -l.ru(d.izu,:)-i*l.rv(d.izu,:)+d.weight(d.izu,:)*0,2);
  ud=medianan(d.ru(d.izd,:)+i*d.rv(d.izd,:)...
             -l.ru(d.izd,:)-i*l.rv(d.izd,:)+d.weight(d.izd,:)*0,2);
  clear l

 ii=find(~isfinite(uu+ud));
 uu(ii)=0;
 ud(ii)=0;
 uoff=(ud-uu)*p.offsetup2down; 
 uoffm(d.izu,:)=meshgrid(uoff/2,d.izu);
 uoffm(d.izd,:)=meshgrid(-uoff/2,d.izd);
 d.ru=d.ru+real(uoffm); 
 d.rv=d.rv+imag(uoffm); 
 if existf(d,'bvel')
  d.bvel(:,1)=d.bvel(:,1)+real(-uoff/2)';
  d.bvel(:,2)=d.bvel(:,2)+imag(-uoff/2)';
 end
 if existf(d,'uoff')
  d.uoff=d.uoff+uoff;
 else
  d.uoff=uoff;
 end
 % estimate tilt error that could explain difference
 % beam velocity ~ W
 wbeam=(abs(d.wm)*cos(d.down.Beam_angle*pi/180));
 ii=find(wbeam<0.4);
 wbeam(ii)=NaN;
 % uerror for 1 degree tilt error
 du=wbeam/sin((d.down.Beam_angle+0.5)*pi/180)-...
    wbeam/sin((d.down.Beam_angle-0.5)*pi/180);
 % linear projection to estimate tilt error
 d.tilterr=abs(d.uoff)./abs(du);

 figure(10), clf
 orient tall
 subplot(311)
 plot(real(d.uoff))
 uoffav=boxav(d.uoff',20);
 hold on
 plot([0.5:length(uoffav)]*20,real(uoffav),'-r')
 grid
 axis tight
 ax=axis;
 ax(3:4)=[-1 1]*max(medianan(abs(real(uoffav)))*6,0.1);
 axis(ax)

 title('U offset [m/s]')
 subplot(312)
 plot(imag(d.uoff))
 hold on
 plot([0.5:length(uoffav)]*20,imag(uoffav),'-r')
 grid
 axis tight
 ax=axis;
 ax(3:4)=[-1 1]*max(medianan(abs(imag(uoffav)))*6,0.1);
 axis(ax)
 title('V offset [m/s]')

 subplot(313)
 plot(d.tilterr)
 hold on
 tilterrav=boxav(d.tilterr',20);
 plot([0.5:length(tilterrav)]*20,tilterrav,'-r')
 axis tight
 ax=axis;
 ax(3:4)=[0 1]*medianan(d.tilterr)*6;
 if isfinite(sum(ax))
  ax(4)=max(3,ax(4));
  axis(ax)
 end
 grid
 title('Tilt error [degree] consistent with offset')
  
 streamer([p.name,' Figure 10']);
 pause(0.01)

 disp(' adjusted for velocity offset in up and down looking ADCP')

end

if existf(p,'drot')
  if ~isfinite(p.drot)
   warn=[' magnetic deviation given is NAN '];
   p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
  end
else
   warn=[' NO magnetic deviation given '];
   p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
end
  

 % check if rotated already

% rotate upward looking ADCP to downward looking ADCP
if (p.rotup2down~=0 & length(d.zd)~=length(d.ru(:,1)))

  da=angle(u1d)-angle(u1u*exp(-i*hoff*pi/180));
  ii=find(da>pi); da(ii)=da(ii)-2*pi; ii=find(da<-pi); da(ii)=da(ii)+2*pi;
  d.diff_hdg=da;

  u1uc=exp(-i*(d.hdg(2,:)-hoff)*pi/180);
  hrotcomp=angle(u1uc./u1d)*180/pi;
  d.rot_comp=hrotcomp;

  % get mean heading offset from velocities
  if exist('dr','var')
   disp(' remove first guess ocean velocity from raw data')
   % OCEAN velocity
   [ib,it]=size(d.ru);
   z=-d.izm+d.ru*0;
   ii=find(z>=min(dr.z) & z<=max(dr.z));
   l.uoce(:,1)=interp1q(dr.z,dr.u,z(ii));
   l.uoce(:,2)=interp1q(dr.z,dr.v,z(ii));
   [prof,bin]=meshgrid(1:it,1:ib);
   l.ru=full(sparse(bin(ii),prof(ii),l.uoce(:,1)));
   l.rv=full(sparse(bin(ii),prof(ii),l.uoce(:,2)));
   l.ru(ib,it)=NaN;
   l.rv(ib,it)=NaN;
   ii=find(~(z>=min(dr.z) & z<=max(dr.z)));
   l.ru(ii)=NaN;
   l.rv(ii)=NaN;

   uu=medianan(d.ru(d.izu,:)+i*d.rv(d.izu,:)...
              -l.ru(d.izu,:)-i*l.rv(d.izu,:)+d.weight(d.izu,:)*0,2);
   ud=medianan(d.ru(d.izd,:)+i*d.rv(d.izd,:)...
              -l.ru(d.izd,:)-i*l.rv(d.izd,:)+d.weight(d.izd,:)*0,2);
   clear l

  else
   iz=1:4;
   uu=medianan(d.ru(d.izu(iz),:)+i*d.rv(d.izu(iz),:)+d.weight(d.izu(iz),:),1);
   ud=medianan(d.ru(d.izd(iz),:)+i*d.rv(d.izd(iz),:)+d.weight(d.izd(iz),:),1);
  end
%  try to take speed into account
  hrotvel=angle(uu./ud)*180/pi;
  d.rot_vel=hrotvel;

%  decide which to use 
  if p.rotup2down==1
   hrot=hrotcomp;
   disp(' rot up2down use mean up/down compass')
  elseif p.rotup2down==2
   hrot=hrotvel;
   disp(' rot up2down use velocities')
  elseif p.rotup2down==3
   hrot=hrotcomp;
   disp(' rot up2down use down compass')
  elseif p.rotup2down==4
   hrot=hrotcomp;
   disp(' rot up2down use up compass')
  else
   hrot=d.ru(1,:)*0;
   disp(' not sure what you want rot=0')
  end  

% plot what I am about to do
  figure(5)
  clf
  orient tall
  plot([hrotcomp-90;hrotvel-180;hrot]')
  hold on
  ax=axis;
  plot(ax(1:2),[0 0],'-k',ax(1:2),[0 0]-90,'-k',ax(1:2),[0 0]-180,'-k') 
  text(1,20,'rotation used ')
  text(1,-60,'rotation compass')
  text(1,-150,'rotation velocity')
  ax(3)=-250;
  axis(ax)
  ylabel(' heading correction for uplooking')
  streamer([p.name,' Figure 5']);
  pause(.01)

 if existf(d,'hrot')
  disp(' rotated earlier, use difference ')
  oldrot=d.hrot;
 else
  oldrot=0;
 end

 % apply rotation
  hrotm=meshgrid(hrot-oldrot,d.izd);
   % dont rotate NaNs
  ii=find(~isfinite(hrotm));
  hrotm(ii)=meannan(hrot-oldrot);
 % rotate both by half difference
  if p.rotup2down==4
   [ru,rv]=uvrot(d.ru(d.izd,:),d.rv(d.izd,:),-hrotm);
   d.ru(d.izd,:)=ru;
   d.rv(d.izd,:)=rv;
   [bu,bv]=uvrot(d.bvel(:,1),d.bvel(:,2),-(hrot-oldrot)');
   d.bvel(:,1)=bu;
   d.bvel(:,2)=bv;
  elseif p.rotup2down~=3
   [ru,rv]=uvrot(d.ru(d.izd,:),d.rv(d.izd,:),-hrotm/2);
   d.ru(d.izd,:)=ru;
   d.rv(d.izd,:)=rv;
   [bu,bv]=uvrot(d.bvel(:,1),d.bvel(:,2),-(hrot-oldrot)'/2);
   d.bvel(:,1)=bu;
   d.bvel(:,2)=bv;
  end
 
  hrotm=meshgrid(hrot-oldrot,d.izu);
   % dont rotate NaNs
  ii=find(~isfinite(hrotm));
  hrotm(ii)=meannan(hrot);
 % rotate both by half difference
  if p.rotup2down==4
   [ru,rv]=uvrot(d.ru(d.izu,:),d.rv(d.izu,:),hrotm);
   d.ru(d.izu,:)=ru;
   d.rv(d.izu,:)=rv;
  elseif p.rotup2down~=3
   [ru,rv]=uvrot(d.ru(d.izu,:),d.rv(d.izu,:),hrotm/2);
   d.ru(d.izu,:)=ru;
   d.rv(d.izu,:)=rv;
  end

  d.hrot=hrot;

end


% sound speed correction
if p.soundcorr==1
   if existf(d,'p')~=1
     d.p=press(abs(d.z));
   end
   if existf(d,'ss')~=1
    disp(' make soundspeed based on pressure and ADCP temp')
    d.ss=sounds(d.p,d.temp(1,:),34.5);
   end
   if d.soundc==0
    disp(' correct velocities for sound speed ')
    sc=meshgrid(d.ss./d.sv(1,:),d.izd);
    d.ru(d.izd,:)=d.ru(d.izd,:).*sc;
    d.rv(d.izd,:)=d.rv(d.izd,:).*sc;
    d.rw(d.izd,:)=d.rw(d.izd,:).*sc;
    if length(d.zd)~=length(d.ru(:,1))
     sc=meshgrid(d.ss./d.sv(1,:),d.izu);
     d.ru(d.izu,:)=d.ru(d.izu,:).*sc;
     d.rv(d.izu,:)=d.rv(d.izu,:).*sc;
     d.rw(d.izu,:)=d.rw(d.izu,:).*sc;
    end
    d.soundc=1;
   else
    disp(' will not correct for sound speed twice')
   end
end
 

if (isnan(p.avdz) | p.avdz<=0 ) & (isnan(p.avens) | p.avens < 2)
 disp(' avdz=NAN  => No pre-averaging done !!!')
 di.ru=d.ru;
 di.rv=d.rv;
 di.ruvs=d.ru*0+p.superens_std_min;
 di.rw=d.rw;
 di.re=d.re;
 di.ts=d.ts;
 di.tg=d.tg;
 di.weight=d.weight;
 di.bvel=d.bvel';
 di.hbot=d.hbot;
% make up std
 di.bvels=d.bvel'*0+p.single_ping_accuracy;
 di.hdg=d.hdg;
 di.pit=d.pit;
 di.rol=d.rol;
 di.temp=d.temp;
 di.tsd=d.ts(d.izd(2),:);
 di.tsd_out=d.ts(d.izd(end),:);
 di.dtiv=d.z*0+1;
 di.time_jul=d.time_jul;
 di.z=d.z;
 di.izm=d.izm;
 di.slat=d.slat;
 di.slon=d.slon;
 di.izd=d.izd;
 di.izu=d.izu;

else

% default reference velocity bin range bin 2 + 3
 if length(d.izd)>2
   izr=[min(d.izd)+[1:2]];
 end
 if length(d.izu)>2
  izr=[izr,max(d.izu)-[1:2]];
 end
 d=setdefv(d,'izr',izr);
 di.izr=d.izr;
 izr=di.izr;

% remove reference velocity and then average ensembles
  disp(' remove reference velocity and average ensembles ')
  ilast=1;
  il=length(d.izm);
  im=0;
  ibin=1:p.nbins;

% big loop
  while ilast<il 
   im=im+1;
   i=ilast;
   if p.avens>0

%  fixed number of ensembles
    i1=ilast+[1:p.avens];
   else

%  set up index for ensembles to average based on depth
    ii=find(abs(d.izm(1,(ilast+1):il)-d.izm(1,ilast))>p.avdz);
    if length(ii)<1, ii=il-ilast; end
    i1=ilast+[1:ii(1)];
   end

   i1l=length(i1)/2*p.oversample;
   i1=round(mean(i1)+[-i1l:i1l]);
   ii=find(i1<1); i1(ii)=[];
   ii=find(i1>il); i1(ii)=[];
   if length(i1)==1, i1=[i1 i1]; end
   ilast=max(i1);
  

   w=d.weight(izr,i1)*0+1;
   w2=d.weight(:,i1)*0+1;

% U
   ur=medianan(d.ru(izr,i1).*w);
   ruav=meannan(ur);
   i3=find(isnan(ur));
   ur(i3)=i3*0;
   iav=round(length(ur)/200*p.avpercent);
   ur=meshgrid(ur,ibin);
   di.ru(:,im)=medianan([d.ru(:,i1).*w2-ur]',iav)'+ruav;
   rus=stdnan([d.ru(:,i1).*w2]')';
% V
   vr=medianan(d.rv(izr,i1).*w);
   rvav=meannan(vr);
   i3=find(isnan(vr));
   vr(i3)=i3*0;
   iav=round(length(vr)/200*p.avpercent);
   vr=meshgrid(vr,ibin);
   di.rv(:,im)=medianan([d.rv(:,i1).*w2-vr]',iav)'+rvav;
% estimate mean STD of U and V
   di.ruvs(:,im)=sqrt(rus.^2+stdnan([d.rv(:,i1).*w2]')'.^2);
% W
   wr=medianan(d.rw(izr,i1).*w);
   rwav=meannan(wr);
   i3=find(isnan(wr));
   wr(i3)=i3*0;
   iav=round(length(wr)/200*p.avpercent);
   wr=meshgrid(wr,ibin);
   di.rw(:,im)=medianan([d.rw(:,i1).*w2-wr]',iav)'+rwav;

%EA
   di.re(:,im)=meannan(d.re(:,i1)')';

%TS
   di.ts(:,im)=meannan(d.ts(:,i1)')';
   di.tg(:,im)=meannan(d.tg(:,i1)')';

% weight
   di.weight(:,im)=meannan(d.weight(:,i1)')';

% bottom track
   di.bvel(:,im)=meannan(d.bvel(i1,:))';
   bvel=d.bvel(i1,:);
%  remove mean vertical velocity from bottom track w prior to STD 
   bvel(:,3)=bvel(:,3)-wr(1,:)';
   di.bvels(:,im)=stdnan(bvel)';
%  distance of bottom
   di.hbot(im)=meannan(d.hbot(i1));

% bin depth
   di.izm(:,im)=mean(d.izm(:,i1)')';

% heading
   di.hdg(1,im)=-angle(mean(u1d(i1)))*180/pi;
   if length(d.izu)>1
     di.hdg(2,im)=-angle(mean(u1u(i1)))*180/pi;
   end

% pitch and roll
   di.pit(:,im)=mean(d.pit(:,i1),2);
   di.rol(:,im)=mean(d.rol(:,i1),2);

% target strength
   di.tsd(im)=mean(d.ts(d.izd(2),i1),2);
   di.tsd_out(im)=mean(d.ts(d.izd(end),i1),2);

% target strength
   di.temp(im)=mean(d.temp(i1));

% ships position
   di.slon(im)=medianan(d.slon(i1));
   di.slat(im)=medianan(d.slat(i1));

% number of ensembles
   di.dtiv(im)=length(i1);

% time
   di.time_jul(im)=meannan(d.time_jul(i1));

% depth
   di.z(im)=meannan(d.z(i1));
% end of big loop
end

% adjust compass to RDI definition [0-360]
di.hdg=di.hdg+(di.hdg<0)*360;

% remove outlier
di.izd=d.izd;
di.izu=d.izu;
[di,p]=outlier(di,p);

% remove bottom track data with single ping or large wstd
if isfinite(p.zbottom)
 % remove single ping bottom track ensembles
 ii=find(prod(di.bvels(1:3,:))==0);
 disp([' found ',int2str(length(ii)),' bottom track std==0 set to 0.1 m/s'])
 di.bvels(:,ii)=0.1;

 ii=find(di.bvels(3,:)>0);
 if length(ii)>0
  disp([' found ',int2str(length(ii)),' finite bottom track ensembles'])
  p=setdefv(p,'btrk_wstd',median(di.bvels(3,ii))*2);
  ii=find(di.bvels(3,:)>p.btrk_wstd | di.bvels(3,:)==0);
  disp([' discarded ',int2str(length(ii)),...
    ' bottom tracks velocities because of wstd  > ',num2str(p.btrk_wstd)])
  di.bvel(:,ii)=NaN;
  di.bvels(:,ii)=NaN;
 else
  ii=find(isfinite(d.bvel(3,:)));
  disp([' found no valid bottom track ensemble from ',int2str(length(ii)),...
        ' finite raw bottom tracks '])
 end 
end


% remove ensembles without velocity data
ii=find(~isfinite(maxnan(di.ru)));

if length(ii)>0
 disp([' removed ',int2str(length(ii)),...
        ' non finite super ensembles'])
 di.rw(:,ii)=[];
 di.ru(:,ii)=[];
 di.rv(:,ii)=[];
 di.ruvs(:,ii)=[];
 di.bvel(:,ii)=[];
 di.bvels(:,ii)=[];
 di.weight(:,ii)=[];
 di.izm(:,ii)=[];
 di.z(ii)=[];
 di.hdg(:,ii)=[];
 di.pit(:,ii)=[];
 di.rol(:,ii)=[];
 di.time_jul(ii)=[];
 di.dtiv(ii)=[];
 di.hbot(ii)=[];
 di.tsd(ii)=[];
 di.temp(ii)=[];
 di.slon(ii)=[];
 di.slat(ii)=[];
end

% check for small std's
% blank out data using weight 
di.ruvs=di.ruvs+di.weight*0;
ii=find(di.ruvs==0);
di.weight(ii)=nan;
disp([' set ',int2str(length(ii)),' weight values to nan  because super ensemble std =0 '])
di.ruvs=di.ruvs+di.weight*0;
ii=find(di.ruvs<p.superens_std_min);
di.ruvs(ii)=p.superens_std_min;
disp([' set ',int2str(length(ii)),' values to minimum super ensemble std ',...
     num2str(p.superens_std_min)]);

disp([' reduced profile length = ',int2str(length(di.z)),' super-ensemble bins'])
if length(di.z)<5
	error('not enough data to process station ')
end
end

% time difference in seconds
dtt=diff(di.time_jul)*24*3600;
di.dt=mean([dtt([1,1:end]);dtt([1:end,end])]);

% set positons
p=setdefv(p,'poss',[0 0 0 0]);
p=setdefv(p,'pose',[0 0 0 0]);

slat = p.poss(1)+p.poss(2)/60.0;
elat = p.pose(1)+p.pose(2)/60.0;
slon = p.poss(3)+p.poss(4)/60.0;
elon = p.pose(3)+p.pose(4)/60.0;

% compute ships drift
sjul =julian(p.time_start);
ejul =julian(p.time_end);

p.dt_profile = (ejul-sjul)*24*3600;

dlat = elat - slat;
dlon = elon - slon;

p.lat= ( slat + elat ) /2.0;
p.lon= ( slon + elon ) /2.0;

p.xdisp = dlon * cos(p.lat*pi/180) * 60.0 * 1852.0;
p.ydisp = dlat * 60.0 * 1852.0;

p.uship= p.xdisp / p.dt_profile;
p.vship= p.ydisp / p.dt_profile;


%--------------------------------------------------
function  hoff=compoff(u1,u2)
% compute mean compass offset 

h1=-angle(u1)*180/pi;
h1=h1+(h1<0)*360;
nhead=36;
dhead=360/2/nhead;
h0=linspace(5,355,nhead);

for i=1:nhead
 ii=find(abs(h1-h0(i))<=dhead);
 if length(ii)>1
  u1a(i)=mean(u1(ii));
  u2a(i)=mean(u2(ii));
 else
  u1a(i)=NaN;
  u2a(i)=NaN;
 end
end

ii=find(isfinite(u1a+u2a));

if length(ii)>0
 hoff=angle(u1a(ii)/u2a(ii))*180/pi;
else
 hoff=0;
end

return

%-----------------------------------------------------
function y=boxav(x,n)
% boxaverage and subsample
[ld,lv]=size(x);
in=fix(ld/n);

for i=1:lv
 xm=reshape(x(1:(in*n),i),n,in);
 y(:,i)=meannan(xm)';
end
