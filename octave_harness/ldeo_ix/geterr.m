%======================================================================
%                    G E T E R R . M 
%                    doc: Wed Jun 30 23:24:51 2004
%                    dlm: Thu Dec  7 09:50:52 2017
%                    (c) 2004 ladcp@
%                    uE-Info: 23 47 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% MODIFICATIONS BY ANT:
%  Jun 30, 2004: - BUG: bin numbering was wrong for asymmetric up/down
%	                bin setup
%  Jul  5, 2004: - added comments to debug depth mapping
%  Jul 16, 2004: - added global variable skip_figure_3 to workaround
%		   linux matlab bug
%  Oct  7, 2008: - extensively modified procfig 3 for version IX_6
%  Jun 29, 2011: - removed skp_figure_3
%		 - added ps.fig3_colormap, ps.fig3_err_y_axis, ps.fig3_avgerr
%  Jun 30, 2011: - fixed fig.3 middle column plot title for median plot
%  Jul  6, 2001: - fixed plot title
%  Jan 25, 2015: - separated uc/dc in bin-averaged residual plots
%  Jan 28, 2015: - BUG: figure legend typo
%  Dec  7, 2017: - BUG: btmi was set to nan for P6 station 94; fixed with
%			symptomatic work-around

function l=geterr(ps,dr,d,iplot)
% function l=geterr(dr,d,iplot)
% returns predicitons of U_ocean and
% U_ctd on the raw data grid
% 
% CTD velocity
if nargin<4, iplot=1; end

ps=setdefv(ps,'fig3_colormap',2);	% 1: jet	2: polar
ps=setdefv(ps,'fig3_err_y_axis',2);	% 1: bin#	2: depth
ps=setdefv(ps,'fig3_avgerr',2');	% 1: mean	2: median 

tim=dr.tim;
tim(1)=-1e30;
tim(end)=1e30;

uctd=-interp1q(tim',dr.uctd',d.time_jul');
vctd=-interp1q(tim',dr.vctd',d.time_jul');

[ib,it]=size(d.ru);

wm=medianan(d.rw,3);
wz=gradient(-d.z,d.time_jul*24*3600);
l.ru_ctd=meshgrid(uctd,1:ib)+d.weight*0;
l.rv_ctd=meshgrid(vctd,1:ib)+d.weight*0;
l.rw_ctd=meshgrid(wm,1:ib)+d.weight*0;
l.rw_ctd_z=meshgrid(wz,1:ib)+d.weight*0;
if existf(d,'wctd')
 l.rw_ctd_p=meshgrid(d.wctd,1:ib)+d.weight*0;
end

% OCEAN velocity

z=-d.izm+d.ru*0;
dz=diff(d.izm(:,1))';

ii=find(z>=min(dr.z) & z<=max(dr.z));

uoce=interp1q(dr.z,dr.u,z(ii));
voce=interp1q(dr.z,dr.v,z(ii));

[prof,bin]=meshgrid(1:it,1:ib);

l.ru_oce=full(sparse(bin(ii),prof(ii),uoce));
l.rv_oce=full(sparse(bin(ii),prof(ii),voce));
l.ru_oce(ib,it)=NaN;
l.rv_oce(ib,it)=NaN;
l.ru_oce=l.ru_oce+d.weight*0;
l.rv_oce=l.rv_oce+d.weight*0;
ii=find(l.ru_oce==0 & l.rv_oce==0);
l.ru_oce(ii)=NaN;
l.rv_oce(ii)=NaN;

% ocean velocity as a function of depth and time

					% ib is number of bins
					% it is number of times (super ensembles)
itm=meshgrid(1:it,1:ib);		% each of ib rows of itm contains 1:it

					% d.izm contains for each time (colums),
					% list of absolute depths for each bin
dzdo=mean(abs(diff(d.izm(d.izd,1))));	% dzdo contains sound-speed corrected
					% mean bin length of downlooker at surface
					% NB: at depth, bins are smaller, because
					%     of increased soundspeed!

if length(d.izu)>1			% uplooker bin length
 dzup=mean(abs(diff(d.izm(d.izu,1))));
else
 dzup=dzdo;
end
dz=min([dzdo dzup]);			% dz is min bin length near surface
iz=-(d.izm/dz);				% iz is d.izm with depth coordinate given
					% as distance from surface, measured in 
					% near-surface bin lengths 

					% d.ru contains super-ensemble velocities
					% dr.z contains output depth grid
ii=find(isfinite(d.ru) & iz>0 & iz<max(dr.z)/dz);
					% ii contains indices (valid for d.ru,
					% d.izm, iz, ...) with valid velocities,
					% inside the output depth grid

ij=find( iz>0 & iz<max(dr.z)/dz);	% ij contains same as ii but also for
					% invalid velocities

if abs(dzup-dzdo)>dzup*0.1
 disp([' sorry dz not constant loop ',int2str(length(ii)),' elements'])
 for j=1:length(ii)
  iiz=ceil(iz(ii(j)));
  iit=itm(ii(j));
  l.u_oce(iiz,iit)=d.ru(ii(j))-l.ru_ctd(ii(j));
  l.u_err(iiz,iit)=d.ru(ii(j))-l.ru_oce(ii(j))-l.ru_ctd(ii(j));
  l.v_oce(iiz,iit)=d.rv(ii(j))-l.rv_ctd(ii(j));
  l.v_err(iiz,iit)=d.rv(ii(j))-l.rv_oce(ii(j))-l.rv_ctd(ii(j));
  l.w_oce(iiz,iit)=d.rw(ii(j))-l.rw_ctd(ii(j));
  l.w_oce_z(iiz,iit)=d.rw(ii(j))-l.rw_ctd_z(ii(j));
  if existf(l,'rw_ctd_p')
   l.w_oce_p(iiz,iit)=d.rw(ii(j))-l.rw_ctd_p(ii(j));
  end
  if existf(d,'tg')
   l.tg_oce(iiz,iit)=d.tg(ii(j));
  end

  l.u_ocean(iiz,iit)=l.ru_oce(ii(j));
  l.v_ocean(iiz,iit)=l.rv_oce(ii(j));

  l.u_adcp(iiz,iit)=d.ru(ii(j));
  l.v_adcp(iiz,iit)=d.rv(ii(j));
 end
else % uplooker and downlooker bin sizes are equal
 l.u_oce=full(sparse(ceil(iz(ii)),itm(ii),d.ru(ii)-l.ru_ctd(ii)));
 l.u_err=full(sparse(ceil(iz(ii)),itm(ii),d.ru(ii)-l.ru_oce(ii)-l.ru_ctd(ii)));
 l.v_oce=full(sparse(ceil(iz(ii)),itm(ii),d.rv(ii)-l.rv_ctd(ii)));
 l.v_err=full(sparse(ceil(iz(ii)),itm(ii),d.rv(ii)-l.rv_oce(ii)-l.rv_ctd(ii)));
 l.w_oce=full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd(ii)));
 l.w_oce_z=full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd_z(ii)));
 if existf(l,'rw_ctd_p')
  l.w_oce_p=full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd_p(ii)));
 end
 if existf(d,'tg')
  l.tg_oce=full(sparse(ceil(iz(ij)),itm(ij),d.tg(ij)));
 end

 l.u_ocean=full(sparse(ceil(iz(ii)),itm(ii),l.ru_oce(ii)));
 l.v_ocean=full(sparse(ceil(iz(ii)),itm(ii),l.rv_oce(ii)));

 l.u_adcp=full(sparse(ceil(iz(ii)),itm(ii),d.ru(ii)));
 l.v_adcp=full(sparse(ceil(iz(ii)),itm(ii),d.rv(ii)));
end

ik=find(l.u_oce==0 & l.v_oce==0);
l.u_oce(ik)=NaN;
l.v_oce(ik)=NaN;
l.w_oce(ik)=NaN;
l.w_oce_z(ik)=NaN;
if existf(l,'rw_ctd_p')
 l.w_oce_p(ik)=NaN;
end
l.u_adcp(ik)=NaN;
l.v_adcp(ik)=NaN;
if existf(d,'tg')
 ik=find(l.tg_oce==0);
 l.tg_oce(ik)=NaN;
end

[lz,lt]=size(l.u_oce);
l.itv=1:lt;

l.z_oce=([1:lz]-.5)*dz;
l.u_oce_m=meannan(l.u_oce');
l.v_oce_m=meannan(l.v_oce');

l.u_oce_s=stdnan(l.u_oce');
l.v_oce_s=stdnan(l.v_oce');

l.ru_err=d.ru-l.ru_oce-l.ru_ctd;
l.rv_err=d.rv-l.rv_oce-l.rv_ctd;

l.izm=d.izm;

[lz,lt]=size(l.ru_err);
l.itv2=1:lt;

if iplot

% blank out shallow/deep estimates
ii=find(iz<0 | iz>max(dr.z)/dz);
d.ru(ii)=nan;
d.rv(ii)=nan;

   figure(3)
   clf
   orient landscape
   
% find downcast/upcast separation
% Dec 2017
%	- P6 station 94 bombs because the first statement sets bmti to nan
%	- checking the code revelaed that btmi is never used iwth l.u_ocea but, rather
%	  with l.ru_err and l.rv_err
%	- therefore I added the if statement as a symptomatic cop-out
%	- the 2nd statement may be more appropriate than the first one but because
%	  the software has worked for many years with the first statement I kept 
%	  it like that
  btmi = fix(median(find(isfinite(l.u_oce(end,:)))));
  if ~isfinite(btmi)
      btmi = fix(median(find(isfinite(l.ru_err(end,:)))));
  end

% define color map
   if ps.fig3_colormap == 2
     colormap(polarmap(21));
   else
     col=jet(128);
     col=([[1 1 1]; col]);
     colormap(col)
   end

   ib=1:size(l.ru_err,1);
   ib=ib-length(d.izu);

   subplot(231)
   if ps.fig3_err_y_axis == 2
     if ps.fig3_colormap == 2
       tmp = l.u_err; tmp(isnan(tmp)) = 0;
       pcolorn(l.itv,-l.z_oce,tmp), shading flat
     else
       pcolorn(l.itv,-l.z_oce,l.u_err), shading flat
     end
     ylabel('Depth [m]');
   else
     if ps.fig3_colormap == 2
       tmp = l.ru_err; tmp(isnan(tmp)) = 0;
       pcolorn(l.itv2,-ib,tmp), shading flat
     else
       pcolorn(l.itv2,-ib,l.ru_err), shading flat
     end
     ylabel('Bin #');
   end
   fac=meannan(l.u_oce_s);
   fac=max([fac,1e-2]);
   caxis([-3 3]*fac)
   colorbar
   xlabel('Super Ensemble #');
   title(sprintf('U-err std: %.03f',meannan(stdnan(l.ru_err'))))
   
   subplot(232)
   if ps.fig3_avgerr == 2
     plot(medianan(l.ru_err(:,1:btmi)')',-ib,'r',medianan(l.ru_err(:,btmi:end)')',-ib,'b');
     title('median(U-err) [r/b: down-/up-cast]')
   else
     plot(meanan(l.ru_err(:,1:btmi)')',-ib,'r',meanan(l.ru_err(:,btmi:end)')',-ib,'b');
     title('mean(U-err) [r/b: down-/up-cast]')
   end
   set(gca,'XLim',[-0.05 0.05]);
   set(gca,'Ylim',[-ib(end) -ib(1)]);
   set(gca,'Xtick',[-0.04:0.02:0.04]);
   grid
   xlabel('Residual [m/s]');
   ylabel('Bin #');
   
   subplot(233)
   if ps.fig3_colormap == 2
     tmp = l.u_oce; tmp(isnan(tmp)) = 0;
     pcolorn(l.itv,-l.z_oce,tmp), shading flat
   else
     pcolorn(l.itv,-l.z_oce,l.u_oce), shading flat
   end
   ca = caxis;
   if abs(ca(1)) > abs(ca(2))
    caxis([-abs(ca(1)) abs(ca(1))]);
   else 
    caxis([-abs(ca(2)) abs(ca(2))]);
   end
   if existf(dr,'zbot')
    hold on
    plot(-d.z+d.hbot,'.k')
    ax=axis;
    ax(4)=maxnan([-d.z+d.hbot,ax(4)]);
    axis(ax);
   end
   colorbar
   xlabel('Ensemble #');
   ylabel('Depth [m]');
   title('U_{oce}')
   
   subplot(234)
   if ps.fig3_err_y_axis == 2
     if ps.fig3_colormap == 2
       tmp = l.v_err; tmp(isnan(tmp)) = 0;
       pcolorn(l.itv,-l.z_oce,tmp), shading flat
     else
       pcolorn(l.itv,-l.z_oce,l.v_err), shading flat
     end
     ylabel('Depth [m]');
   else
     if ps.fig3_colormap == 2
       tmp = l.rv_err; tmp(isnan(tmp)) = 0;
       pcolorn(l.itv2,-ib,tmp), shading flat
     else
       pcolorn(l.itv2,-ib,l.rv_err), shading flat
     end
     ylabel('Bin #');
   end
   fac=meannan(l.v_oce_s);
   fac=max([fac,1e-2]);
   caxis([-3 3]*fac)
   colorbar
   xlabel('Super Ensemble #');
   ylabel('Bin #');
   title(sprintf('V-err std: %.03f',meannan(stdnan(l.rv_err'))))
   
   subplot(235)
   if ps.fig3_avgerr == 2
     plot(medianan(l.rv_err(:,1:btmi)')',-ib,'r',medianan(l.rv_err(:,btmi:end)')',-ib,'b');
     title('median(V-err) [r/b: down-/up-cast]')
   else
     plot(meanan(l.rv_err(:,1:btmi)')',-ib,'r',meanan(l.rv_err(:,btmi:end)')',-ib,'b');
     title('mean(V-err) [r/b: down-/up-cast]')
   end
   set(gca,'XLim',[-0.05 0.05]);
   set(gca,'Ylim',[-ib(end) -ib(1)]);
   set(gca,'Xtick',[-0.04:0.02:0.04]);
   grid
   xlabel('Residual [m/s]');
   ylabel('Bin #');

   subplot(236)
   if ps.fig3_colormap == 2
     tmp = l.v_oce; tmp(isnan(tmp)) = 0;
     pcolorn(l.itv,-l.z_oce,tmp), shading flat
   else
     pcolorn(l.itv,-l.z_oce,l.v_oce), shading flat
   end
   ca = caxis;
   if abs(ca(1)) > abs(ca(2))
    caxis([-abs(ca(1)) abs(ca(1))]);
   else 
    caxis([-abs(ca(2)) abs(ca(2))]);
   end
   if existf(dr,'zbot')
    hold on
    plot(-d.z+d.hbot,'.k')
    ax=axis;
    ax(4)=maxnan([-d.z+d.hbot,ax(4)]);
    axis(ax);
   end
   colorbar
   xlabel('Ensemble #');
   ylabel('Depth [m]');
   title('V_{oce}')
   
   streamer([dr.name,'  Figure 3']);

% reset colormap
figure(11)
colormap(jet(128))

end


%======================================================================
%                    P O L A R M A P . M 
%                    doc: Tue Oct  7 11:03:28 2008
%                    dlm: Tue Oct  7 11:13:04 2008
%                    (c) 2008 A.M. Thurnherr
%                    uE-Info: 21 51 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

function map = polarmap(n)

if nargin<1, n=129, end;	% same as for jet()

map = ones(n,3);

firstred  = ceil(n/2) + 1;
lastblue = floor(n/2);

map([1:lastblue],1) = [0:lastblue-1]'/lastblue;
map([1:lastblue],2) = [0:lastblue-1]'/lastblue;
map([firstred:end],2) = [lastblue-1:-1:0]'/lastblue;
map([firstred:end],3) = [lastblue-1:-1:0]'/lastblue;
