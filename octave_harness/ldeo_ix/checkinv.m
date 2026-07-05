%======================================================================
%                    C H E C K I N V . M 
%                    doc: Tue Jun 13 14:29:15 2017
%                    dlm: Tue Jun 13 14:33:27 2017
%                    (c) 2017 A.M. Thurnherr
%                    uE-Info: 10 77 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%	Jun 13, 2017: Matlab incompatibility (legend had extraneous -1 final arg)

function p=checkinv(dr,di,de,der,p,ps)
% function p=checkinv(dr,di,de,der,p,ps)
% check inversion for consistency 
% Martin Visbeck 2004

uerrm=median(dr.uerr);
unoise=meannan(stdnan([der.ru_err ; der.rv_err])) ;

disp('CHECKINV all values are given in [m/s] ')
 disp([' Velocity profile error: ',num3str(uerrm,6,3),...
       '  should be about noise: ',num3str(unoise,6,3)])

% check bottom track
if existf(dr,'zbot') 
 if existf(de,'bvels') 

  % std of missmatch between bottom track and U_ctd
  ubso=([real(de.bvel)+dr.uctd]);
  vbso=([imag(de.bvel)+dr.vctd]);
  uvbso=meannan(sqrt(ubso.^2+vbso.^2));

 % std of bottom track 
  uvbsi=meannan(de.bvels);
 else
  disp(' do not know how uncertain the bottom track is')
  uvbso=nan;
  uvbsi=nan;
 end

 disp([' Check bottom track rms: ',num3str(uvbso,6,3),...
       '  should be smaller than ',num3str(uvbsi,6,3),' / ',num3str(ps.botfac,6,3)])

end

% check SADCP
if existf(dr,'z_sadcp')

 % std of missmatch between SADCP and LADCP velocity profile
 us=interp1(dr.z,dr.u,dr.z_sadcp);
 vs=interp1(dr.z,dr.v,dr.z_sadcp);
 usso=([us-dr.u_sadcp]);
 vsso=([vs-dr.v_sadcp]);
 uvsso=meannan(sqrt(usso.^2+vsso.^2));
 
 % given std of SADCP profile
 uvssi=medianan(dr.uerr_sadcp);

 disp([' Check SADCP        rms: ',num3str(uvsso,6,3),...
       '  should be smaller than ',num3str(uvssi,6,3),' / ',num3str(ps.sadcpfac,6,3)])

end

% check GPS
if existf(dr,'uship')

 % difference of mean ships drift from LADCP and GPS 
 dtiv=gradient(dr.tim);
 us=sum(dtiv.*dr.uship)/sum(dtiv);
 vs=sum(dtiv.*dr.vship)/sum(dtiv);
 uvgso=sqrt((us-p.uship).^2 + (vs-p.vship).^2);
 
 % computed uncertainty of GPS derived ships speed
 uvgsi=ps.barvelerr;

 disp([' GPS-LADCP ship spd diff:',num3str(uvgso,6,3),...
       '  should be smaller than ',num3str(uvgsi,6,3),' / ',num3str(ps.barofac,6,3)])
end

% plot inverse solution weights

if existf(de,'type_constraints')
 
 ic=find(sum([de.ctd_constraints,de.ocean_constraints]')~=0);  
 
 figure(12)
 clf
 orient tall

 ic2=3*length(ic);
 col=jet(ic2);
 ic3=[1:length(ic)]*3;
 ic3=[ic3(1:2:end),fliplr(ic3(2:2:end))];
 colormap(col(ic3,:));
 
 subplot(211)
 bar(dr.z,de.ocean_constraints(ic,:)','stack') 
% legend(de.type_constraints(ic,:),-1)
 legend(de.type_constraints(ic,:))
 ylabel('sum of weights')
 title('ocean velocity constraints')
 xlabel('depth [m]')
 axis tight

 subplot(212)
 bar(de.ctd_constraints(ic,:)','stack') 
% legend(de.type_constraints(ic,:),-1)
 legend(de.type_constraints(ic,:))
 title('CTD velocity constraints')
 ylabel('sum of weights')
 xlabel('super ensemble')
 axis tight
 
 streamer([p.name,'  Figure 12']);
end
