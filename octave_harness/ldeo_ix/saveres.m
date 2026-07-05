%======================================================================
%                    S A V E R E S . M 
%                    doc: Fri Sep 26 15:55:18 2014
%                    dlm: Wed Mar 29 12:33:09 2017
%                    (c) 2014 A.M. Thurnherr
%                    uE-Info: 22 51 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

function []=saveres(dr,p,ps,f,d,att,da)
% function []=saveres(dr,p,ps,f,d,att,da)
%
% store LADCP result in RODB format
%

% MODIFICTIONS BY ANT:
%	Sep 26, 2014: - added support for p.orig (patch by Dan Torres)
%	Mar 29, 2017: - added att and da as requested by UH for archiving

if p.orig == 1
	eval(['save ',[f.res,'.mat'],' dr p ps f d da att'])
else
	eval(['save ',[f.res,'.mat'],' dr p ps f da att'])
end

% version 0.1	last change 28.6.2000
% open file
fid = fopen([f.res,'.lad'],'wt');

fprintf(fid,['Filename    = %s\n'],f.res);
fprintf(fid,['Date        = %4d/%2d/%2d\n'],p.time_start(1:3));
fprintf(fid,['Start_Time  = %2d:%2d:%2d\n'],p.time_start(4:6));
[lats,lons] = pos2str([p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60]);
fprintf(fid,['Start_Lat   = %s\n'],lats);
fprintf(fid,['Start_Lon   = %s\n'],lons);
fprintf(fid,['Deviation   = %f\n'],p.drot);
fprintf(fid,['Columns     = z:u:v:ev\n'],[]);
if existf(dr,'uerr')~=1
 dr.uerr=dr.u*NaN;
end
fprintf(fid,['%6.1f %6.3f %6.3f %6.3f \n'],[dr.z,dr.u,dr.v,dr.uerr]');

fclose(fid);

if existf(dr,'ubot')==1

% save bottom track data
% open file
fid = fopen([f.res,'.bot'],'wt');

fprintf(fid,['Filename    = %s\n'],f.res);
fprintf(fid,['Date        = %4d/%2d/%2d\n'],p.time_start(1:3));
fprintf(fid,['Start_Time  = %2d:%2d:%2d\n'],p.time_start(4:6));
[lats,lons] = pos2str([p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60]);
fprintf(fid,['Start_Lat   = %s\n'],lats);
fprintf(fid,['Start_Lon   = %s\n'],lons);
fprintf(fid,['Deviation   = %f\n'],p.drot);
fprintf(fid,['Bottom depth= %d\n'],fix(p.zbottom));
fprintf(fid,['Columns     = z:u:v:err\n'],[]);
fprintf(fid,['%6.1f %6.3f %6.3f %6.3f\n'],...
            [dr.zbot,dr.ubot,dr.vbot,dr.uerrbot]');

fclose(fid);

end
