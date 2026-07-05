%======================================================================
%                    S A V E P R O T . M 
%                    doc: Wed Jan  7 16:52:50 2009
%                    dlm: Wed Jan  7 16:53:04 2009
%                    (c) 2009 A.M. Thurnherr
%                    uE-Info: 11 0 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jan  7, 2009: - tightened use of exist()

% save a processing protocol sheet for LADCP profiles

disp(' ')
disp(['Writing processing protocol ',f.res,'.txt'])

% open file
fid = fopen([f.res,'.txt'],'wt');

fprintf(fid,'    LADCP processing protocol \n',[]);
fprintf(fid,' \n',[]);
fprintf(fid,'Station name         : %s \n',p.name);
fprintf(fid,' \n',[]);
fprintf(fid,'LADCP input files    : \n',[]);
if isfield(f,'ladcpdo')
  fprintf(fid,'down                 : %s \n',f.ladcpdo);
else
  fprintf(fid,'down                 : %s \n','none');
end
if isfield(f,'ladcpup')
  fprintf(fid,'up                   : %s \n',f.ladcpup);
else
  fprintf(fid,'up                   : %s \n','none');
end
fprintf(fid,' \n',[]);
fprintf(fid,'Command file info    : %s \n',[]);
fprintf(fid,' \n',[]);
if exist('pinfo','var')
  fprintf(fid,'Start time (CTD)     : %4d %2d %2d %2d %2d %2d \n',...
	pinfo.ladcpprofstart);
  fprintf(fid,'End time (CTD)       : %4d %2d %2d %2d %2d %2d \n',...
	pinfo.ladcpprofend);
  fprintf(fid,'Rec start (LADCP)    : %4d %2d %2d %2d %2d %2d \n',...
	pinfo.ladcpstart);
  fprintf(fid,'Rec end (LADCP)      : %4d %2d %2d %2d %2d %2d \n',...
	pinfo.ladcpend);
end
fprintf(fid,'Prof start (LADCP)   : %4d %2d %2d %2d %2d %2d \n',...
	p.time_start);
fprintf(fid,'Prof end (LADCP)     : %4d %2d %2d %2d %2d %2d \n',...
	p.time_end);
[strlat,strlon] = pos2str([p.poss(1)+p.poss(2)/60,p.poss(3)+p.poss(4)/60]);
fprintf(fid,'Start position       : %s   ',strlat);
fprintf(fid,'%s \n',strlon);
[strlat,strlon] = pos2str([p.pose(1)+p.pose(2)/60,p.pose(3)+p.pose(4)/60]);
fprintf(fid,'End position         : %s   ',strlat);
fprintf(fid,'%s \n',strlon);
fprintf(fid,'Magnetic deviation   : %f \n',p.drot);
if isfield(p,'timoff')
  fprintf(fid,'Time offset          : %f \n',p.timoff);
else
  fprintf(fid,'Time offset          : 0.0 \n',[]);
end
fprintf(fid,'Depths               : %d %d %d \n',p.zpar);
fprintf(fid,' \n',[]);
fprintf(fid,'Processing thresholds  \n',[]);
fprintf(fid,' \n',[]);
fprintf(fid,'min percent good     : %f \n',p.pglim);
fprintf(fid,'max error velocity   : %f \n',p.elim);
fprintf(fid,'max w variability    : %f \n',p.wlim);
fprintf(fid,' \n',[]);
fprintf(fid,'Processing results  \n',[]);
fprintf(fid,' \n',[]);
fprintf(fid,'u ship speed   (m/s) : %f \n',p.uship);
fprintf(fid,'v ship speed   (m/s) : %f \n',p.vship);
fprintf(fid,'x displacement   (m) : %f \n',p.xdisp);
fprintf(fid,'y displacement   (m) : %f \n',p.ydisp);
fprintf(fid,'bottom found at  (m) : %f \n',p.zbottom);
fprintf(fid,'maximum instr dep(m) : %f \n',p.maxdepth);
fprintf(fid,'maximum prof dep (m) : %f \n',max(dr.z));
fprintf(fid,'bottom track mode    : %f \n',p.btrk_mode);
fprintf(fid,'bottom track thrsd   : %f \n',p.btrk_ts);
fprintf(fid,'bottom track range   : %f \n',p.btrk_range);
if p.btrk_used == 1
  fprintf(fid,'bottom track source  : RDI \n',[]);
elseif p.btrk_used == 2
  fprintf(fid,'bottom track source  : normal data \n',[]);
elseif p.btrk_used == -1
  fprintf(fid,'bottom track source  : not enough data \n',[]);
elseif p.btrk_used == 0
  fprintf(fid,'bottom track source  : none \n',[]);
end
fprintf(fid,'profile ensembles    : %d to %d \n',d.firstlastindx);
fprintf(fid,' \n',[]);
fprintf(fid,'Computer : %s \n',computer);
fprintf(fid,'Total processing time: %4.0f seconds \n',toc);

fclose(fid);
