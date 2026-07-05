function [] = mkSADCP(CODAS_CONTOUR_DIR,SADCP_OUTPUT_FILE)
% mkSADCP(CODAS_CONTOUR_DIR,SADCP_OUTPUT_FILE)
%
%	  CODAS_CONTOUR_DIR path with CODAS output files 
%			    [contour_xy.mat] and [contour_uv.mat]
% 	  SADCP_OUTPUT_FILE output file (including path)

%======================================================================
%                    M K S A D C P . M 
%                    doc: Mon Mar 15 06:00:43 2004
%                    dlm: Fri Jan  6 12:28:45 2017
%                    (c) 2004 A.M. Thurnherr
%                    uE-Info: 31 61 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================
%
% pre-processing script to get UH/CODAS-processed SADCP data (usually
% available in near-real time) into a format that can be loaded with
% [loadsadcp.m].
%
% NOTES:
%     - this script needs only be called after re-processing the SADCP
%	data with the UH/CODAS software
%     - during ship-board processing, this script will be called near the 
%	beginning of [set_cast_params.m] to update the SADCP matrix
%     - during post-cruise processing this script will only have to be
%	called once, when the finalized SADCP data have become available
%
% HISTORY:
%  Mar 15, 2004: - created
%  Jan  5, 2007: - documented and added to IX_4 distribution
%  Jan  6, 2016: - added time of first SADCP record to output

%======================================================================

eval(['load ' CODAS_CONTOUR_DIR '/contour_xy']);
eval(['load ' CODAS_CONTOUR_DIR '/contour_uv']);

u_sadcp   = uv(:,1:2:end);
v_sadcp   = uv(:,2:2:end);
lat_sadcp = xyt(2,:);
lon_sadcp = xyt(1,:)-360;
dday      = xyt(3,:);
z_sadcp   = zc;

badi             = find(isnan(dday) | isnan(lat_sadcp) | isnan(lon_sadcp));
dday(badi)       = [];
lat_sadcp(badi)  = [];
lon_sadcp(badi)  = [];
u_sadcp(:, badi) = [];
v_sadcp(:, badi) = [];
tim_sadcp        = julian([year_base 1 1 0 0 0])+dday;

lat_sadcp = lat_sadcp(:);
lon_sadcp = lon_sadcp(:);
ii=find(lon_sadcp<-180);
lon_sadcp(ii)=lon_sadcp(ii)+360;
z_sadcp   = z_sadcp(:);
tim_sadcp = tim_sadcp(:);

disp(sprintf('first SADCP processing date %d-%d-%d %d:%d:%d',...
	gregoria(tim_sadcp(1))));
disp(sprintf(' last SADCP processing date %d-%d-%d %d:%d:%d',...
	gregoria(tim_sadcp(end))));
disp([' Time is now: ',datestr(now)])

eval(['save ' SADCP_OUTPUT_FILE ' tim_sadcp lat_sadcp lon_sadcp u_sadcp v_sadcp z_sadcp']);
