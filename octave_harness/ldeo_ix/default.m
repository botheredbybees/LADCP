%======================================================================
%                    D E F A U L T . M 
%                    doc: Sat Jun 26 06:10:09 2004
%                    dlm: Fri Sep 14 18:19:40 2018
%                    (c) 2004 ladcp@
%                    uE-Info: 40 50 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%	Jan  7, 2009: - tightened use of exist()
%	Jul 22, 2009: - typo in ps.smallfac
%				  - changed default of ps.smallfac
%	Apr 26, 2012: - removed diffusivity calculation
%	May  4, 2012: - updated version to IX_8
%	Jun 18, 2013: - updated version to IX_9
%	Jun 24, 2013: - updated version to IX_10beta
%	Mar  5, 2014: - updated version to IX_10
%	Mar 21, 2014  - updated version to IX_11beta
%	Jun 11, 2014: - changed savecdf default to 0
%	Nov 25, 2015: - changed version to IX_11
%	Feb 18, 2016: - changed version to IX_12 and published
%				  - changed version to IX_13beta
%	Mar 29, 2017: - added saveplot_pdf
%	Jan 17, 2018: - changed ersion to IX_13 and published

% LADCP processing software
% M. Visbeck. LDEO/2003
% http://www.ldeo.columbia.edu/ladcp
%
% set default values for parameter in LADCP processing
% 
%  INPUT is provided by three structures.
%  some of this may get changed by the software
%
% structure f.??? contains file names
% structure p.??? contains parameter relevant to reading and preparing
%             the data
% structure ps.??? contains parameter for the solution
% structure att.??? contains attributes
p.software='LDEO LADCP software: Version IX_14beta';

% file names
% f.ladcpdo  is the ONLY required input
% f=setdefv(f,'ladcpdo',' ');
%
% you can try also our fuzzy file finder
% f=getfile('/data/kN001',f);
% which searches for "BB-ADCP" file names and picks the one with the most data 


% file for up looking ADCP data
f=setdefv(f,'ladcpup',' ');

% file for CTD time series data
f=setdefv(f,'ctd',' ');

% file for CTD profile data
f=setdefv(f,'ctdprof',' ');

% file for ships navigation time series
f=setdefv(f,'nav',' ');

% file for SADCP velocity data
f=setdefv(f,'sadcp',' ');

% file name for results (extensions will be added by software)
%  *.bot            bottom referenced ASCII data
%  *.lad            profile ASCII data
%  *.mat            MATLAB  format >> dr p ps f
%  *.cdf            NETCDF  (binary) LADCP data format 
%  *.log            ASCII log file of processing
%  *.txt            ASCII short log file
%  *.ps             post-script figure of result 
% default no file output
f=setdefv(f,'res',' ');

%--------------------------------------------------------------
% Parameter for loading and primary error check  p.* structure
%--------------------------------------------------------------


% p. structure has a lot of important information, some given, some computed

% set to 1 or higher if you want the processing software to pause frequently
p=setdefv(p,'debug',0);

% name of data processer
% on LINUX computer can get login ID otherwise set to your name
% p=setdefv(p,'whoami',' Joe Blow');
p=setdefv(p,'whoami',whoami);

% give serial number of up and down instrument
p=setdefv(p,'down_sn',NaN);
p=setdefv(p,'up_sn',NaN);

% STATION name
% p=setdefv(p,'name','unknown');

% give LADCP station number
p=setdefv(p,'ladcp_station',NaN);

% give LADCP cast number
p=setdefv(p,'ladcp_cast',1);   

% give Cruise id
p=setdefv(p,'cruise_id', 'unknown');

% TIME range
% take from data if not given
%  [year, month, day, hour, minute, second]
%  i.e  [1997 8 6 9 22 45]
% p=setdefv(p,'time_start',[0 0 0 0 0 0]);
% p=setdefv(p,'time_end',[0 0 0 0 0 0]);

% how to get depth from W integration
% getdepth=1 use plain integral of W (mfile getdpth)
% getdepth=2 use inverse method to use bottom reflection
% and integral of W (mfile getdpthi) [default]
%	NB (IX_6): based on a single profile that I've processed,
% 			   as well as on comments by Gerd & Martin,
%			   it may well be that getdepth=1 works better
%			   with shallow stations
 p =setdefv(p,'getdepth',2);

% restrict time range to profile and disregard data close to surface
% p.cut = 0 dont restrict
% p.cut > 0 restrict time to adcp depth below a depth of p.cut
% p=setdefv(p,'cut',10);


% POSITION of the start and end point
% [degree lat, minute lat, degree lon, minute lon]
%  i.e. [-59 -30.5697 -44 -22.4986]
% p=setdefv(p,'pose',[0 0 0 0]);
% p=setdefv(p,'poss',[0 0 0 0]);
% navigation error in m
% p=setdefv(p,'nav_error',30);
% average navigation from nav file over a certain fration of days (2 minutes)
% p=setdefv(p,'navtime_av',2/60/24);

% SUPER ENSEMBLES 
% 	are calculated in prepinv.m to reduce the number of raw profiles
% 	The ides is to obtain one average profile for each vertical dz=const
% 	that the CTD traveled through. As a result a constant number of super
% 	ensembles are obtained for the up and down cast.
%   but a fixed number of ensembles can also be averaged

% p.avdz sets the depth interval between adjacent super-ensembles
% default one bin length
% p=setdefv(p,'avdz',medianan(abs(diff(d.izm(:,1)))));

% p.avens overrided p.avdz and sets a fixed number of ensembles to average
% default NAN not used
%	NB (IX_6): When p.avens == 1, p.single_ping_accuracy has to be set!
%			   Otherwise, the software cannot determine the weight of the
%			   BT constraint.
% p=setdefv(p,'avens',NaN);

% p.oversample increases the number of individual profiles used for
% each super ensemble by that factor (not recomended)
% p=setdefv(p,'oversample',1.0);

% BOTTOM TRACK
% 	The are several options to get bottom track data
% 
% mode = 1 :   use only RDI bottom track
%        2 :   use only own bottom track
%        3 :   use RDI, if existent, own else (default)
%        0 :   use not bottom track at all
% p=setdefv(p,'btrk_mode',3);

% p.btrk_ts is in dB to detect bottom above bin1 level (for own btm track)
% p=setdefv(p,'btrk_ts',10);

% p.btrk_below gives binoffset used below target strength maximum
% to make bottom track velocity
% p=setdefv(p,'btrk_below',1);

% p.btrk_range gives minumum / maximum distance for bottom track
% p=setdefv(p,'btrk_range',[300 50]);

% p.btrk_wstd gives maximum accepted wstd for super ensemble averages
% p=setdefv(p,'btrk_wstd',0.1);

% maximum allowed difference between reference layer W and W bottom track
% p=setdefv(p,'btrk_wlim',0.05);

% force to recalculate bottom distance using target strenght
p=setdefv(p,'bottomdist',0);

% p.surfdist = 1 use surface reflections of up looking ADCP to get start
% depth
p=setdefv(p,'surfdist',1);

% MAGNETIC deviation in degree
p=setdefv(p,'drot',NaN);

% COMPASS manipulation
% experts only 
% fix_compass:1 means hdg_offset gets added
% fix_compass:2 means up looker gets down compass + hdg_offset
% fix_compass:3 means down looker gets up compass + hdg_offset
% 
p=setdefv(p,'fix_compass',0);

% give compass offset in addition to declination (1) for down (2) for up
% p=setdefv(p,'hdg_offset',[0 0]);

% COMPASS: 
% how to best adjust compass to best match 
% if 1 rotate up-looking and down-looking instrument to mean heading
%    2 rotate up-looking and down-looking velocities to match up velocities
%        (not really recommended)
%    3 rotate up-looking velocities to down heading
%        (use if suspect the up heading is bad
%    4 rotate down-looking velocities to up heading
%        (use if suspect the down heading is bad
p=setdefv(p,'rotup2down',1);

% OFFset correction
% if 1 remove velocity offset between up and down looking ADCP
% this will correct errors due to tilt biases etc.
p=setdefv(p,'offsetup2down',1);


% DEPTH of the start, bottom and end of the profile
p=setdefv(p,'zpar',[0 NaN 0]);

% blank data below the  BOTTOM  
% default use first 2 then one bin length
% p=setdefv(p,'dzbelow',[2 1]*medianan(abs(diff(d.zd))));

% maximum number of bins to be used
% 0 mean all get used 
% p=setdefv(p,'maxbinrange',0);

% set ctdtime=1 if CTD data are to be used to limit range of good data
% p=setdefv(p,'ctdtime',1);

% set ctdmaxlag.=100 to the maximu pings that the ADCP data can be shifted to 
% best match W calculated from CTD pressure time series (loadctd)
% If you have good times set it to 10... if your time base is questionable
% you can set it 100 or more
% p=setdefv(p,'ctdmaxlag',100);

% minimum correlation to be accepted for shifting time of ADCP using 
% W calculated from CTD pressure and W from ADCP
% p=setdefv(p,'ctdmincorr',0.9-90/length(wctd));
% p=setdefv(p,'ctdmincorr',0.8);

% produce much more RAW data output
% set to 1 for more data
% p=setdefv(p,'orig',0);
 
% save individual target strength p.ts_save=[1 2 3 4]
% p=setdefv(p,'ts_save',0);
% save individual correlation p.cm_save=[1 2 3 4]
% p=setdefv(p,'cm_save',0);
% save individual percent good pings p.pg_save=[1 2 3 4]
% p=setdefv(p,'pg_save',0);


%OUTLIER detection is called twice once to clean the raw data
%	and a second time to clean the super ensembles
%        [n1 n2 n3 ...] the length gives the number of scans and
%	each value the maximum allowed departure from the mean in std
%	applied for the u,v,w fields for each bin over blocks 
%   of p.outlier_n profiles
% p=setdefv(p,'outlier',[4.0  3.0]);
% default for p.outlier_n is number of profiles in 5 minutes
% p=setdefv(p,'outlier_n',100);
% minimum std for horizontal velocities of super ensemble
% p=setdefv(p,'superens_std_min',0.01);


%SPIKES
% 	maximum value for abs(V-error) velocity
% p=setdefv(p,'elim',0.5);
% 	maximum value for horizontal velocity 
% p=setdefv(p,'vlim',2.5);
% 	minimum value for %-good
% p=setdefv(p,'pglim',0); 
%	maximum value for W difference between the mean W and actual
%        W(z) for each profile. 
% p=setdefv(p,'wlim',0.20);

% TILT  flag data with large tilt or tilt differences as bad
% [22  (max tilt allowed) 
%  4 (maximum tilt difference between pings allowd)]
% p=setdefv(p,'tiltmax',[22 4]);

% TILT  reduce weight for large tilts
% p=setdefv(p,'tilt_weight',10]);

% fix TIME of the ADCP in days
% p=setdefv(p,'timoff',0);

% usually bin 1 is not very good thus you can reduce its weight
% p=setdefv(p,'weighbin1',0.5);

% apply sound speed correction
% p=setdefv(p,'soundcorr',1);

% apply tilt correction
% tiltcor(1)=down-pitch bias
% tiltcor(2)=down-rol bias
% tiltcor(3)=up-pitch bias
% tiltcor(4)=up-rol bias
 p=setdefv(p,'tiltcor',0);

% Give bin number for the best W to compute depth of the ADCP
%	default uses bin 2-3 but be careful when up/down instruments
%	are used. The good bins are in the middle! 
% p=setdefv(p,'wizr',[2 3]);

% SET ambiguity velocity used [m/s]
  p = setdefv(p,'ambiguity',2.5);

% Give single ping accuracy;
%	NB: this is only used in [prepinv.m] to set the stddevs of the BT
%	    velocities if no super-ensemble-averaging is carried out
  p = setdefv(p,'single_ping_accuracy',NaN);

% Give warning when transmit current drops below minimum
  p = setdefv(p,'xmv_min',0);

% Write matlab file
  p = setdefv(p,'savemat',0);

% Write netcdf file
  p = setdefv(p,'savecdf',0);

% Collect warnings with regards to LADCP processing
  p.warnp=' LADCP processing warnings: ';

% Save Plots 
% Save figure numbers to ps file
%    1 : Summary Plot
%    2 : Engineering Data
%    3 : Data Quality
%    4 : Depth
%    5 : Heading Corrections 
%    6 : Up/Down Differences
%    7 : CTD Position
%    8 : Shear
%    9 : SADCP U, V 
%   10 : U, V Offsets, Tilt Error
%   11 : Processing Warnings
%   12 : Inversion Constraints
%   13 : Bottom Track detail
%	14 : Echo Amplitude (data editing)

p = setdefv(p,'saveplot',1);
p = setdefv(p,'saveplot_png',[]);
p = setdefv(p,'saveplot_pdf',[]);
  
%--------------------------------------------------------------
% Parameter for inversion   ps.* structure
%--------------------------------------------------------------

% Process data using shear based method
% compute shear based solution
% ps.shear=2  ; use super ensemble
% ps.shear=1  ; use raw data

if exist('ps','var')
  ps=setdefv(ps,'shear',1);
else
  ps.shear=1;
end

% decide how to weight data 
% 1 : use super ensemble std 
% 0 : use correlation based field
ps=setdefv(ps,'std_weight',1);

% Weight for the barotropic constraint
% ps=setdefv(ps,'barofac',1);

% Weight for the bottom track constraint
% ps=setdefv(ps,'botfac',1); 

% Process up and down cast seperately
% ps=setdefv(ps,'down_up',1);

% Depth resolution for final profile
%	default one bin length
% ps=setdefv(ps,'dz',medianan(abs(diff(di.izm(:,1)))));

% Smoothing of the final profile
% ps=setdefv(ps,'smoofac',0.01);

% comment this out to request that shears are small  (experts only)
ps=setdefv(ps,'smallfac',[1 0]);

% weight bottom track data with distance of bottom
%  use Gaussian with offset (btrk_weight_nblen(1) * bin)
%  and width (btrk_weight_nblen(2) * bin) 
%  one might set this to [15 5] to reduce the weight of close by bottom track data
ps=setdefv(ps,'btrk_weight_nblen',[0 0]);

% Weight for SADCP data
% ps.sadcpfac=1 about equal weight for SDACP profile
% ps=setdefv(ps,'sadcpfac',3);

% average over data within how many standart deviations
% ps=setdefv(ps,'shear_stdf',2);

% the minimum weight a bin must have to be accepted for shear
% ps=setdefv(ps,'shear_weightmin',0.1);

% restrict inversion to one instrument only 1: up+dn, 2:dn only  3:up only
ps=setdefv(ps,'up_dn_looker',1);

% super ensemble velocity error
% try to use the scatter in W to get an idea of the "noise"
% in the velocity measurement
%% This is a bit of code used in GETINV.m
%% nmax=min(length(di.izd),7);
%% sw=stdnan(di.rw(di.izd(1:nmax),:)); ii=find(sw>0);
%% sw=medianan(sw(ii))/tan(p.beamangle*pi/180);
%% ps=setdefv(ps,'velerr',max([sw,0.02]));
%
% ps=setdefv(ps,'velerr',0.02);


% How to solve the inverse
%     ps.solve = 0  Cholseky transform
%              = 1  Moore Penrose Inverse give error for solution
% ps=setdefv(ps,'solve',1); 

% Threshold for minimum weight, data with smaller weights
%  	will be ignored
% ps=setdefv(ps,'weightmin',0.05);

% Change the weights by 
%	weight=weight^ps.weightpower 
% ps=setdefv(ps,'weightpower',1); 

% Change remove 1% of outlier after solve
% ps.outlier times
 ps=setdefv(ps,'outlier',1); 

% set ps.down_up=1 if up/down cast should be solved seperately
% ps=setdefv(ps,'down_up',1);setdefv(ps,'weightpower',1);

% Weight for the cable drag constraint
% only for experts
% ps=setdefv(ps,'dragfac',0); 
% ps=setdefv(ps,'drag_tilt_vel',15);
% ps=setdefv(ps,'drag_lagmax',15);
% ps=setdefv(ps,'drag_zmax',2000);

% Set fixed range for velocity plots
% ps=setdefv(ps,'urange',ur);
% ps=setdefv(ps,'zrange',ax(3:4));

% ------------------------------------------------------------------
% OUTPUT data files
% dr structure has the main output from the LADCP processing
%
% name: 'demo-4'                   | file name
% date: [1997 8 5 12 8 54]         | mean time 
% lat: -58.7044                    | mean position N
% lon: -44.5215                    | mean position E
% zbot: [21x1 double]              | bottom referenced profile depth [m]
% ubot: [21x1 double]              | bottom referenced profile U [m/s]
% vbot: [21x1 double]              | bottom referenced profile V [m/s]
% uerrbot: [21x1 double]           | bottom referenced velocity error[m/s]
% z_sadcp: [37x1 double]           | SADCP profile depth [m]
% u_sadcp: [37x1 double]           | SADCP profile U [m/s]
% v_sadcp: [37x1 double]           | SADCP profile V [m/s]
% uerr_sadcp: [37x1 double]        | SADCP profile error [m/s]
% z: [153x1 double]                | LADCP profile depth [m]
% u: [153x1 double]                | LADCP profile U [m/s]
% v: [153x1 double]                | LADCP profile V [m/s]
% uerr: [153x1 double]             | LADCP profile erroe [m/s]
% nvel: [153x1 double]             | LADCP number of ensembles per bin
% ubar: -0.0974                    | LADCP U barotropic [m/s]
% vbar: -0.0041                    | LADCP V barotropic [m/s]
% tim: [1x246 double]              | Station time series [Julian Days]
% shiplon: [1x246 double]          | Ships position time series E
% shiplat: [1x246 double]          | Ships position time series N
% xship: [1x246 double]            | Ships position relative to start E [m]
% yship: [1x246 double]            | Ships position relative to start N [m]
% uship: [1x246 double]            | Ships velocity U [m/s]
% vship: [1x246 double]            | Ships velocity V [m/s]
% zctd: [1x246 double]             | Depth of CTD [m]
% wctd: [1x246 double]             | CTD velocity W [m/s]
% uctd: [1x246 double]             | CTD velocity U [m/s]
% vctd: [1x246 double]             | CTD velocity V [m/s]
% uctderr: [246x1 double]          | CTD velocity error [m/s]
% xctd: [1x246 double]             | CTD position relative to start E [m]
% yctd: [1x246 double]             | CTD position relative to start E [m]
% range: [153x1 double]            | ADCP total range of data [m]
% range_do: [153x1 double]         | ADCP down looking range of data [m]
% range_up: [153x1 double]         | ADCP up looking range of data [m]
% ts: [1x153 double]               | ADCP echo amplitude profile bin 1
% ts_out: [1x153 double]           | ADCP echo amplitude profile last down bin 
% u_do: [153x1 double]             | LADCP down only profile U [m/s]
% v_do: [153x1 double]             | LADCP down only profile V [m/s]
% u_up: [153x1 double]             | LADCP up only profile U [m/s]
% v_up: [153x1 double]             | LADCP up only profile V [m/s]
% p: [153x1 double]                | LADCP profile pressure [dBar]
% ensemble_vel_err: [153x1 double] | ADCP ensemble velocity error [m/s]
% u_shear_method: [153x1 double]   | LADCP shear method profile U [m/s]
% v_shear_method: [153x1 double]   | LADCP shear method profile V [m/s]
% ctd_t:          [153x1 double]   | CTD profile temperature [^oC]
% ctd_s:          [153x1 double]   | CTD profile salinity
% ctd_ss:         [153x1 double]   | CTD profile sound speed [m/s] 
% ctd_N2:         [153x1 double]   | CTD profile stability [1/s^2]

% 
% der structure contains error analysis
% de  structure contains inverse matricies.
% 

% Attribute List 
% For variable dr.var1, we will have att.var1.long_name ... etc

% Attributes For Structure dr
att.name.long_name      	  = 'Cast ID';
att.date.long_name      	  = 'Date';
att.date.units          	  = 'Y M D H M S';
att.lat.long_name       	  = 'Latitude';
att.lat.units           	  = 'Degree North';
att.lon.long_name       	  = 'Longitude';
att.lon.units           	  = 'Degree East';
att.zbot.long_name      	  = 'Bottom Referenced Profile Depth';
att.zbot.units          	  = 'm';
att.ubot.long_name      	  = 'Bottom Referenced Profile U';
att.ubot.units          	  = 'm/s';
att.vbot.long_name      	  = 'Bottom Referenced Profile V';
att.vbot.units          	  = 'm/s';
att.uerrbot.long_name   	  = 'Bottom Referenced Profile Velocity Error';
att.uerrbot.units       	  = 'm/s';
att.z_sadcp.long_name   	  = 'SADCP Profile Depth';
att.z_sadcp.units   		  = 'm';
att.u_sadcp.long_name   	  = 'SADCP Profile U';
att.u_sadcp.units   		  = 'm/s';
att.v_sadcp.long_name   	  = 'SADCP Profile V';
att.v_sadcp.units   		  = 'm/s';
att.uerr_sadcp.long_name   	  = 'SADCP Profile Velocity Error';
att.uerr_sadcp.units   		  = 'm/s';
att.z.long_name                   = 'Depth';
att.z.units                       = 'Meters';
att.u.long_name                   = 'U';
att.u.units                       = 'm/s';
att.v.long_name                   = 'V';
att.v.units                       = 'm/s';
att.uerr.long_name                = 'Velocity Error';
att.uerr.units                    = 'm/s';
att.nvel.long_name                = 'LADCP number of ensembles per bin';
att.ubar.long_name                = 'LADCP U Barotropic';
att.ubar.units              	  = 'm/s';
att.vbar.long_name                = 'LADCP V Barotropic';
att.vbar.units              	  = 'm/s';
att.tim.long_name 		  = 'Station Time Series';
att.tim.units 			  = 'Julian Days';
att.tim_hour.long_name 		  = 'Station Time Series';
att.tim_hour.units 		  = 'Hour of Day';
att.shiplon.long_name 		  = 'Longitude';
att.shiplon.units     		  = 'Degree East';
att.shiplat.long_name 		  = 'Latitude';
att.shiplat.units     		  = 'Degree North';
att.xship.long_name   		  = 'Ship Position E';
att.xship.units       		  = 'Meters East';
att.yship.long_name   		  = 'Ship Position N';
att.yship.units       		  = 'Meters North';
att.uship.long_name   		  = 'Ship Velocity U';
att.uship.units       		  = 'm/s';
att.vship.long_name   		  = 'Ship Velocity V';
att.vship.units       		  = 'm/s';
att.zctd.long_name    		  = 'Depth of CTD';
att.zctd.units        		  = 'm';
att.wctd.long_name    		  = 'CTD Velocity W';
att.wctd.units        		  = 'm/s';
att.uctd.long_name    		  = 'CTD Velocity U';
att.uctd.units        		  = 'm/s';
att.vctd.long_name    		  = 'CTD Velocity V';
att.vctd.units        		  = 'm/s';
att.uctderr.long_name 		  = 'CTD Velocity Error';
att.uctderr.units     		  = 'm/s';
att.xctd.long_name    		  = 'CTD Position Relative to Start E';
att.xctd.units        		  = 'm';
att.yctd.long_name    		  = 'CTD Position Relative to Start N';
att.yctd.units        		  = 'm';
att.range.long_name      	  = 'ADCP total range of data';
att.range.units          	  = 'm';
att.range_do.long_name   	  = 'ADCP down looking range of data';
att.range_do.units       	  = 'm';
att.range_up.long_name   	  = 'ADCP up looking range of data';
att.range_up.units       	  = 'm';
att.ts.long_name         	  = 'ADCP echo amplitude profile bin 1';
att.ts.units                      = 'dB';
att.ts_out.long_name              = ['ADCP echo amplitude ' ...
		                          'profile last down bin'];
att.ts_out.units                  = 'dB';
att.u_do.long_name       	  = 'LADCP down only profile U';
att.u_do.units           	  = 'm/s';
att.v_do.long_name       	  = 'LADCP down only profile V';
att.v_do.units           	  = 'm/s';
att.u_up.long_name       	  = 'LADCP up only profile U';
att.u_up.units           	  = 'm/s';
att.v_up.long_name       	  = 'LADCP up only profile V';
att.v_up.units           	  = 'm/s';
att.p.long_name          	  = 'Pressure';
att.p.units              	  = 'dBar';
att.ensemble_vel_err.long_name    = 'ADCP ensemble velocity error';
att.ensemble_vel_err.units        = 'm/s';
att.u_shear_method.long_name      = 'LADCP shear method profile U';
att.u_shear_method.units 	  = 'm/s';
att.v_shear_method.long_name      = 'LADCP shear method profile V';
att.v_shear_method.units 	  = 'm/s';
att.ctd_t.long_name      	  = 'CTD profile temperature';
att.ctd_t.units          	  = 'Degree C';
att.ctd_s.long_name      	  = 'CTD profile salinity';
att.ctd_s.units          	  = 'psu';
att.ctd_ss.long_name     	  = 'CTD profile sound speed';
att.ctd_ss.units         	  = 'm/s';
att.ctd_N2.long_name     	  = 'CTD profile stability';
att.ctd_N2.units         	  = '1/s^2';

