%======================================================================
%                    S A V E A R C H . M 
%                    doc: Wed Jan  7 16:51:58 2009
%                    dlm: Thu Aug 15 12:17:03 2013
%                    (c) 2009 A.M. Thurnherr
%                    uE-Info: 217 16 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jan  7, 2009: - tightened use of exist()
%   Aug 15, 2013: - adapted to new [ladcp2cdf.m]

function [da]=savearch(dr,d,p,ps,f,att)
% function [da]=savearch(dr,d,p,ps,f,att)
%
% store LADCP result LADCP archive format
%
p=setdefv(p,'ladcp_station',NaN);
p=setdefv(p,'ladcp_cast',1);
g=gregoria(d.time_jul(1));
p=setdefv(p,'ref_year',g(1));
year0=julian([p.ref_year,0,0,0,0,0]);


 da.GEN_Velocity_Units                = 'm/s';
 da.GEN_LADCP_station                 = p.ladcp_station;
 da.GEN_LADCP_cast                    = p.ladcp_cast;
 da.GEN_Profile_start_decimal_day     = d.time_jul(1)-year0;

[m,ii]=min(d.z);
 da.GEN_Profile_bottom_decimal_day    = d.time_jul(ii)-year0; 
 da.GEN_Profile_end_decimal_day       = d.time_jul(end)-year0; 

 da.GEN_Profile_start_longitude       = p.poss(3)+p.poss(4)/60; 
 da.GEN_Profile_end_longitude         = p.pose(3)+p.pose(4)/60; 
 da.GEN_Profile_start_latitude        = p.poss(1)+p.poss(2)/60; 
 da.GEN_Profile_end_latitude          = p.pose(1)+p.pose(2)/60; 

 da.GEN_Ocean_depth_m                 = round(p.zbottom);
 da.GEN_Profile_max_depth_m           = round(p.maxdepth);
 da.GEN_Magnetic_deviation_deg        = p.drot;

 da.BAR_ref_U                     = dr.ubar; 
 da.BAR_ref_V                     = dr.vbar;
 da.BAR_ref_error                 = 2*p.nav_error/p.dt_profile; 
 da.BAR_tide_U                    = NaN;
 da.BAR_tide_V                    = NaN;
 da.INPUT_SADCP_profile_avail           = existf(dr,'u_sadcp'); 
 da.INPUT_Pegasus_profile_avail         = 0; 
 da.INPUT_Bottom_track_profile_avail    = (p.btrk_used>0 & existf(dr,'zbot'));
 da.INPUT_Nav_time_series_avail         = length(f.nav)>1 ; 
 da.INPUT_CTD_time_series_avail         = length(f.ctd)>1 ; 

 da.BAR_ref_descr                 = [];
 if ps.barofac>0, da.BAR_ref_descr= [da.BAR_ref_descr,'[Navigation]']; end
 if ps.botfac>0 & existf(dr,'zbot'), 
        da.BAR_ref_descr= [da.BAR_ref_descr,'[Bottom-track]']; 
 end
 if ps.dragfac>0, da.BAR_ref_descr= [da.BAR_ref_descr,'[CTD-drag]']; end
 if ps.sadcpfac>0, da.BAR_ref_descr= [da.BAR_ref_descr,'[SADCP]']; end
 if length(da.BAR_ref_descr) < 1; da.BAR_ref_descr='[NA]'; end
 da.BAR_tide_model=                 '[NA]';

% down instrument information

 da.LADCP_dn_hard_type=['[RDI-',int2str(d.down.Frequency),'BB]'];
 if round(d.down.Frequency)==300
  da.LADCP_dn_hard_type=['[RDI-',int2str(d.down.Frequency),'WH]'];
 end
 if existf(d.down,'NarrowBand')
  da.LADCP_dn_hard_type=['[RDI-',int2str(d.down.Frequency),'NB]'];
 end
 da.LADCP_dn_hard_freq_kHz            = d.down.Frequency;
 da.LADCP_dn_hard_SN                  = p.down_sn;
 da.LADCP_dn_hard_FV                  = d.down.Firm_Version;
 da.LADCP_dn_hard_TNO                 = '[convex4]';
 da.LADCP_dn_hard_beam_ang_deg        = d.down.Beam_angle;
 da.LADCP_dn_hard_comp_type           = '[RDI]';
 da.LADCP_dn_hard_general_comments    = '  ';   

 da.LADCP_dn_conf_blank_intvl_m       = d.down.Blank/100;
 da.LADCP_dn_conf_bin_len_m           = d.down.Cell_length/100; 
 da.LADCP_dn_conf_pulse_len_m         = d.down.Pulse_length/100;
 da.LADCP_dn_conf_number_bins         = length(d.izd); 
 da.LADCP_dn_conf_ping_stagr          = '[NA]';
 da.LADCP_dn_conf_ping_trns_intvl_sec = d.down.Time_Pings;
 da.LADCP_dn_conf_number_pings        = d.down.Pings_per_Ensemble;
 da.LADCP_dn_conf_vel_ambiguity       = p.ambiguity;
 da.LADCP_dn_conf_single_ping_acc     = d.down.Single_Ping_Err;
 da.LADCP_dn_xmit_cur            = p.xmc(1);
 da.LADCP_dn_xmit_vol            = p.xmv(1);
 da.LADCP_dn_xmit_pings          = p.nping_total(1);
 da.LADCP_dn_beam_range          = p.dn_range;
 if d.down.Coordinates==3
  da.LADCP_dn_conf_coord_system='[earth]';
 elseif d.down.Coordinates==1
  da.LADCP_dn_conf_coord_system='[beam]';
 else
  da.LADCP_dn_conf_coord_system='[unknown]';
 end
 da.LADCP_dn_conf_bottom_trkr=p.btrk_used;
 if isfinite(p.zbottom) & p.btrk_used>0 & d.down.Up==0
  da.LADCP_dn_btrk_u_bias = p.btrk_u_bias;
  da.LADCP_dn_btrk_v_bias = p.btrk_v_bias;
  da.LADCP_dn_btrk_u_std =  p.btrk_u_std;
  da.LADCP_dn_btrk_v_std =  p.btrk_v_std;
 end
 da.LADCP_dn_conf_general_comments='   ';

if length(f.ladcpup)>1
% up instrument information

 da.LADCP_up_hard_type=['[RDI-',int2str(d.up.Frequency),'BB]'];
 if round(d.up.Frequency)==300
  da.LADCP_up_hard_type=['[RDI-',int2str(d.up.Frequency),'WH]'];
 end
 da.LADCP_up_hard_freq_kHz            = d.up.Frequency;
 da.LADCP_up_hard_SN                  = p.up_sn;
 da.LADCP_up_hard_FV                  = d.up.Firm_Version;
 da.LADCP_up_hard_TNO                 = '[convex4]';
 da.LADCP_up_hard_beam_ang_deg        = d.up.Beam_angle;
 da.LADCP_up_hard_comp_type           = '[RDI]';
 da.LADCP_up_hard_general_comments    ='  ';   

 da.LADCP_up_conf_blank_intvl_m       = d.up.Blank/100;
 da.LADCP_up_conf_bin_len_m           = d.up.Cell_length/100; 
 da.LADCP_up_conf_pulse_len_m         = d.up.Pulse_length/100;
 da.LADCP_up_conf_number_bins         = length(d.izd); 
 da.LADCP_up_conf_ping_stagr          = '[NA]';
 da.LADCP_up_conf_ping_trns_intvl_sec = d.up.Time_Pings;
 da.LADCP_up_conf_number_pings        = d.up.Pings_per_Ensemble;
 da.LADCP_up_conf_vel_ambiguity       = p.ambiguity;
 da.LADCP_up_conf_single_ping_acc     = d.up.Single_Ping_Err;
 da.LADCP_up_xmit_cur                 = p.xmc(2);
 da.LADCP_up_xmit_vol                 = p.xmv(2);
 da.LADCP_up_xmit_pings               = p.nping_total(2);
 da.LADCP_up_beam_range               = p.up_range;

 if p.rotup2down==2
  da.LADCP_up_compass='[velocity-match]';
 elseif p.rotup2down==1
  da.LADCP_up_compass='[down-compass]';
 else
  da.LADCP_up_compass='[up-compass]';
 end
 if d.up.Coordinates==3
  da.LADCP_up_conf_coord_system='[earth]';
 elseif d.up.Coordinates==1
  da.LADCP_up_conf_coord_system='[beam]';
 else
  da.LADCP_up_conf_coord_system='[unknown]';
 end
 da.LADCP_up_conf_general_comments='   ';
end

 da.GEN_LADCP_ensemble_time_mean_sec=mean(diff(d.time_jul*24*3600));
 da.GEN_LADCP_ensemble_time_std_sec=std(diff(d.time_jul*24*3600));
 da.GEN_conf_general_comments ='  ';

 da.GEN_Matlab_version=version;
 da.GEN_Processing_personnel= p.whoami;
 da.GEN_Processing_date=date;
 da.GEN_Proc_methodology= '[inverse]';
 da.GEN_Software_orig= p.software;
 % da.Sound_sp_calc= (choose one) <T> <T-P> <T-P-S> <NA> <unconfirmed> 
 if d.soundc==1
  if existf(d,'ctdprof_ss') | existf(d,'ctd_ss')
   da.GEN_Sound_sp_calc= '[CTD]';
  else
   da.GEN_Sound_sp_calc= '[T-P]';
  end
 else
  da.GEN_Sound_sp_calc= '[NA]';
 end
 % da.Depth_source:  (choose one) <w> <w&Pmax> <measured P (CTD)>
 %                <measured P (other)> <NA> <unconfirmed>
 da.GEN_Depth_source=  '[w]'; 
 if p.ladcpdepth==2
  da.GEN_Depth_source=  '[w&surface&bottom]';
 end
 if isfinite(p.zpar(2))
  da.GEN_Depth_source(end)=[];
  da.GEN_Depth_source=  [da.GEN_Depth_source,'&Pmax]']; 
 end
 if p.ctddepth==1
  da.GEN_Depth_source=  '[measured P (CTD)]';
 end

jok = cumprod(size(find(~isnan(d.rw))));
j = cumprod(size(find(isnan(d.re) & ~isnan(d.rw))));
 da.GEN_Percent_3beam=round(j*100/jok); 
 da.GEN_Editing_parm_descr=p.outlier;
 da.GEN_Inverse_weight_bottom=ps.botfac;
 da.GEN_Inverse_weight_navigation=ps.barofac;
 da.GEN_Inverse_weight_smooth=ps.smoofac;
 da.GEN_Proc_general_comments='  ';


if exist([f.res,'.log'],'file')
 diary off
 id=fopen([f.res,'.log']);
 da.LOG_Inverse_log=setstr(fread(id))';
 fclose(id);
end

% save to matlab file 
if p.savemat
 disp(['save ',[f.res,'.ladcp.mat'],' da dr att'])
 eval(['save ',[f.res,'.ladcp.mat'],' da dr att'])
end

% save to netcdf file if you have a modern-enough version of matlab

if p.savecdf
  if exist('ncwrite') == 2
   dr.tim=dr.tim-year0;
   disp([' save results in netcdf file: ',f.res,'.nc'])
   ladcp2cdf([f.res,'.nc'],dr,da,p,ps,f,att)
  else
   disp('upgrade your matlab version to get the netcdf LADCP archive output')
  end 
end
