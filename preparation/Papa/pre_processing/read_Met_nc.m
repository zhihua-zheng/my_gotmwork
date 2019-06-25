%% read_Met_nc
%
% This the script used to read the meterological and subsurface observation
% data from OCSPapa mooring, the MAT file is saved as 'Met_IC_ocsp.mat'.

% by Zhihua Zheng (UW/APL), updated on Mar. 25 2019


%% General setting

clear
data_dir = '~/GDrive/UW/Research/Data/OCSP/Mooring/20070608_20190616/';

WINDname = fullfile(data_dir,'w50n145w_hr.cdf');
SSTname  = fullfile(data_dir,'sst50n145w_hr.cdf');
SSSname  = fullfile(data_dir,'sss50n145w_hr.cdf');
SSDname  = fullfile(data_dir,'ssd50n145w_hr.cdf');
RHname   = fullfile(data_dir,'rh50n145w_hr.cdf');
RAINname = fullfile(data_dir,'rain50n145w_hr.cdf');
RADname  = fullfile(data_dir,'rad50n145w_hr.cdf');
LWname   = fullfile(data_dir,'lw50n145w_hr.cdf');
BPname   = fullfile(data_dir,'bp50n145w_hr.cdf');
AIRTname = fullfile(data_dir,'airt50n145w_hr.cdf');

% profiles
Tname    = fullfile(data_dir,'t50n145w_hr.cdf');
Sname    = fullfile(data_dir,'s50n145w_hr.cdf');
Dname    = fullfile(data_dir,'d50n145w_hr.cdf');
CURname  = fullfile(data_dir,'cur50n145w_hr.cdf');

%% read Met. variables 

% -------- hours since 2007-06-08 04:00:00 --------------------------------
lat    =  ncread(WINDname,'lat');
lon    =  ncread(WINDname,'lon');
w_u    =  ncread(WINDname,'WU_422');
w_v    =  ncread(WINDname,'WV_423');
z_wind = -ncread(WINDname,'depth');
w_dir  =  ncread(WINDname,'WD_410'); % [degree]
w_gust =  ncread(WINDname,'WG_430'); % [m/s]
w_spd  =  ncread(WINDname,'WS_401');
t_w    =  ncread(WINDname,'time');

rh     =  ncread(RHname,'RH_910'); % RELATIVE HUMIDITY [%]
z_rh   = -ncread(RHname,'depth');
t_rh   =  ncread(RHname,'time');

rain   =  ncread(RAINname,'RN_485'); % RAIN RATE [mm/hr]
t_rain =  ncread(RAINname,'time');

Rl     =  ncread(LWname,'Ql_136'); % LONGWAVE RADIATION [W/m^2]
t_Rl   =  ncread(LWname,'time');

P      =  ncread(BPname,'BP_915'); % BAROMETRIC PRESSURE [hPa]
t_P    =  ncread(BPname,'time');

Ta     =  ncread(AIRTname,'AT_21'); % AIR TEMPERATURE [C]
z_Ta   = -ncread(AIRTname,'depth');
t_Ta   =  ncread(AIRTname,'time'); 

% -------- hours since 2007-06-08 00:00:00 -------------------------------- 
sst    =  ncread(SSTname,'T_25'); % [C]
t_sst  =  ncread(SSTname,'time');

sss    =  ncread(SSSname,'S_41'); % [PSU]
t_sss  =  ncread(SSSname,'time');

ssd    =  ncread(SSDname,'STH_71'); % SIGMA-THETA [kg/m^3], potden
t_ssd  =  ncread(SSDname,'time');

% -------- hours since 2007-06-08 05:00:00 --------------------------------
Rs     = ncread(RADname,'RD_495'); % DOWNGOING SHORTWAVE RADIATION [W/m^2]
t_Rs   = ncread(RADname,'time');

%% read profiles from netCDF files

% hours since 2007-06-08 00:00:00
depth_t = ncread(Tname,'depth');
tprof   = ncread(Tname,'T_20');

depth_s = ncread(Sname,'depth');
sprof   = ncread(Sname,'S_41');

depth_d = ncread(Dname,'depth');
dprof   = ncread(Dname,'STH_71');

% minutes since 2007-06-08 00:00:00
depth_cur = ncread(CURname,'depth');
uprof     = ncread(CURname,'U_320'); % [cm/s]
vprof     = ncread(CURname,'V_321'); % [cm/s]

% all the T-S-D profile data use one timestamp t_tsd
t_tsd = ncread(Tname,'time');
t_cur = ncread(CURname,'time');

%% Reduce single dimension, convert class to double

w_u    = double(squeeze(w_u));
w_v    = double(squeeze(w_v));
w_spd  = double(squeeze(w_spd));
w_dir  = double(squeeze(w_dir));
w_gust = double(squeeze(w_gust));
sss    = double(squeeze(sss));
sst    = double(squeeze(sst));
ssd    = double(squeeze(ssd));
Ta     = double(squeeze(Ta));
Rs     = double(squeeze(Rs));
Rl     = double(squeeze(Rl));
rain   = double(squeeze(rain));
rh     = double(squeeze(rh));
P      = double(squeeze(P));

lon    = double(lon);
lat    = double(lat);
z_rh   = double(z_rh);
z_wind = double(z_wind);
z_Ta   = double(z_Ta);

t_w    = double(t_w);
t_rh   = double(t_rh);
t_rain = double(t_rain);
t_Ta   = double(t_Ta);
t_P    = double(t_P);
t_Rl   = double(t_Rl);
t_Rs   = double(t_Rs);
t_ssd  = double(t_ssd);
t_sss  = double(t_sss);
t_sst  = double(t_sst);
t_tsd  = double(t_tsd);
t_cur  = double(t_cur);

% ------ absolute time ----------------------------------------------------
t_w    = datenum(2007,6,8,4,0,0) + t_w/24;
t_rh   = datenum(2007,6,8,4,0,0) + t_rh/24;
t_rain = datenum(2007,6,8,4,0,0) + t_rain/24;
t_Rl   = datenum(2007,6,8,4,0,0) + t_Rl/24;
t_P    = datenum(2007,6,8,4,0,0) + t_P/24;
t_Ta   = datenum(2007,6,8,4,0,0) + t_Ta/24;

t_sst  = datenum(2007,6,8,0,0,0) + t_sst/24;
t_sss  = datenum(2007,6,8,0,0,0) + t_sss/24;
t_ssd  = datenum(2007,6,8,0,0,0) + t_ssd/24;
t_tsd  = datenum(2007,6,8,0,0,0) + t_tsd/24;
t_cur  = datenum(2007,6,8,0,0,0) + t_cur/24/60;

t_Rs   = datenum(2007,6,8,5,0,0) + t_Rs/24;
% -------------------------------------------------------------------------

sprof  = double(squeeze(squeeze(sprof)));
tprof  = double(squeeze(squeeze(tprof)));
dprof  = double(squeeze(squeeze(dprof)));
uprof  = double(squeeze(squeeze(uprof)));
vprof  = double(squeeze(squeeze(vprof)));

depth_s   = double(depth_s);
depth_t   = double(depth_t);
depth_cur = double(depth_cur);
depth_sd  = depth_s; clear depth_s depth_d % depth_d is the same as depth_s

%% Unify the timestamp for Met. variables

% set the first timestamp: 2007-06-08 05:00:00
time   = t_Rs;          

w_u    = w_u(2:end);
w_v    = w_v(2:end);
w_spd  = w_spd(2:end);
w_dir  = w_dir(2:end);
w_gust = w_gust(2:end);   clear t_w
Ta     = Ta(2:end);       clear t_Ta
rh     = rh(2:end);       clear t_rh
P      = P(2:end);        clear t_P
Rl     = Rl(2:end);       clear t_Rl
rain   = rain(2:end);     clear t_rain

inx = find(t_sst==t_Rs(1));
sst = sst(inx:end);       clear t_sst
sss = sss(inx:end);       clear t_sss
ssd = ssd(inx:end);       clear t_ssd

tprof = tprof(:,inx:end);
sprof = sprof(:,inx:end);
dprof = dprof(:,inx:end); clear t_tsd
uprof = uprof(:,inx:end);
vprof = vprof(:,inx:end); clear t_cur

% some time series are shorter, fill with outlier
rain = [rain; ones(length(rh)-length(rain),1)*1e35];
sss  = [sss;  ones(length(rh)-length(sss), 1)*1e35];
ssd  = [ssd;  ones(length(rh)-length(ssd), 1)*1e35];

%% Save

clear data_dir *name inx t_Rs
save('~/GDrive/UW/Research/Data/OCSP/Mooring/Met_IC_ocsp.mat');
