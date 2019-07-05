
%% read_mooring_min
%
% This the script used to read the meterological and subsurface observation
% data from OCSPapa mooring, the MAT file is saved as 'ocsp_2011high.mat'.

% by Zhihua Zheng (UW/APL), updated on July 1 2019

%% General setting

clear
data_dir = '~/GDrive/UW/Research/Data/OCSP/Mooring/2011_high_res/';

WINDname = fullfile(data_dir,'w50n145w_10m.cdf');
SSTname  = fullfile(data_dir,'sst50n145w_10m.cdf');
RHname   = fullfile(data_dir,'rh50n145w_10m.cdf');
RAINname = fullfile(data_dir,'rain50n145w_10m.cdf');
RADname  = fullfile(data_dir,'rad50n145w_1m.cdf');
LWname   = fullfile(data_dir,'lw50n145w_1m.cdf');
BPname   = fullfile(data_dir,'bp50n145w_10m.cdf');
AIRTname = fullfile(data_dir,'airt50n145w_10m.cdf');

% profiles
Tname    = fullfile(data_dir,'t50n145w_10m.cdf');

%% read variables 

% -------- minutes since 2011-02-09 00:00:00 --------------------------------
lat    =  ncread(WINDname,'lat');
lon    =  ncread(WINDname,'lon');
w_u    =  ncread(WINDname,'WU_422');
w_v    =  ncread(WINDname,'WV_423');
zw     = -ncread(WINDname,'depth');
w_dir  =  ncread(WINDname,'WD_410'); % [degree]
w_gust =  ncread(WINDname,'WG_430'); % [m/s]
w_spd  =  ncread(WINDname,'WS_401');
t_w    =  ncread(WINDname,'time');

rh     =  ncread(RHname,'RH_910'); % RELATIVE HUMIDITY [%]
zrh    = -ncread(RHname,'depth');

rain   =  ncread(RAINname,'RN_485'); % RAIN RATE [mm/hr]

P      =  ncread(BPname,'BP_915'); % BAROMETRIC PRESSURE [hPa]

Ta     =  ncread(AIRTname,'AT_21'); % AIR TEMPERATURE [C]
zTa    = -ncread(AIRTname,'depth');

sst    =  ncread(SSTname,'T_25'); % [C]

Rs     =  ncread(RADname,'RD_495'); % DOWNGOING SHORTWAVE RADIATION [W/m^2]
t_R    =  ncread(RADname,'time');

Rl     =  ncread(LWname,'Ql_136'); % LONGWAVE RADIATION [W/m^2]

depth_t = ncread(Tname,'depth');
tprof   = ncread(Tname,'T_20');

%% Reduce single dimension, convert class to double

w_u    = double(squeeze(w_u));
w_v    = double(squeeze(w_v));
w_spd  = double(squeeze(w_spd));
w_dir  = double(squeeze(w_dir));
w_gust = double(squeeze(w_gust));
sst    = double(squeeze(sst));
Ta     = double(squeeze(Ta));
Rs     = double(squeeze(Rs));
Rl     = double(squeeze(Rl));
rain   = double(squeeze(rain));
rh     = double(squeeze(rh));
P      = double(squeeze(P));

lon    = double(lon);
lat    = double(lat);
zrh    = double(zrh);
zw     = double(zw);
zTa    = double(zTa);

tprof   = double(squeeze(squeeze(tprof)));
depth_t = double(depth_t);

t_w = double(t_w);
t_R = double(t_R);

% ------ absolute time ----------------------------------------------------
t_w  = datenum(2011,2,9,0,0,0) + t_w/60/24;
t_R  = datenum(2011,2,9,0,0,0) + t_R/60/24;
time = t_w; clear t_w

%% Pre-processing

w_u   = interp1(time(w_u<1e5),w_u(w_u<1e5),time);
w_v   = interp1(time(w_v<1e5),w_v(w_v<1e5),time);
Ta    = interp1(time(Ta<1e5), Ta(Ta<1e5),  time);
rh    = interp1(time(rh<1e5), rh(rh<1e5),  time);
sst   = interp1(time(sst<1e5),sst(sst<1e5),time);
P     = interp1(time(P<1e5),  P(P<1e5),    time);
Rs    = interp1(t_R(Rs<1e5),  Rs(Rs<1e5),  t_R);
Rl    = interp1(t_R(Rl<1e5),  Rl(Rl<1e5),  t_R);

w_gust = interp1(time(w_gust<1e5),        w_gust(w_gust<1e5),      time);
rain   = interp1(time(rain>=0 & rain<1e5),rain(rain>=0 & rain<1e5),time);
w_spd  = sqrt(w_u.^2+w_v.^2);

%% Save

clear data_dir *name
save('~/GDrive/UW/Research/Data/OCSP/Mooring/ocsp_2011high.mat');
