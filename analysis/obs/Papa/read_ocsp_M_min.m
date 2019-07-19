
%% read_ocsp_M_min
%
% This the script used to read the meterological and subsurface observation
% data from OCSPapa mooring, the MAT file is saved as 'ocspM_2011min.mat'.

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
Sname    = fullfile(data_dir,'s50n145w_hr.cdf');

%% read variables 

% -------- minutes since 2011-02-09 00:00:00 --------------------------------
lat    =  ncread(WINDname,'lat');
lon    =  ncread(WINDname,'lon');
wu     =  ncread(WINDname,'WU_422');
wv     =  ncread(WINDname,'WV_423');
zw     = -ncread(WINDname,'depth');
wdir   =  ncread(WINDname,'WD_410'); % [degree]
wgust  =  ncread(WINDname,'WG_430'); % [m/s]
t_w    =  ncread(WINDname,'time');

rh     =  ncread(RHname,'RH_910'); % RELATIVE HUMIDITY [%]
zrh    = -ncread(RHname,'depth');

rain   =  ncread(RAINname,'RN_485'); % RAIN RATE [mm/hr]

P      =  ncread(BPname,'BP_915'); % BAROMETRIC PRESSURE [hPa]

Ta     =  ncread(AIRTname,'AT_21'); % AIR TEMPERATURE [C]
zTa    = -ncread(AIRTname,'depth');

sst    =  ncread(SSTname,'T_25'); % [C]

Rs1    =  ncread(RADname,'RD_495'); % DOWNGOING SHORTWAVE RADIATION [W/m^2]
t_R    =  ncread(RADname,'time');

Rl1    =  ncread(LWname,'Ql_136'); % LONGWAVE RADIATION [W/m^2]

depth_t = ncread(Tname,'depth');
tprof   = ncread(Tname,'T_20');
tq      = ncread(Tname,'QT_5020');
depth_s = ncread(Sname,'depth');
sprof   = ncread(Sname,'S_41');

%% Reduce single dimension, convert class to double

wu    = double(squeeze(wu));
wv    = double(squeeze(wv));
wdir  = double(squeeze(wdir));
wgust = double(squeeze(wgust));
sst   = double(squeeze(sst));
Ta    = double(squeeze(Ta));
Rs1   = double(squeeze(Rs1));
Rl1   = double(squeeze(Rl1));
rain  = double(squeeze(rain));
rh    = double(squeeze(rh));
P     = double(squeeze(P));

zw  = double(zw);
zrh = double(zrh);
zTa = double(zTa);
lat = double(lat);
lon = double(lon);

t_w = double(t_w);
t_R = double(t_R);

tprof   = double(squeeze(tprof));
tq      = double(squeeze(tq));
sprof   = double(squeeze(sprof));
depth_t = double(depth_t);
depth_s = double(depth_s);

% ------ absolute time ----------------------------------------------------
t_w  = datenum(2011,2,9,0,0,0) + t_w/60/24;
t_R  = datenum(2011,2,9,0,0,0) + t_R/60/24;
time = t_w; clear t_w

%% Pre-processing

wu    = interp1(time(wu<1e5), wu(wu<1e5),  time);
wv    = interp1(time(wv<1e5), wv(wv<1e5),  time);
Ta    = interp1(time(Ta<1e5), Ta(Ta<1e5),  time);
rh    = interp1(time(rh<1e5), rh(rh<1e5),  time);
sst   = interp1(time(sst<1e5),sst(sst<1e5),time);
P     = interp1(time(P<1e5),  P(P<1e5),    time);
Rs1   = interp1(t_R(Rs1<1e5), Rs1(Rs1<1e5),  t_R);
Rl1   = interp1(t_R(Rl1<1e5), Rl1(Rl1<1e5),  t_R);

Rs  = interp1(t_R,Rs1,time);
Rl  = interp1(t_R,Rl1,time);

wgust = interp1(time(wgust<1e5),         wgust(wgust<1e5),       time);
rain  = interp1(time(rain>=0 & rain<1e5),rain(rain>=0 & rain<1e5),time);

datm = datetime(time,'ConvertFrom','datenum');
AQ   = timetable(datm,wu,wv,Ta,rh,sst,P,Rs,Rl,rain);
clear datm wu wv Ta rh sst P Rs Rl rain

%% Save

clear data_dir *name
save('~/GDrive/UW/Research/Data/OCSP/Mooring/ocspM_2011min.mat');
