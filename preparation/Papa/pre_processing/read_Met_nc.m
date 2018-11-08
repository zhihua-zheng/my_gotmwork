%% read_Met_nc.m
%
% This the script used to read the information from meterological focing
% files and the generated workspace is named as 'met_forcing_p2007.mat'.

% by Zhihua Zheng (UW/APL), updated on Jun. 15 2018

%% read variables from netCDF files

clear
data_dir = '~/Documents/Study/GraduteResearch/Data/OCS_P/Met_Forcing/original_data_2007/';

lat = ncread([data_dir,'w50n145w_hr.cdf'],'lat');
lon = ncread([data_dir,'w50n145w_hr.cdf'],'lon');
w_u = ncread([data_dir,'w50n145w_hr.cdf'],'WU_422');
w_v = ncread([data_dir,'w50n145w_hr.cdf'],'WV_423');
z_wind = - ncread([data_dir,'w50n145w_hr.cdf'],'depth');
w_dir = ncread([data_dir,'w50n145w_hr.cdf'],'WD_410');
w_gust = ncread([data_dir,'w50n145w_hr.cdf'],'WG_430');
w_spd = ncread([data_dir,'w50n145w_hr.cdf'],'WS_401');
sst = ncread([data_dir,'sst50n145w_hr.cdf'],'T_25');
sss = ncread([data_dir,'sss50n145w_hr.cdf'],'S_41');
ssd = ncread([data_dir,'ssd50n145w_hr.cdf'],'STH_71');
rh = ncread([data_dir,'rh50n145w_hr.cdf'],'RH_910');
z_rh = - ncread([data_dir,'rh50n145w_hr.cdf'],'depth');
rain = ncread([data_dir,'rain50n145w_hr.cdf'],'RN_485');
Rs = ncread([data_dir,'rad50n145w_hr.cdf'],'RD_495');
Rl = ncread([data_dir,'lw50n145w_hr.cdf'],'Ql_136');
P = ncread([data_dir,'bp50n145w_hr.cdf'],'BP_915');
t_air = ncread([data_dir,'airt50n145w_hr.cdf'],'AT_21');
z_ta = - ncread([data_dir,'airt50n145w_hr.cdf'],'depth');
time = ncread([data_dir,'rad50n145w_hr.cdf'],'time'); % hours since 2007-06-08 05:00:00
time_rain = ncread([data_dir,'rain50n145w_hr.cdf'],'time'); % hours since 2007-06-08 04:00:00
dn2007 = datenum(2007,6,8,5,0,0);
dn2007_rain = datenum(2007,6,8,4,0,0);


%% read profiles from netCDF files

depth_s = ncread([data_dir,'s50n145w_hr.cdf'],'depth');
sprof = ncread([data_dir,'s50n145w_hr.cdf'],'S_41');
depth_t = ncread([data_dir,'t50n145w_hr.cdf'],'depth');
tprof = ncread([data_dir,'t50n145w_hr.cdf'],'T_20');
depth_rho = ncread([data_dir,'d50n145w_hr.cdf'],'depth');
rhoprof = ncread([data_dir,'d50n145w_hr.cdf'],'STH_71');

% uprof_adcp = ncread('adcp50n145w_mon.cdf'],'u_1205');
% vprof_adcp = ncread('adcp50n145w_mon.cdf'],'v_1206');
% depth_adcp = ncread('adcp50n145w_mon.cdf'],'depth');
% time_adcp = ncread('adcp50n145w_mon.cdf'],'time'); % days since 2010 Jun 16 12:00:00
% 
% uprof_sen = ncread('adcp_sentinel50n145w_mon.cdf'],'u_1205');
% vprof_sen = ncread('adcp_sentinel50n145w_mon.cdf'],'v_1206');
% depth_sen = ncread('adcp_sentinel50n145w_mon.cdf'],'depth');
% time_sen = ncread('adcp_sentinel50n145w_mon.cdf'],'time'); % days since 2010 Jun 16 12:00:00

cur_u = ncread([data_dir,'profiles/cur50n145w_mon.cdf'],'U_320');
cur_v = ncread([data_dir,'profiles/cur50n145w_mon.cdf'],'V_321');
cur_spd = ncread([data_dir,'profiles/cur50n145w_mon.cdf'],'CS_300');
cur_dir = ncread([data_dir,'profiles/cur50n145w_mon.cdf'],'CD_310');
depth_cur = ncread([data_dir,'profiles/cur50n145w_mon.cdf'],'depth');

time_cur = ncread([data_dir,'profiles/cur50n145w_mon.cdf'],'time'); % days since 2007-06-16 12:00:00
dn2007_cur = datenum(2007,6,16,12,0,0);


%% Reduce single dimension, convert class to double, unify beginning time

w_u = double(squeeze(squeeze(w_u)));
w_v = double(squeeze(squeeze(w_v)));
w_spd = double(squeeze(squeeze(w_spd)));
w_dir = double(squeeze(squeeze(w_dir)));
w_gust = double(squeeze(squeeze(w_gust)));
sss = double(squeeze(squeeze(sss)));
sst = double(squeeze(squeeze(sst)));
ssd = double(squeeze(squeeze(ssd)));
t_air = double(squeeze(squeeze(t_air)));
Rs = double(squeeze(squeeze(Rs)));
Rl = double(squeeze(squeeze(Rl)));
rain = double(squeeze(squeeze(rain)));
rh = double(squeeze(squeeze(rh)));
P = double(squeeze(squeeze(P)));
lon = double(lon);
lat = double(lat);
z_rh = double(z_rh);
z_wind = double(z_wind);
z_ta = double(z_ta);

sprof = double(squeeze(squeeze(sprof)));
tprof = double(squeeze(squeeze(tprof)));
rhoprof = double(squeeze(squeeze(rhoprof)));
cur_u = double(squeeze(squeeze(cur_u)));
cur_v = double(squeeze(squeeze(cur_v)));
cur_spd = double(squeeze(squeeze(cur_spd)));
cur_dir = double(squeeze(squeeze(cur_dir)));

time = double(time);
time_rain = double(time_rain);
time_cur = double(time_cur);
depth_s = double(depth_s);
depth_t = double(depth_t);
depth_rho = double(depth_rho); % same as 'depth_s'
depth_cur = double(depth_cur);

% unify the beginning time of the all the time series, as 'time(1)'

% time series one hour earlier than 'dn2007'
P = P(2:end);
rh = rh(2:end);
Rl = Rl(2:end);
t_air = t_air(2:end);
w_dir = w_dir(2:end);
w_gust = w_gust(2:end);
w_spd = w_spd(2:end);
w_u = w_u(2:end);
w_v = w_v(2:end);
rain = rain(2:end); 

% time series 5 hours earlier than 'dn2007'
sss = sss(6:end);
sst = sst(6:end);
ssd = ssd(6:end);
sprof = sprof(:,6:end);
tprof = tprof(:,6:end);
rhoprof = rhoprof(:,6:end);

%% Timestamp

time = time/24+dn2007;
time_cur = time_cur+dn2007_cur;
time_rain = time_rain/24+dn2007_rain;
time_rain = time_rain(2:end); % 'dn2007_rain' is one hour earlier than 'dn2007'

date = datestr(time, 'yyyy/mm/dd HH:MM:SS');
date = string(date);

date_rain = datestr(time_rain, 'yyyy/mm/dd HH:MM:SS');
date_rain = string(date_rain);

date_cur = datestr(time_cur, 'yyyy/mm/dd HH:MM:SS');
date_pcur = string(date_cur);

%% Save

clear data_dir

save('met_forcing_p2007.mat');
