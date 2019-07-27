%% read_ocsp_TS_min

% Read original 10 min time series of T profile and 1 hour S profile from OCSP mooring.
% 10 min T profile is averaged in an hour to get hourly profile.

% Zhihua Zheng, UW-APL, July 17 2019

%% General setting

clear
data_dir = '~/GDrive/UW/Research/Data/OCSP/Mooring/L2009_2019/';

Tname = fullfile(data_dir,'t50n145w_10m.cdf');
Sname = fullfile(data_dir,'s50n145w_hr.cdf');

%% Read variables

time_t  = ncread(Tname,'time');
depth_t = ncread(Tname,'depth');
tprof   = ncread(Tname,'T_20');
qtprof  = ncread(Tname,'QT_5020');
lon     = ncread(Tname,'lon');
lat     = ncread(Tname,'lat');

time_s  = ncread(Sname,'time');
depth_s = ncread(Sname,'depth');
sprof   = ncread(Sname,'S_41');
qsprof  = ncread(Sname,'QS_5041');

%% Class

time_t  = double(time_t);
depth_t = double(depth_t);
time_s  = double(time_s);
depth_s = double(depth_s);
lon     = double(lon);
lat     = double(lat);

tprof  = double(squeeze(tprof));
sprof  = double(squeeze(sprof));
qtprof = double(squeeze(qtprof));
qsprof = double(squeeze(qsprof));

%% Pre-processing

time_t = datenum(2007,6,7,23,0,0) + time_t/60/24;
datm_t = datetime(time_t,'ConvertFrom','datenum');

time_s = datenum(2009,6,13,0,0,0) + time_s/24;
datm_s = datetime(time_s,'ConvertFrom','datenum');

tprof(qtprof==0 | qtprof==4 | qtprof==5) = NaN;
sprof(qsprof==0 | qsprof==4 | qsprof==5) = NaN;

% Seabird SBE37 MicroCAT at 1m, 10m has temperature resolution of 0.001C
% Seabird SBE39 at 5m has temperature resolution of 0.0001C
% tprof here has already been filtered by 13 points Hanning window

[TM,Z] = meshgrid(time_t,-depth_t);figure;
contourf(TM,Z,tprof,'linestyle','none');datetick('x','yyyy');ylim([-100 0])

tprof  = tprof';
datmE  = datm_s - minutes(30);

Tprof  = timetable(datm_t,tprof);
THprof = retime(Tprof,datmE,'mean'); % hourly average
THprof.datm_t = datm_s;
% test   = retime(Tprof,datm_s,'mean');

sproff  = vert_fill(sprof,depth_s,time_s);

SAprof  = gsw_SA_from_SP(sprof,depth_s,lon,lat)';
SAproff = gsw_SA_from_SP(sproff,depth_s,lon,lat)';
SAproff = interp1(depth_s,SAproff',depth_t,'linear','extrap')';

PTprof  = gsw_pt0_from_t(SAproff,THprof.tprof,depth_t);
datm    = datm_s;

[TM,Z] = meshgrid(time_s,-depth_t);figure;
contourf(TM,Z,PTprof','linestyle','none');datetick('x','yyyy');ylim([-100 0])

%% Buoyancy

g      = 9.81;
CTprof = gsw_CT_from_pt(SAproff,PTprof);
PDprof = gsw_sigma0(SAproff,CTprof);
contourf(TM,Z,PDprof','linestyle','none');datetick('x','yyyy')

ALPHAprof = gsw_alpha(SAproff,CTprof,depth_t);
BETAprof  = gsw_beta(SAproff,CTprof,depth_t);

Bprof   = g*(ALPHAprof.*CTprof - BETAprof.*SAproff);
NSQprof = center_diff(Bprof,-depth_t',2,'mid'); % buoyancy frequency squared [1/s^2]

%% Timestable

PROF = timetable(datm,PTprof,SAprof,PDprof,Bprof,NSQprof);

%% Save

save('~/GDrive/UW/Research/Data/OCSP/Mooring/ocsp_prof_hrBox.mat',...
     'depth_t','depth_s','PROF');
