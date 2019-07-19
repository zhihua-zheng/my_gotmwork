%% read_ocsp_TS_hr

% Read hourly TS profile data from OCSP mooring. The original temporal resolutin for S profile is
% 1 hour but for T profile is 10 min. The hourly T profile downloaded from OCSPcwebsite results 
% from a 13 point Hanning filter applied to the original 10 min time series.

% Zhihua Zheng, UW-APL, July 15 2019

%% General setting

clear
data_dir = '~/GDrive/UW/Research/Data/OCSP/Mooring/L2009_2019/';

Tname = fullfile(data_dir,'t50n145w_hr.cdf');
Sname = fullfile(data_dir,'s50n145w_hr.cdf');

%% Read variables

time    = ncread(Tname,'time');
depth_t = ncread(Tname,'depth');
tprof   = ncread(Tname,'T_20');
qtprof  = ncread(Tname,'QT_5020');
lon     = ncread(Tname,'lon');
lat     = ncread(Tname,'lat');

depth_s = ncread(Sname,'depth');
sprof   = ncread(Sname,'S_41');
qsprof  = ncread(Sname,'QS_5041');

%% Class

time    = double(time);
depth_t = double(depth_t);
depth_s = double(depth_s);
lon     = double(lon);
lat     = double(lat);

tprof  = double(squeeze(tprof));
sprof  = double(squeeze(sprof));
qtprof = double(squeeze(qtprof));
qsprof = double(squeeze(qsprof));

%% Pre-processing

time = datenum(2009,6,13,0,0,0) + time/24;
datm = datetime(time,'ConvertFrom','datenum');

tprof(qtprof==0 | qtprof==4 | qtprof==5) = NaN;
sprof(qsprof==0 | qsprof==4 | qsprof==5) = NaN;

% Seabird SBE37 MicroCAT at 1m, 10m has temperature resolution of 0.001C
% Seabird SBE39 at 5m has temperature resolution of 0.0001C
% tprof here has already been filtered by 13 points Hanning window

[TM,Z] = meshgrid(time,-depth_t);figure;
contourf(TM,Z,tprof,'linestyle','none');datetick('x','yyyy')

sproff = vert_fill(sprof,depth_s,time);
SAprof = gsw_SA_from_SP(sproff,depth_s,lon,lat);
SAprof = interp1(depth_s,SAprof,depth_t,'linear','extrap')';
PTprof = gsw_pt0_from_t(SAprof,tprof',depth_t);

figure;contourf(TM,Z,PTprof','linestyle','none');datetick('x','yyyy')

%% Timestable

PROF = timetable(datm,PTprof,SAprof);

%% Save

save('~/GDrive/UW/Research/Data/OCSP/Mooring/ocsp_prof_hr.mat',...
     'depth_t','depth_s','PROF');
 
