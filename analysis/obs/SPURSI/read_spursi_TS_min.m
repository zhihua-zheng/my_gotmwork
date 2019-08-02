
%% read_spursi_TS_min
%
% Script to read subsurface observational data from SPURS-I mooring
% MAT file is saved as 'spursi_prof_hrBox.mat'.

% Zhihua Zheng, UW-APL, July 17 2019

%% General setting

clear
data_dir = '~/GDrive/UW/Research/Data/SPURSI/Mooring/';

TSname = fullfile(data_dir,'SPURS_2012_D_TS.nc'); % 5 min TS profiles 
Mname  = fullfile(data_dir,'SPURS_2012_D_M.nc');  % 1 min Meteorology

%% read variables 

% -------- days since 1950-01-01T00:00:00Z --------------------------------
lat     = ncread(TSname,'LATITUDE');
lon     = ncread(TSname,'LONGITUDE');
time    = ncread(TSname,'TIME');

depth_t = ncread(TSname,'DEPTH');
tprof   = ncread(TSname,'TEMP');
depth_s = ncread(TSname,'DEPTH');
sprof   = ncread(TSname,'PSAL');
instM   = ncread(TSname,'INST_MODEL')';

sss1 = ncread(Mname,'PSAL');
sst1 = ncread(Mname,'TEMP');
tss  = ncread(Mname,'TIME');
depth_ss = ncread(Mname,'DEPTH');

%% Adjust timestamp

% ------ absolute time ----------------------------------------------------
time = datenum(1950,1,1,0,0,0) + time;
tss  = datenum(1950,1,1,0,0,0) + tss;
datm = datetime(time,'ConvertFrom','datenum');
dss  = datetime(tss,'ConvertFrom','datenum');

%% Pre-processing

% SBE37 - accuracy 0.002C 0.0003S/m - resolution 0.0001C 0.00001S/m [sst,sss]
% SBE39 - accuracy 0.002C           - resolution 0.0001C
% SBE37 - accuracy 0.002C 0.0003S/m - resolution 0.0001C 0.00001S/m

% tprof(1,:) = round(tprof(1,:),4);
% tprof(2,:) = round(tprof(2,:),4);
% tprof(3,:) = round(tprof(3,:),4);

depth_s = [depth_ss; depth_s];
depth_t = [depth_ss; depth_t];
sss = interp1(tss,sss1,time)';
sst = interp1(tss,sst1,time)';
sprof = [sss; sprof]; % includes sss
tprof = [sst; tprof]; % includes sst

sproff = vert_fill(sprof,depth_s,time); % fill some gaps to preserve temp data points

SAproff = gsw_SA_from_SP(sproff,depth_s,lon,lat)';
SAprof  = gsw_SA_from_SP(sprof,depth_s,lon,lat)';

% filled salinity is used to help derive other variables, as only subtle
% difference exists between SAprof and SAproff.
PTprof  = gsw_pt0_from_t(SAproff,tprof',depth_t);

%% Buoyancy

g      = 9.81;
rho0   = 1025;
CTprof = gsw_CT_from_pt(SAproff,PTprof);
PDprof = gsw_sigma0(SAproff,CTprof);

ALPHAprof = gsw_alpha(SAproff,CTprof,depth_t);
BETAprof  = gsw_beta(SAproff,CTprof,depth_s);
Bprof     = g*(ALPHAprof.*CTprof - BETAprof.*SAproff);

% mPDprof    = mean(PDprof,1);
% PDprof_bar = interp1(-depth_t(~isnan(mPDprof)),mPDprof(~isnan(mPDprof)),-depth_t');
% Bprof   = -g/rho0*(1000+PDprof);

NSQprof = center_diff(Bprof,-depth_t',2,'mid'); % buoyancy frequency squared [1/s^2]

%% Timetable

Prof = timetable(datm,PTprof,SAprof,PDprof,Bprof,NSQprof);
PROF = retime(Prof,'hourly','mean');

datm_hr = (datm(1):hours(1):datm(end))' + minutes(30);
PROF.datm = datm_hr;

%% Save

save('~/GDrive/UW/Research/Data/SPURSI/Mooring/spursi_prof_hrBox.mat',...
     'depth_t','depth_s','PROF','Prof');
