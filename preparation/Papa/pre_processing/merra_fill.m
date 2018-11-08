%% merra_fill

% subrountine used to fill the large gaps in meteorological measurements
% from PMEL ocean climate station Papa mooring, using reanalysis product
% from MERRA2 (https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/data_access/)

% ocean rainfall and open water skin temperature are from Ocean Surface 
% Diagnostics dataset 'tavg1_2d_ocn_Nx'.

% sea level pressure is from Single-Level Diagnostics dataset
% 'inst1_2d_asm_Nx'.

% by Zhihua Zheng (UW/APL), updated on Jun. 21 2018

%% slp

% large gap in barometric pressure after pre_day (Mar. 2010 - Nov.2010)
bad_P = find(P>10000);
% the index of outliers in P time series


ncvars =  {'lat', 'lon', 'SLP', 'time'};
merra_dir = ['~/Documents/Study/GraduteResearch/Data/OCS_P/',...
    'Met_Forcing/original_data_2007/MERRA/slp'];
dinfo = dir(fullfile(merra_dir,'*.nc'));
num_files = length(dinfo);
filenames = fullfile(merra_dir,{dinfo.name});

lat_merra = ncread(filenames{1}, ncvars{1});
lon_merra = ncread(filenames{1}, ncvars{2});

slp = zeros(length(lat_merra),length(lon_merra),24*num_files)*NaN;
time_slp = zeros(24*num_files,1)*NaN;

for K = 1:24:(24*num_files-23)
    
  this_file = filenames{ceil(K/24)};
  slp(:,:,K:K+23) = ncread(this_file, ncvars{3});
  
  % get the reference time
  t_ref = ncreadatt(this_file,'time','units'); % get the attribute 'units' for 'time'
  t_ref = t_ref(15:end); % truncate to the time string
  t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS');
  
  % date numbers for the timestamp
  time_slp(K:K+23) = double(ncread(this_file, ncvars{4}))/(60*24) + t_ref;
  
end

% date_slp = string(datestr(time_slp,'yyyy/mm/dd HH:MM:SS'));
[X_merra, Y_merra, Z_merra] = meshgrid(lon_merra,lat_merra,time_slp);

% interpolate to get the sea level pressure time series at buoy site
slp_papa = interp3(X_merra,Y_merra,Z_merra,slp,...
    (lon-360)*ones(size(bad_P)),lat*ones(size(bad_P)),time(bad_P))/100;
% the unit for sea level pressure in merra is Pa

P_r = P;
P_r(bad_P) = slp_papa;


%% sst

bad_sst = find(sst>10000);


ncvars =  {'lat', 'lon', 'TSKINWTR', 'time'};
merra_dir = ['~/Documents/Study/GraduteResearch/Data/OCS_P', ...
    '/Met_Forcing/original_data_2007/MERRA/rain_skin_t'];
dinfo = dir(fullfile(merra_dir,'*.nc'));
num_files = length(dinfo);
filenames = fullfile(merra_dir,{dinfo.name});

lat_merra = ncread(filenames{1}, ncvars{1});
lon_merra = ncread(filenames{1}, ncvars{2});

ts = zeros(length(lat_merra),length(lon_merra),24*num_files)*NaN; % skin temperature
time_sst = zeros(24*num_files,1)*NaN;

for K = 1:24:(24*num_files-23)
    
  this_file = filenames{ceil(K/24)};
  ts(:,:,K:K+23) = ncread(this_file, ncvars{3});
  
  % get the reference time
  t_ref = ncreadatt(this_file,'time','units'); % get the attribute 'units' for 'time'
  t_ref = t_ref(15:end); % truncate to the time string
  t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS');
  
  % date numbers for the timestamp
  time_sst(K:K+23) = double(ncread(this_file, ncvars{4}))/(60*24) + t_ref;
  
end

[X_merra, Y_merra, Z_merra] = meshgrid(lon_merra,lat_merra,time_sst);

% interpolate to get the sea level pressure time series at buoy site
ts_papa = interp3(X_merra,Y_merra,Z_merra,ts,...
    (lon-360)*ones(size(bad_sst)),lat*ones(size(bad_sst)),time(bad_sst))-273.15;
% the unit for surface skin temperature in merra is Kelvin

sst_r = sst;
sst_r(bad_sst) = ts_papa;

%% rain

bad_rain = find(rain>10000);

bad_rain = [bad_rain; (length(time_rain)+1:length(time))'];
% incoporate the period after the last observation record of rain rate

ncvars =  {'lat', 'lon', 'RAINOCN', 'time'};
merra_dir = ['~/Documents/Study/GraduteResearch/Data/OCS_P', ...
    '/Met_Forcing/original_data_2007/MERRA/rain_skin_t'];
dinfo = dir(fullfile(merra_dir,'*.nc'));
num_files = length(dinfo);
filenames = fullfile(merra_dir,{dinfo.name});

lat_merra = ncread(filenames{1}, ncvars{1});
lon_merra = ncread(filenames{1}, ncvars{2});

rainfall = zeros(length(lat_merra),length(lon_merra),24*num_files)*NaN; % rain rate
time_rainfall = zeros(24*num_files,1)*NaN;

for K = 1:24:(24*num_files-23)
    
  this_file = filenames{ceil(K/24)};
  rainfall(:,:,K:K+23) = ncread(this_file, ncvars{3});
  
  % get the reference time
  t_ref = ncreadatt(this_file,'time','units'); % get the attribute 'units' for 'time'
  t_ref = t_ref(15:end); % truncate to the time string
  t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS');
  
  % date numbers for the timestamp
  time_rainfall(K:K+23) = double(ncread(this_file, ncvars{4}))/(60*24) + t_ref;
  
end

[X_merra, Y_merra, Z_merra] = meshgrid(lon_merra,lat_merra,time_rainfall);

% interpolate to get the sea level pressure time series at buoy site
rainfall_papa = interp3(X_merra,Y_merra,Z_merra,rainfall,...
    (lon-360)*ones(size(bad_rain)),lat*ones(size(bad_rain)),time(bad_rain))*3600;
% the unit for rain rate in merra is kg/(m^2*s) ~ mm/s (no salt in rain)

rain_r = [rain; ones(length(time_rain)-length(time))];
rain_r(bad_rain) = rainfall_papa;

%% clean

clear t_ref

