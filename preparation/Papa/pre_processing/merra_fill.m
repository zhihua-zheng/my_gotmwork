function [P_r,sst_r,rain_r] = merra_fill(P,sst,rain,lat,lon,time)

% merra_fill
%==========================================================================
%
% USAGE:
%  [P_r,sst_r,rain_r] = merra_fill(P,sst,rain,lat,lon,time)
%
% DESCRIPTION:
%  Function to fill the large gaps in meteorological measurements
%    from PMEL OCSPapa mooring, using reanalysis product
%    MERRA2 (https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/data_access/)
%
%  ocean_rainfall and open_water_skin_temperature are from Ocean Surface
%  Diagnostics dataset 'tavg1_2d_ocn_Nx'. [-145.625,49.5,-144.375,50.5]
%  (https://disc.gsfc.nasa.gov/datasets/M2T1NXOCN_V5.12.4/summary?keywords=%22MERRA-2%22)
%
%  sea_level_pressure is from Single-Level Diagnostics dataset
%  'tavg1_2d_slv_Nx'. [-145.625,49.5,-144.375,50.5]
%  (https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_V5.12.4/summary?keywords=%22MERRA-2%22)
%
% INPUT:
%
%  P    - 1-D vector, time series of observed air pressure [hPa]
%  sst  - 1-D vector, time series of observed sea surface temp. [C]
%  rain - 1-D vector, time series of observed rain rate [mm/hr]
%  lat  - scalar, latitude  of buoy site (-90 ~ 90 )
%  lon  - scalar, longitude of buoy site (  0 ~ 360)
%  time - 1-D vector, timestamp of observation time series
%
% OUTPUT:
%
%  P_r    - 1-D vector, air pressure with gap filled [hPa]
%  sst_r  - 1-D vector, sea surface temp. with gap filled [C]
%  rain_r - 1-D vector, rain rate with gap filled [mm/hr]
%
% AUTHOR:
%  Jun. 21 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%  Mar. 26 2019. converted to function
%==========================================================================

%% General setting

merra_dir = '~/GDrive/UW/Research/Data/OCSP/MERRA';

P_r    = P;
sst_r  = sst;
rain_r = rain;

%% Sea level pressure

bad_P = find(P>10000);

ncvars =  {'lat','lon','SLP','time'};
dinfo = dir(fullfile(merra_dir,'/slp/*.nc'));
nFile = length(dinfo);
filenames = fullfile({dinfo.folder},{dinfo.name});

lat_merra = ncread(filenames{1}, ncvars{1});
lon_merra = ncread(filenames{1}, ncvars{2});

slp      = nan(length(lat_merra),length(lon_merra),24,nFile);
t_merra  = nan(24,nFile);

for K = 1:nFile

  this_file = filenames{K};
  slp(:,:,:,K) = ncread(this_file, ncvars{3});

  % get the reference time
  t_ref = ncreadatt(this_file,'time','units'); % get the attribute 'units' for 'time'
  t_ref = t_ref(15:end); % truncate to the time string
  t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS');

  % date numbers for the timestamp
  t_merra(:,K) = double(ncread(this_file, ncvars{4}))/60/24 + t_ref;

end

[X_merra, Y_merra, T_merra] = meshgrid(lon_merra,lat_merra,t_merra(:));
slp = slp(:,:,:);

% interpolate to get the barometric pressure at buoy site
slp_papa = interp3(X_merra,Y_merra,T_merra,slp,...
    (lon-360)*ones(size(bad_P)),lat*ones(size(bad_P)),time(bad_P))/100;
% unit for sea level pressure in merra is Pa

P_r(bad_P) = slp_papa;

%% open_water_skin_temperature and ocean_rain_fall

bad_sst  = find(sst>10000);
bad_rain = find(rain>10000);

ncvars =  {'lat', 'lon', 'TSKINWTR', 'RAINOCN', 'time'};
dinfo = dir(fullfile(merra_dir,'/rain_skinT/*.nc'));
nFile = length(dinfo);
filenames = fullfile({dinfo.folder},{dinfo.name});

lat_merra = ncread(filenames{1}, ncvars{1});
lon_merra = ncread(filenames{1}, ncvars{2});

ts       = nan(length(lat_merra),length(lon_merra),24,nFile); % skin temp.
rainfall = nan(length(lat_merra),length(lon_merra),24,nFile); % rain rate
t_merra = nan(24,nFile);

for K = 1:nFile

  this_file = filenames{K};
  ts(:,:,:,K)       = ncread(this_file, ncvars{3});
  rainfall(:,:,:,K) = ncread(this_file, ncvars{4});

  % get the reference time
  t_ref = ncreadatt(this_file,'time','units'); % get the attribute 'units' for 'time'
  t_ref = t_ref(15:end); % truncate to the time string
  t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS');

  % timestamp
  t_merra(:,K) = double(ncread(this_file, ncvars{5}))/60/24 + t_ref;

end

[X_merra, Y_merra, T_merra] = meshgrid(lon_merra,lat_merra,t_merra(:));
ts       = ts(:,:,:);
rainfall = rainfall(:,:,:);

% interpolate to get the SST and rain rate at buoy site
ts_papa   = interp3(X_merra,Y_merra,T_merra,ts,...
    (lon-360)*ones(size(bad_sst)),lat*ones(size(bad_sst)),...
    time(bad_sst))-273.15;
% unit for surface skin temperature in MERRA is Kelvin

rain_papa = interp3(X_merra,Y_merra,T_merra,rainfall,...
    (lon-360)*ones(size(bad_rain)),lat*ones(size(bad_rain)),...
    time(bad_rain),'linear',0)*3600; % [mm/hr]
% unit for rain rate in MERRA is kg/(m^2*s) ~ mm/s (no salt in rain)

sst_r(bad_sst)   = ts_papa;
rain_r(bad_rain) = rain_papa;

end

