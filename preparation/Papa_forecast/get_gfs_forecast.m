function gfs = get_gfs_forecast(lon,lat,t_0)

% get_gfs_forecast
%==========================================================================
%
% USAGE:
%  gfs = get_gfs_forecast(lon,lat,t_0)
%
% DESCRIPTION:
%  Function to retrieve meteorological forecast data from Global Forecast 
%  System (GFS).
%
% INPUT:
%
%  lon - longitude for interested station ([0~360),decimal degree)
%  lat - latitude for interested station ([-90~90],decimal degree)
%  t_0 - starting time for forecast time series (YYYYMMDD)
%
% OUTPUT:
%
%  Struct containing meteorological forcing info for given location
%
% AUTHOR:
%  September 15 2018. Eric D'Asaro                       
%  September 17 2018. Zhihua Zheng                     [ zhihua@uw.edu ]   

gfs = struct();

%% set the source and get time

% OWS P nominal location 50.1°N, 144.9°W
% lon0 =  360-(144+25/60);
% lat0 =  50+22/60;
lon0 = lon;
lat0 = lat;

% GrADS data server for GFS forecast, using 0.25 degree resolution 
url_fmt = 'http://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs%s/gfs_0p25_1hr_00z'; % there are also 06z, 12z, 18z...
% chosse 00z since the T-S profile data at OCSPapa is published as 00z

% get the forecast since t_0 (could also grab a few days back for comparison):
gfs.url =  sprintf(url_fmt,t_0);
disp(['Get the GFS forecast since ', t_0, ' ...']);

%nc = ncinfo(gfs.url); % don't really need it once we know variable names

t = ncread(gfs.url,'time');
% Mtime = t-736953+datenum(2018,9,15); % their date reference is somewhat strange, this fixes it
t_ref = ncreadatt(gfs.url,'time','units'); % get the attribute 'unit' for 'time'
disp(['The unit for time in GFS netCDF file is: ', t_ref])
gfs.time = 365 + t;
gfs.yd = gfs.time - datenum('01-Jan-2018') + 1; 
gfs.date = datestr(gfs.time,'yyyy-mm-dd HH:MM:SS');

%% find the closest grid point indices:
lon = ncread(gfs.url,'lon');
lat = ncread(gfs.url,'lat');
[~,ilon] = min(abs(lon-lon0));
[~,ilat] = min(abs(lat-lat0));
% disp([ilon lon(ilon)-360 ilat lat(ilat) ])

% ... or just hard-wire for OWS-P location (50,-145)
% ilon = 861;
% ilat = 561;
% return

%% wind speed...
disp('get wind speed ...')
gfs.u = squeeze(ncread(gfs.url,'ugrd10m',[ilon ilat 1],[1 1 inf])); %** 10 m above ground u-component of wind [m/s] 
gfs.v = squeeze(ncread(gfs.url,'vgrd10m',[ilon ilat 1],[1 1 inf])); %** 10 m above ground v-component of wind [m/s] 
W = complex(gfs.u,gfs.v);

%% wind stress...  Switch sign to oceanographic convention
disp('get wind stress ...');

% sign switched to represent surface stress
gfs.tau_x = -squeeze(ncread(gfs.url,'uflxsfc',[ilon ilat 1],[1 1 inf])); % ** surface momentum flux, u-component [n/m^2] 
gfs.tau_y = -squeeze(ncread(gfs.url,'vflxsfc',[ilon ilat 1],[1 1 inf])); % ** surface momentum flux, v-component [n/m^2] 
% tau = complex(gfs.tau_x,gfs.tau_y);

%% Fluxes
disp('get heat fluxes ...')
latent= -squeeze(ncread(gfs.url,'lhtflsfc',[ilon ilat 1],[1 1 inf])); % sign switched
sens=   -squeeze(ncread(gfs.url,'shtflsfc',[ilon ilat 1],[1 1 inf])); % sign switched

disp('get radiations ...')
swd=squeeze(ncread(gfs.url,'dswrfsfc',[ilon ilat 1],[1 1 inf]));
swu=squeeze(ncread(gfs.url,'uswrfsfc',[ilon ilat 1],[1 1 inf]));
lwd=squeeze(ncread(gfs.url,'dlwrfsfc',[ilon ilat 1],[1 1 inf]));
lwu=squeeze(ncread(gfs.url,'ulwrfsfc',[ilon ilat 1],[1 1 inf]));
rain=squeeze(ncread(gfs.url,'pratesfc',[ilon ilat 1],[1 1 inf])); % kg/m^2/s
rain=rain/1000; % m/s
gfs.rain = rain*1000/3600; % mm/hr

%% Flux processing and save data

% double check: downward - upward; only has downward obs. in Papa mooring
gfs.sw=swd-swu;
lw=lwd-lwu;

% double check: short wave radiation is not included in total heat flux (GOTM)
Q=lw+gfs.sw+latent+sens;  % CHECKED SIGNS WITH NOAA MOORING
gfs.hf = lw + latent + sens;

t_stamp=datestr(now,'yymmddHHMMSS');

% disp(['save ' t_stamp '.mat'])
% save(['GFS/' t_stamp],'yd','time','latent','sens','sw','lw','Q','rain','u','v','tau_u','tau_v','url');

%% visualization

disp('plot ...')

figure('position',[ 520   221   656   677 ]);

subplot(3,1,1);
plot(gfs.yd,gfs.sw,gfs.yd,lw,gfs.yd,latent,gfs.yd,sens);
hold on
plot(gfs.yd,Q,'k-','LineWidth',2);
legend('sw','lw','latent','sens','Total');
title(gfs.url,'interpreter','none');

subplot(3,1,2);
plot(gfs.time,gfs.u,gfs.time,gfs.v);
hold on
plot(gfs.time,abs(W),'k-','LineWidth',2);
grid on
hold on
datetick x keeplimits
ylabel('Wind speed (m/s)');
legend('u','v','Sp10');
title('Wind Speed');

subplot(3,1,3)
quiver(gfs.time,zeros(size(gfs.time)),-gfs.tau_x,-gfs.tau_y,2,'Marker','.','Color','k','ShowArrowHead','off');
axis equal
datetick x keeplimits
title('Stress');

%% save the fig
print([ 'figs/' t_stamp '.png' ],'-dpng');
disp('done!')

end
