

%% read_spursi_flux

% Read computed air-sea fluxes for SPURS-I, rainfall data is missing

% Zhihua Zheng, UW-APL, July 17 2019

%% General setting

clear

flux_dir = '~/GDrive/UW/Research/Data/SPURSI/Mooring/flux_UOP/';
M_dir    = '~/GDrive/UW/Research/Data/SPURSI/Mooring/';

SFname   = fullfile(flux_dir,'SPURS_2012_D_F_1hr.nc');
Mname    = fullfile(M_dir,   'SPURS_2012_D_M_1hr.nc');

%% Read variables

time  = ncread(SFname,'TIME');

tau   = ncread(SFname,'TAUMAG'); % [N/m^2]
taud  = ncread(SFname,'TAUDIR'); % [degree] 'to' direction, clockwise from N
nsw   = ncread(SFname,'QS');     % [W/m^2] into the ocean
nlw   = ncread(SFname,'QL');     % [W/m^2] into the ocean
hlb   = ncread(SFname,'QH');     % [W/m^2] into the ocean
hsb   = ncread(SFname,'QB');     % [W/m^2] into the ocean
sst   = ncread(SFname,'TSKIN');  % [C] skin temperature

sss   = ncread(Mname,'PSAL');    % [PSU]
rain  = ncread(Mname,'RAIN');    % [mm/hr]
lon   = ncread(Mname,'LONGITUDE');
lat   = ncread(Mname,'LATITUDE');

ssSA = gsw_SA_from_SP(sss,0,lon,lat); % [g/kg]
Le   = gsw_latentheat_evap_t(ssSA,sst); % heat of evaporation for seawater [J/kg]

rho_fw = 1000;
evap   = - hlb ./ Le / rho_fw; % [m/s], mostly positive
evap   = evap*1000*3600;       % [mm/hr]

tau_x  = tau .* sind(taud);
tau_y  = tau .* cosd(taud);

%% Timestamp

time = datenum(1950,1,1,0,0,0) + time;
datm = datetime(time,'ConvertFrom','datenum');

tlower  = dateshift(datm(1),'start','hour');
tupper  = dateshift(datm(end),'start','hour');
datm_hr = (tlower:hours(1):tupper)' + minutes(30);

SF   = timetable(datm,tau_x,tau_y,tau,nsw,nlw,hlb,hsb,evap,rain,sst,sss);
SF.datm = datm_hr;

%% Save

save('~/GDrive/UW/Research/Data/SPURSI/Mooring/spursi_flux_hrUOP.mat','SF');
