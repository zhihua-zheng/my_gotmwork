%% spursi_forcing
%
% make forcing files for SPURS-I simulation in GOTM
%
% Zhihua Zheng, UW-APL, Mar. 7 2019

clear

%% Load Met. forcing, T-S profiles and Wave data

mooring_dir = '~/GDrive/UW/Research/Data/SPURSI/Mooring/';

MFname  = fullfile(mooring_dir,'SPURS_2012_D_M_1hr.nc'); % Meterological Forcing
TSname  = fullfile(mooring_dir,'SPURS_2012_D_TS.nc'); % T-S profiles
WDS_dir = '~/GDrive/UW/Research/Data/SPURSI/WDS/'; % Wave Directional Spectrum

lat_MF  = ncread(MFname,'LATITUDE');
lon_MF  = ncread(MFname,'LONGITUDE');
time    = ncread(MFname,'TIME'); % days since t_ref

% get the reference time
t_ref = ncreadatt(MFname,'TIME','units'); % attribute 'units' for 'TIME'
t_ref = t_ref(12:end); % truncate to get the time string
t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS');
time_MF = t_ref + time; 
clear time t_ref

% z_TS_s = ncread(MFname,'DEPTH'); % depth of sea water temp. and sal. sensor
% z_Rl   = ncread(MFname,'HEIGHT_LW'); % height of longwave radiation sensor
% z_Rs   = ncread(MFname,'HEIGHT_SW'); % height of shortwave radiation sensor
% z_Rain = ncread(MFname,'HEIGHT_RAIN'); % height of rain gauge sensor
% z_Pa   = ncread(MFname,'HEIGHT_ATMS'); % height of barometric pressure sensor
z_RhTa = ncread(MFname,'HEIGHT_RHAT'); % height of Rh and T_a sensor
z_W    = ncread(MFname,'HEIGHT_WND'); % height of wind sensor
    
T_a  = ncread(MFname,'AIRT');
Rh   = ncread(MFname,'RELH'); % [%]
P_a  = ncread(MFname,'ATMS'); % [millibars]
Rl   = ncread(MFname,'LW');
Rs   = ncread(MFname,'SW');
Rain = ncread(MFname,'RAIN'); % [mm/hour]
W_e  = ncread(MFname,'UWND'); % east wind
W_n  = ncread(MFname,'VWND'); % north wind
T_s  = ncread(MFname,'TEMP'); % sea water temperature
S_s  = ncread(MFname,'PSAL'); % sea water salinity

%--------------------------------------------------------------------------
% lat_TS   = ncread(TSname,'LATITUDE');
% lon_TS   = ncread(TSname,'LONGITUDE');
time  = ncread(TSname,'TIME');

% get the reference time
t_ref = ncreadatt(TSname,'TIME','units'); % attribute 'units' for 'TIME'
t_ref = t_ref(12:end); % truncate to get the time string
t_ref = datenum(t_ref, 'yyyy-mm-ddTHH:MM:SSZ');
time_TS = t_ref + time; 
clear time t_ref

depth_TS = ncread(TSname,'DEPTH');
temp     = ncread(TSname,'TEMP');
salt     = ncread(TSname,'PSAL');

%--------------------------------------------------------------------------

% frequency bins for spectra
load([WDS_dir,'freq']); % 1st bin is for noise, last 3 are dummy 0's

WDSfiles = dir([WDS_dir,'/41061_*.mat']);
n_WDS    = length(WDSfiles);
n_f      = length(df);
g        = 9.81;

time_WDS = zeros(1,  n_WDS);
spec_h2  = zeros(n_f,n_WDS); 
xcmp     = zeros(n_f,n_WDS); % fraction for x-direction component
ycmp     = zeros(n_f,n_WDS); % fraction for y-direction component

for j = 1:n_WDS
    
    WDS = load([WDS_dir,WDSfiles(j).name]);
    
    time_WDS(j) = WDS.mday; % matlab datenumber
    
    % wave mean direction, switched to degree counter-clockwise from East,
    % pointing to where the waves are propagating toward
    wave_mdir   = 90 - (180 + WDS.alpha1); 
    xcmp(:,j)   = cosd(wave_mdir);
    ycmp(:,j)   = sind(wave_mdir);
    
%     dU_St = 16*pi^3/g * WDS.c11 .* WDS.R1 .* fctr.^3 .* xcmp(:,j) .* df;
%     dV_St = 16*pi^3/g * WDS.c11 .* WDS.R1 .* fctr.^3 .* ycmp(:,j) .* df;

    % wave spectrum input for GOTM
    spec_h2(:,j) = WDS.c11 .* WDS.R1 .* df;

end

fctr    = fctr(1:end-3);
xcmp    = xcmp(1:end-3,:);
ycmp    = ycmp(1:end-3,:);
spec_h2 = spec_h2(1:end-3,:);

%-- Note: all time in UTC

%% Calculation of fluxes

W_spd = sqrt(W_e.^2 + W_n.^2);

% call COARE (3.5) bulk formula
A = coare35vn(W_spd,z_W,T_a,z_RhTa,Rh,z_RhTa,P_a,T_s,Rs,...
              Rl,lat_MF,NaN,Rain,NaN,NaN);

tau = A(:,2); % wind stress [N/m^2]
hsb = A(:,3); % sensible heat flux into ocean [W/m^2]
hlb = A(:,4); % latent heat flux into ocan [W/m^2]
hbb = A(:,5); % surface buoyancy flux into ocean [W/m^2]

% seperate wind stress components
w_cos = W_e ./ W_spd;
w_sin = W_n ./ W_spd;
tau_x = tau .* w_cos;
tau_y = tau .* w_sin;

% compute net short wave heat flux into ocean
dVec_MF = datevec(time_MF);
yd_MF   = date2doy(time_MF) - 1;
nsw     = swhf(yd_MF,dVec_MF(:,1),-lon_MF*ones(size(time_MF)),...
               lat_MF*ones(size(time_MF)),Rs);
          % longitude: west is positive!

% compute net long wave heat flux into ocean
nlw = lwhf(T_s,Rl,Rs);

% net heatflux (+ in) = latent heat flux (+ in) + sensible heat flux (+ in)
%                       + net longwave radiation (+ in)

nhf = hsb + hlb + nlw;

%% TS profiles interpolation 

% [~,hr_inx] = ismember(time_MF,time_TS);
% 
% % get hourly profiles from original 5-mins data
% temp_hr = temp(:,hr_inx);
% salt_hr = salt(:,hr_inx);
% 
% % vertically interpolate to fill gaps
% temp_r = vert_fill(temp_hr,depth_TS);
% salt_r = vert_fill(salt_hr,depth_TS);

%% Project all data to WDS timeframe (shortest)

indx  = find(time_WDS <= time_MF(end),1,'last');
time  = time_WDS(1:indx);

tau_x = interp1(time_MF, tau_x, time);
tau_y = interp1(time_MF, tau_y, time);
nhf   = interp1(time_MF, nhf,   time);
nsw   = interp1(time_MF, nsw,   time);
Rain  = interp1(time_MF, Rain,  time);
T_s   = interp1(time_MF, T_s,   time);
S_s   = interp1(time_MF, S_s,   time);

% Not sure if initial conditions should have same timeframe as forcing,
% but I'll make them the same anyway.

% T-S profiles for writing out, in WDS timeframe
temp_o = zeros(length(depth_TS),length(time));
salt_o = zeros(length(depth_TS),length(time));

for i = 1:length(depth_TS)
    
    temp_o(i,:) = interp1(time_TS, temp(i,:), time);
    salt_o(i,:) = interp1(time_TS, salt(i,:), time);
end

% vertically interpolate to fill gaps
temp_o = vert_fill(temp_o,depth_TS,time);
salt_o = vert_fill(salt_o,depth_TS,time);

%% Write out to files

gotmdata_root = '~/Documents/GitHub/GOTM/gotmwork/data/';
basecase      = [gotmdata_root,'SPURSI/'];

date = string(datestr(time, 'yyyy-mm-dd HH:MM:SS'));

TAUname  = [basecase,'tau_file.dat'];
HFname   = [basecase,'heatflux_file.dat'];
RAINname = [basecase,'precip_file.dat'];
SSSname  = [basecase,'sss_file.dat'];
SSTname  = [basecase,'sst_file.dat'];
SWRname  = [basecase,'swr_file.dat'];
SPname   = [basecase,'sprof_file.dat'];
TPname   = [basecase,'tprof_file.dat'];
WDSname  = [basecase,'spec_file.dat']; % wave spectrum file

write_gotm_flux(TAUname,  [tau_x' tau_y'], date)
write_gotm_flux(HFname,   nhf',            date)
write_gotm_flux(RAINname, Rain',           date)
write_gotm_flux(SSSname,  S_s',            date)
write_gotm_flux(SSTname,  T_s',            date)
write_gotm_flux(SWRname,  nsw',            date)

temp_o3 = reshape(temp_o,length(depth_TS),1,length(time));
salt_o3 = reshape(salt_o,length(depth_TS),1,length(time));

write_gotm_ini(TPname,temp_o3,date,-depth_TS)
write_gotm_ini(SPname,salt_o3,date,-depth_TS)

waveInput  = cat(3,spec_h2,xcmp,ycmp);
waveInput3 = permute(waveInput,[1 3 2]);

write_gotm_ini(WDSname,waveInput3,date,fctr)

gzip(WDSname)
delete(WDSname)
