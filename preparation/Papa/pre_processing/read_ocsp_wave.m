
%% read_ocsp_wave
%
% Script to read directional wave spectrum data from OCSP
% MAT file is saved as 'ocsp_wave.mat'.

% Zhihua Zheng, UW-APL, August 13 2018

%% General setting

clear

% NDBC buoy collects pitch, roll, and heave displacements of the buoy
% CDIP waverider uses Lagrangian translation

cdip_dir = '~/GDrive/UW/Research/Data/OCSP/CDIP/';
NCname     = fullfile(cdip_dir,'166p1_historic.nc');

%% read variables 

wave_time = double(ncread(NCname,'waveTime')); % wave time in UTC
fctr      = double(ncread(NCname,'waveFrequency')); % band center frequency (hertz)
wave_bw   = ncread(NCname,'waveBandwidth'); % frequency bandwidth (hertz)

% spectral energy density (m^2/hertz)
wave_spec = (ncread(NCname,'waveEnergyDensity'))'; % equivalent to C_11

% band a1, b1 angular Fourier coefficients
wave_a1 = (ncread(NCname,'waveA1Value'))';
wave_b1 = (ncread(NCname,'waveB1Value'))';

% Directional spreading function at frequency f:
% D(theta) = 1/pi( 1/2 + sum{ An*cos(n*theta) + Bn*sin(n*theta) } )

% band mean direction (from), in degree clockwise from true North
wave_mdir = (ncread(NCname,'waveMeanDirection'))'; 
% equivalent to theta_1 = atan2(wave_b1,wave_a1)*180/pi + 360; (verified!)

wave_Hs = ncread(NCname,'waveHs'); % [m]

%% timestamp

% get the reference time
t_ref = ncreadatt(NCname,'waveTime','units'); % get the attribute 'units' for 'time'
t_ref = t_ref(15:end); % truncate to get the time string
t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS'); 

wave_time = double(t_ref + wave_time./3600/24);
datm = datetime(wave_time,'ConvertFrom','datenum');

clear t_ref

%% direction conversion

wave_mdir = 90 - (180 + wave_mdir);
% henceforth direction (towards) is in degree counter-clockwise from East

%% compute Stokes drift profiles

wq  = timetable(datm,wave_spec,wave_a1,wave_b1);
zSt = (-30:.5:0)';

nz  = length(zSt);
ntm = length(datm);
uSt = nan(nz,ntm);
vSt = nan(nz,ntm);

% split the array due to MATLAB size limit
jt = (0:13)'*10000;
jt(1) = 1;
jt(end) = ntm + 1;

for j = 1:13
    
    jl = jt(j);
    jr = jt(j+1) - 1;
    [uSt(:,jl:jr),vSt(:,jl:jr)] = St_from_dws(wq(jl:jr,:),fctr,wave_bw,zSt);
end

%% timetable

uSt = uSt';
vSt = vSt';
sd  = timetable(datm,uSt,vSt);

tlower  = dateshift(datm(1),'start','hour');
tupper  = dateshift(datm(end),'end','hour');
datm_hr = (tlower:hours(1):tupper)';
datmE   = datm_hr - minutes(30);

SD = retime(sd,datmE,'mean'); % hourly average
SD.datm = datm_hr;

%% save

save('~/GDrive/UW/Research/Data/OCSP/CDIP/ocsp_wave.mat','SD','zSt');

