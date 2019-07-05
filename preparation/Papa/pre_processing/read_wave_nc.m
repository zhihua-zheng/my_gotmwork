%% read_wave_nc.m
%
% This is the script used to read the information from wave spectrum netCDF
% file and the generated workspace is named as 'Wave_2010.mat'.

% by Zhihua Zheng (UW/APL), updated on Aug. 13 2018

%% General setting

clear

% NDBC buoy collects pitch, roll, and heave displacements of the buoy
% CDIP waverider uses Lagrangian translation

cdip_dir = '~/GDrive/UW/Research/Data/OCSP/CDIP/';
NCname     = fullfile(cdip_dir,'166p1_historic.nc');

%% read variables 

wave_time = double(ncread(NCname,'waveTime')); % wave time in UTC
wave_freq = double(ncread(NCname,'waveFrequency')); % band center frequency (hertz)
wave_bw   = ncread(NCname,'waveBandwidth'); % frequency bandwidth (hertz)

% spectral energy density (m^2/hertz)
wave_spec = ncread(NCname,'waveEnergyDensity'); % equivalent to C_11

% band a1, b1 angular Fourier coefficients
wave_a1 = ncread(NCname,'waveA1Value');
wave_b1 = ncread(NCname,'waveB1Value');

% Directional spreading function at frequency f:
% D(theta) = 1/pi( 1/2 + sum{ An*cos(n*theta) + Bn*sin(n*theta) } )

% band mean direction (from), in degree clockwise from true North
wave_mdir = ncread(NCname,'waveMeanDirection'); 
% equivalent to theta_1 = atan2(wave_b1,wave_a1)*180/pi + 360; (verified!)

wave_Hs = ncread(NCname,'waveHs'); % [m]

%% timestamp

% get the reference time
t_ref = ncreadatt(NCname,'waveTime','units'); % get the attribute 'units' for 'time'
t_ref = t_ref(15:end); % truncate to get the time string
t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS'); 

wave_time = double(t_ref + wave_time./3600/24);
wave_date = string(datestr(wave_time,'yyyy-mm-dd HH:MM:SS'));

clear t_ref

%% direction conversion

wave_mdir = 90 - (180 + wave_mdir);
% henceforth direction (towards) is in degree counter-clockwise from East

%% save

clear cdip_dir NCname
save('~/GDrive/UW/Research/Data/OCSP/CDIP/Wave_2010.mat');

