%% read_spursi_wave
%
% Script to read directional wave spectrum data from SPURS-I site
% MAT file is saved as 'spursi_wave.mat'.

% Zhihua Zheng, UW-APL, July 30 2019

%% General setting

clear
data_dir = '~/GDrive/UW/Research/Data/SPURSI/DWS/';

%% Loading

% frequency bins for spectra
% 1st bin is for noise, last 3 are dummy 0's
load([data_dir,'freq'],'fctr','wave_bw');

DWSfiles = dir([data_dir,'/41061_*.mat']);
n_DWS    = length(DWSfiles);
n_f      = length(wave_bw);

wave_time = zeros(n_DWS,1);
wave_spec = zeros(n_DWS,n_f); % spectral energy density (m^2/hertz)
wave_a1   = zeros(n_DWS,n_f);
wave_b1   = zeros(n_DWS,n_f);
wave_mdir = zeros(n_DWS,n_f);

for i = 1:n_DWS
    
    DWS = load([data_dir,DWSfiles(i).name]);
    
    wave_time(i) = DWS.mday; % matlab datenumber
    
    % wave mean direction, clockwise from North [degree, from]
    wave_mdir(i,:) = 90 - (180 + DWS.alpha1);
    
%     xcmp(:,j)   = cosd(wave_mdir);
%     ycmp(:,j)   = sind(wave_mdir);   
%     dU_St = 16*pi^3/g * DWS.c11 .* DWS.R1 .* fctr.^3 .* xcmp(:,j) .* df;
%     dV_St = 16*pi^3/g * DWS.c11 .* DWS.R1 .* fctr.^3 .* ycmp(:,j) .* df;

    wave_spec(i,:) = DWS.c11;
    wave_a1(i,:)   = DWS.R1 .* cosd(DWS.alpha1);
    wave_b1(i,:)   = DWS.R1 .* sind(DWS.alpha1);
end

fctr      = fctr(1:end-3);
wave_bw   = wave_bw(1:end-3);

wave_spec = wave_spec(:,1:end-3);
wave_a1   = wave_a1(:,1:end-3);
wave_b1   = wave_b1(:,1:end-3);
wave_mdir = wave_mdir(:,1:end-3);

datm = datetime(wave_time,'ConvertFrom','datenum');

%% compute Stokes drift profiles

wq  = timetable(datm,wave_spec,wave_a1,wave_b1);
zSt = (-30:.5:0)';

nz  = length(zSt);
ntm = length(datm);

[uSt,vSt] = St_from_dws(wq,fctr,wave_bw,zSt);

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

save('~/GDrive/UW/Research/Data/SPURSI/DWS/spursi_wave.mat','SD','zSt');

