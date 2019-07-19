%% read_ocsp_flux

% Read processed air-sea fluxes for OCSP from PMEL website.
% Negative rainfall is treated as missing data.

% Zhihua Zheng, UW-APL, July 15 2019

%% General setting

clear
data_dir = '~/GDrive/UW/Research/Data/OCSP/Mooring/flux_PMEL/';

TAUname = fullfile(data_dir,'tau50n145w_hr.cdf');
NSWname = fullfile(data_dir,'swnet50n145w_hr.cdf');
NLWname = fullfile(data_dir,'lwnet50n145w_hr.cdf');
QLHname = fullfile(data_dir,'qlat50n145w_hr.cdf');
QSHname = fullfile(data_dir,'qsen50n145w_hr.cdf');
Ename   = fullfile(data_dir,'evap50n145w_hr.cdf');
Pname   = fullfile(data_dir,'rain_wspd_cor50n145w_hr.cdf');
POSname = fullfile(data_dir,'pos50n145w_hr.cdf');

%% Read variables

time  = ncread(TAUname,'time');
tau_x = ncread(TAUname,'TX_442');  % [N/m^2]
tau_y = ncread(TAUname,'TY_443');  % [N/m^2]
tau   = ncread(TAUname,'TAU_440'); % [N/m^2]
qtau  = ncread(TAUname,'QTAU_5440');

nsw   = ncread(NSWname,'SWN_1495'); % [W/m^2] into the ocean
qnsw  = ncread(NSWname,'QSW_5495');

nlw   = -ncread(NLWname,'LWN_1136'); % [W/m^2] into the ocean
qnlw  =  ncread(NLWname,'QLW_5136');

hlb   = -ncread(QLHname,'QL_137');   % [W/m^2] into the ocean
qhlb  =  ncread(QLHname,'QQL_5137');

hsb   = -ncread(QSHname,'QS_138');   % [W/m^2] into the ocean
qhsb  =  ncread(QSHname,'QQS_5138');

evap  = ncread(Ename,'E_250');  % [mm/hr]
qevap = ncread(Ename,'QE_5250');

rain  = ncread(Pname,'RN_485'); % [mm/hr]
qrain = ncread(Pname,'QRN_5485');
train = ncread(Pname,'time');

%% Class

time  = double(time);
train = double(train);
tau_x = double(squeeze(tau_x));
tau_y = double(squeeze(tau_y));
tau   = double(squeeze(tau));
nsw   = double(squeeze(nsw));
nlw   = double(squeeze(nlw));
hlb   = double(squeeze(hlb));
hsb   = double(squeeze(hsb));
evap  = double(squeeze(evap));
rain  = double(squeeze(rain));

qtau  = double(squeeze(qtau));
qnsw  = double(squeeze(qnsw));
qnlw  = double(squeeze(qnlw));
qhlb  = double(squeeze(qhlb));
qhsb  = double(squeeze(qhsb));
qevap = double(squeeze(qevap));
qrain = double(squeeze(qrain));

%% Pre-processing

tau_x(qtau==0) = NaN;
tau_y(qtau==0) = NaN;
tau(qtau==0)   = NaN;
nsw(qnsw==0)   = NaN;
nlw(qnlw==0)   = NaN;
hlb(qhlb==0)   = NaN;
hsb(qhsb==0)   = NaN;
evap(qevap==0) = NaN;
rain(qrain==0) = NaN;
rain(rain<0)   = NaN;

%% Timestamp

time  = datenum(2007,6,8,5,0,0) + time/24;
train = datenum(2007,6,8,5,0,0) + train/24;
datm  = datetime(time,'ConvertFrom','datenum');
drain = datetime(train,'ConvertFrom','datenum');

SF1    = timetable(datm,tau_x,tau_y,tau,nsw,nlw,hlb,hsb,evap);
Precip = timetable(drain,rain);
SF     = synchronize(SF1,Precip);

%% Save

save('~/GDrive/UW/Research/Data/OCSP/Mooring/ocsp_flux_hrPMEL.mat','SF');
