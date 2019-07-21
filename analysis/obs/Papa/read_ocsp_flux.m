%% read_ocsp_flux

% Read processed air-sea fluxes for OCSP from PMEL website.
% Negative rainfall is included.

% Zhihua Zheng, UW-APL, July 15 2019

%% General setting

clear
flux_dir = '~/GDrive/UW/Research/Data/OCSP/Mooring/flux_PMEL/';
M_dir    = '~/GDrive/UW/Research/Data/OCSP/Mooring/L2009_2019/';

TAUname = fullfile(flux_dir,'tau50n145w_hr.cdf');
NSWname = fullfile(flux_dir,'swnet50n145w_hr.cdf');
NLWname = fullfile(flux_dir,'lwnet50n145w_hr.cdf');
QLHname = fullfile(flux_dir,'qlat50n145w_hr.cdf');
QSHname = fullfile(flux_dir,'qsen50n145w_hr.cdf');
Ename   = fullfile(flux_dir,'evap50n145w_hr.cdf');
Pname   = fullfile(flux_dir,'rain_wspd_cor50n145w_hr.cdf');
EMPname = fullfile(flux_dir,'emp50n145w_hr.cdf');
BFname  = fullfile(flux_dir,'bf50n145w_hr.cdf');
TSKname = fullfile(flux_dir,'tsk50n145w_hr.cdf');
POSname = fullfile(flux_dir,'pos50n145w_hr.cdf');

SSTname = fullfile(M_dir,'sst50n145w_hr.cdf');
SSSname = fullfile(M_dir,'sss50n145w_hr.cdf');

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

% emp   = ncread(EMPname,'EMP_251'); % [mm/hr]
% qemp  = ncread(EMPname,'QEMP_5251');
% bf    = ncread(BFname,'BF_191'); % [10^6*kg/(m^2*s)]
% qbf   = ncread(BFname,'QBF_5191');

sst   = ncread(SSTname,'T_25'); % [C]
qsst  = ncread(SSTname,'QT_5025');
tsst  = ncread(SSTname,'time');

sss   = ncread(SSSname,'S_41'); % [PSU]
qsss  = ncread(SSSname,'QS_5041');
tsss  = ncread(SSSname,'time');

%% Class

time  = double(time);
train = double(train);
tsst  = double(tsst);
tsss  = double(tsss);

tau_x = double(squeeze(tau_x));
tau_y = double(squeeze(tau_y));
tau   = double(squeeze(tau));
nsw   = double(squeeze(nsw));
nlw   = double(squeeze(nlw));
hlb   = double(squeeze(hlb));
hsb   = double(squeeze(hsb));
evap  = double(squeeze(evap));
rain  = double(squeeze(rain));
% emp   = double(squeeze(emp));
% bf    = double(squeeze(bf));
% tsk   = double(squeeze(tsk));
sst   = double(squeeze(sst));
sss   = double(squeeze(sss));

qtau  = double(squeeze(qtau));
qnsw  = double(squeeze(qnsw));
qnlw  = double(squeeze(qnlw));
qhlb  = double(squeeze(qhlb));
qhsb  = double(squeeze(qhsb));
qevap = double(squeeze(qevap));
qrain = double(squeeze(qrain));
% qemp  = double(squeeze(qemp));
% qbf   = double(squeeze(qbf));
qsst  = double(squeeze(qsst));
qsss  = double(squeeze(qsss));

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
% rain(rain<0)   = NaN;
% emp(qemp==0)   = NaN;
% bf(qbf==0)     = NaN;

sst(qsst==0 | qsst==5) = NaN;
sss(qsss==0 | qsss==5) = NaN;

%% Timestamp

time  = datenum(2007,6,8,5,0,0) + time/24;
train = datenum(2007,6,8,5,0,0) + train/24;
tsst  = datenum(2007,6,8,0,0,0) + tsst/24;
tsss  = datenum(2007,6,8,0,0,0) + tsss/24;

datm  = datetime(time,'ConvertFrom','datenum');
drain = datetime(train,'ConvertFrom','datenum');
dsst  = datetime(tsst,'ConvertFrom','datenum');
dsss  = datetime(tsss,'ConvertFrom','datenum');

sInx = find(dsst == datm(1));
eInx = find(dsst == datm(end));
sst  = sst(sInx:eInx);

sInx = find(dsss == datm(1));
sss  = sss(sInx:end);
sss  = [sss; nan(length(tau)-length(sss),1)];

sInx = find(drain == datm(1));
rain = rain(sInx:end);
rain = [rain; nan(length(tau)-length(rain),1)];

SF   = timetable(datm,tau_x,tau_y,tau,nsw,nlw,hlb,hsb,evap,rain,sst,sss);

%% Save

save('~/GDrive/UW/Research/Data/OCSP/Mooring/ocsp_flux_hrPMEL.mat','SF');
