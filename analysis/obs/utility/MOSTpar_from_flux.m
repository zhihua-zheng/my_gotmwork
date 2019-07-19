function [Lo,Bf,FS] = MOSTpar_from_flux(idatm,z_inq,casename,sld)
%
% MOSTpar_from_flux
%==========================================================================
%
% USAGE:
%  [Lo,Bf,FS] = MOSTpar_from_flux(idatm,z_inq,casename,sld)
%
% DESCRIPTION:
%  Compute parameters associated with Monin-Obukhov similarity theory for 
%  an inquired hourly period from a time series of hourly surface fluxes.
%
% INPUT:
%
%  idatm - 1-D column vector of MATLAB datetime for the inquired period
%  z_inq - vertical coordinates for evaluation depth  [-, m]
%  casename - string to indicate the parameters for which cases is inquired
%  sld - 
%         
% OUTPUT:
%
%  Lo - Monin-Obukhov lenght for evaluated depth [m]
%  Bf - surface buoyancy forcing [m^2/s^3]
%  FS - struct contains turbulence scale for velocity (Ustar), temperature 
%       (Tstar), salinity (Sstar), and edge/interior indices for night time
%       period (Inighte/Inighti)
%
% AUTHOR:
%  Jul1 29 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Parsing inputs

my_root = '~/GDrive/UW/Research/Data/';

% station longitude is west positive!

switch casename
    
    case 'Papa'
        moor_dir  = [my_root,'OCSP/Mooring/'];
        wave_dir  = [my_root,'OCSP/CDIP/'];
        Fname     = fullfile(moor_dir,'PMEL_flux_hr.mat');
        DWSname   = fullfile(wave_dir,'Wave_2010.mat');
        waterType = 4;
%         lat       = 50.1*ones(size(idatm));
%         lon       = 144.9*ones(size(idatm));
        
    case 'SPURSI'
        moor_dir  = [my_root,'SPURSI/Mooring/'];
        Fname     = fullfile(moor_dir,'spursi_flux_hr.mat');
        waterType = 1;
end

M = load(Fname);

% datm_inq = datetime(t_inq,'ConvertFrom','datenum');
% dvec_inq = datevec(datestr(t_inq));
% yr  = dvec_inq(:,1);
% mon = dvec_inq(:,2);
% da  = dvec_inq(:,3);

%% Constants

rho0    = 1027; % reference density [kg/m^3]
theta0  = 10;   % reference temperature [C]
s0      = 33;   % reference salinity [g/kg]

alpha0 = 1.66e-4; % thermal expansion coefficient [1/C]
beta0  = 7.6e-4;  % saline contraction coefficient [kg/g]
cp0    = gsw_cp0; % specific heat of seawater [J/kg/C]
kappa  = 0.4;     % von Karman's constant
g      = 9.81;

%% Grab fluxes according to inquired time

index = M.SF.datm >= idatm(1) & M.SF.datm <= idatm(end);
SF    = M.SF(index,:);
% SF    = retime(SF,'regular','mean','TimeStep',hours(3));

% find night time indices based on solar radiation, excluding transition
% periods (edges)

mNSW   = mode(SF.nsw);
allN   = find(SF.nsw < mNSW + .3);
dallN  = diff(allN);
Nredge = [find(dallN>1); length(allN)]; % right edge indices
Nledge = [1;          find(dallN>1)+1]; % left  edge indices
Nedge  = union(Nredge,Nledge);
Ninter = setdiff(1:length(allN),Nedge);

FS.Inighte = allN(Nedge);
FS.Inighti = allN(Ninter);

% plot(SF.datm,SF.nsw);hold on; plot(SF.datm(FS.Inighte),SF.nsw(FS.Inighte),'.')

% [rhr,rmin,shr,smin] = sunrise(mon,da,yr,lat,lon);
% srise = datetime(yr,mon,da,rhr-1,0,0);
% sset  = datetime(yr,mon,da,shr+1,0,0);

%% Kinematic fluxes due to turbulent processes

% net surface heat flux Qtur
Qtur = SF.hlb + SF.hsb + SF.nlw; % [W/m^2]

% surface freshwater flux
% SF.rain(isnan(SF.rain)) = 0; % NaNs will dramatically decrease the amount of useful data
Ftur = (SF.rain - SF.evap)/1000/3600; % [m/s]

% surface kinematic fluxes
w_u_0     = -SF.tau_x/rho0;
w_v_0     = -SF.tau_y/rho0;
w_theta_0 = -Qtur/(rho0*cp0); % <w't'> [C*m/s]
w_s_0     =  Ftur.*s0; % <w's'> [(g/kg)*(m/s)]
w_b_0     =  g*(alpha0.*w_theta_0 - beta0.*w_s_0);

%% Fluxes due to non-turbulent processes

depth_SR = 2;

switch depth_SR
    
    case 0
    Iz = 0;
    
    case 2     % 2-band exponential decay of solar radiation
    Iz = get_SRz(SF.nsw,z_inq,2,waterType);
    
    case 9     % 9-band exponential decay of solar radiation
    Iz = get_SRz(SF.nsw,z_inq,9,waterType);
end

w_theta_r = -(SF.nsw - Iz)/(rho0*cp0);
w_b_r     = g*alpha0*w_theta_r;

%% MOST parameters

% surface buoyancy forcing
Bf = -(w_b_0 + w_b_r);

% friction scales
FS.Ustar =  sqrt(SF.tau/rho0);
FS.Tstar = -(w_theta_0 + w_theta_r) ./ FS.Ustar / kappa;
FS.Sstar =  w_s_0 ./ FS.Ustar / kappa;

if exist('sld','var')
    
    DWS     = load(DWSname);
    WQ      = retime(DWS.WQ,datm_inq,'linear');
    zSt     = (-30:.2:0)';
    
    [uSt,vSt] = St_from_dws(WQ,DWS.wave_freq,DWS.wave_bw,zSt);
    FS.SD_sl  = get_St_SL(uSt,vSt,zSt,sld);
    
    Lo = FS.Ustar.^2 .* FS.SD_sl' ./ Bf / kappa; % Obukhov-Langmuir length [m]
else
    Lo = FS.Ustar.^3 ./ Bf / kappa; % Obukhov length [m]
end

end


