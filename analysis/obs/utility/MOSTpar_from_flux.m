function [Lo,zeta,FS] = MOSTpar_from_flux(idatm,zz,casename,sld)
%
% MOSTpar_from_flux
%==========================================================================
%
% USAGE:
%  [Lo,zeta,FS] = MOSTpar_from_flux(idatm,z_inq,casename,sld)
%
% DESCRIPTION:
%  Compute parameters associated with Monin-Obukhov similarity theory for 
%  an inquired hourly period from a time series of hourly surface fluxes.
%
% INPUT:
%
%  idatm - 1-D column vector of MATLAB datetime for the inquired period
%  zz - vertical coordinates for the choosen 2 levels [-, m]
%  casename - string to indicate the parameters for which cases is inquired
%  sld - surface layer depth (optional) [+, m]
%         
% OUTPUT:
%
%  Lo - Monin-Obukhov lenght for evaluated depth [m]
%  zeta - Monin-Obukhov stability parameter
%  FS - struct contains turbulence scale for velocity (Ustar), temperature 
%       (Tstar), salinity (Sstar), buoyancy (Bstar), buoyancy forcing Bf 
%       and edge/interior indices for night time period (Inighte/Inighti)   
%
% AUTHOR:
%  June 29 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Parsing inputs

my_root = '~/GDrive/UW/Research/Data/';

% station longitude is west positive!

switch casename
    
    case 'Papa'
        moor_dir  = [my_root,'OCSP/Mooring/'];
        wave_dir  = [my_root,'OCSP/CDIP/'];
        Fname     = fullfile(moor_dir,'ocsp_flux_hrPMEL.mat');
        Wname     = fullfile(wave_dir,'ocsp_wave.mat');
        waterType = 4;
        lat       = 50.1;
        lon       = 144.9;
        
    case 'SPURSI'
        moor_dir  = [my_root,'SPURSI/Mooring/'];
        wave_dir  = [my_root,'SPURSI/DWS/'];
        Fname     = fullfile(moor_dir,'spursi_flux_hrUOP.mat');
%       Fname     = fullfile(moor_dir,'spursi_flux_hr.mat');
        Wname     = fullfile(wave_dir,'spursi_wave.mat');
        waterType = 1;
        lat       = 24.5811;
        lon       = 38;
end

M = load(Fname);

% dvec_inq = datevec(datestr(t_inq));
% yr  = dvec_inq(:,1);
% mon = dvec_inq(:,2);
% da  = dvec_inq(:,3);

%% Constants

rho0    = 1027; % reference density [kg/m^3]
% theta0  = 10;   % reference temperature [C]
% s0      = 33;   % reference salinity [g/kg]
% 
% alpha0 = 1.66e-4; % constant thermal expansion coefficient [1/C]
% beta0  = 7.6e-4;  % constant saline contraction coefficient [kg/g]
% cp0    = gsw_cp0; % constant specific heat of seawater [J/kg/C]
kappa  = 0.4;     % von Karman's constant
g      = 9.81;

%% Grab fluxes according to inquired time

ntm   = length(idatm);
index = isbetween(M.SF.datm,idatm(1),idatm(ntm));
SF    = M.SF(index,:);

% find night time indices based on solar radiation, excluding transition
% periods (edges)

mNSW   = mode(SF.nsw); % night time nsw
allN   = find(SF.nsw < mNSW + .3);
dallN  = diff(allN);
Nredge = [find(dallN>1); length(allN)]; % right edge indices
Nledge = [1;          find(dallN>1)+1]; % left  edge indices
Nedge  = union(Nredge,Nledge);
Ninter = setdiff(1:length(allN),Nedge);

FS.Inighte = ismember(1:ntm,allN(Nedge))';
FS.Inighti = ismember(1:ntm,allN(Ninter))';
FS.Idayi   = not(ismember(1:ntm,allN))';

% plot(SF.datm,SF.nsw);hold on;
% plot(SF.datm(FS.Idayi),SF.nsw(FS.Idayi),'.')
% plot(SF.datm(FS.Inighte),SF.nsw(FS.Inighte),'.');axis tight

% [rhr,rmin,shr,smin] = sunrise(mon,da,yr,lat,lon);
% srise = datetime(yr,mon,da,rhr-1,0,0);
% sset  = datetime(yr,mon,da,shr+1,0,0);

% season separation
idvec      = datevec(idatm);
FS.Iwinter = idvec(:,2) == 12 | idvec(:,2) == 1  | idvec(:,2) == 2;
FS.Ispring = idvec(:,2) == 3  | idvec(:,2) == 4  | idvec(:,2) == 5;
FS.Isummer = idvec(:,2) == 6  | idvec(:,2) == 7  | idvec(:,2) == 8;
FS.Iautumn = idvec(:,2) == 9  | idvec(:,2) == 10 | idvec(:,2) == 11;

% season vector
FS.Vseason             = zeros(ntm,1);
FS.Vseason(FS.Iautumn) = 1;
FS.Vseason(FS.Iwinter) = 2;
FS.Vseason(FS.Ispring) = 3;
FS.Vseason(FS.Isummer) = 4;

%% T-S dependent coefficients

% ignore the depth-dependence as T-S variation with depth is small

% sea surface SA and CT
ssSA = gsw_SA_from_SP(SF.sss,0,lon,lat); % [g/kg]
ssCT = gsw_CT_from_t(ssSA,SF.sst,0);

% isobaric heat capacity of seawater [J/kg/C]
cp = gsw_cp_t_exact(ssSA,SF.sst,0);

[~,alpha,beta] = gsw_specvol_alpha_beta(ssSA,ssCT,0);

%% Kinematic fluxes due to turbulent processes

% net surface heat flux Qtur
Qtur = SF.hlb + SF.hsb + SF.nlw; % [W/m^2]

% surface freshwater flux into the ocean
Ftur = (SF.rain - SF.evap)/1000/3600; % [m/s]

% surface kinematic fluxes
% w_u_0   = -SF.tau_x/rho0;
% w_v_0   = -SF.tau_y/rho0;
w_theta_0 = -Qtur./cp/rho0; % <w't'> [C*m/s]
w_s_0     =  Ftur.*ssSA;    % <w's'> [(g/kg)*(m/s)]
w_b_0     =  g*(alpha.*w_theta_0 - beta.*w_s_0); % <w'b'> [m^2/s^3]

%% Fluxes due to non-turbulent processes

z1      = zz(1); % upper level
z2      = zz(2); % lower level
zdum    = linspace(z2,z1)'; % dummy variable for integration
band_SR = 9; % solar radiation band model

switch band_SR
    
    case 0
    Iz = 0;
    
    case 2     % 2-band exponential decay of solar radiation
    Iz = get_SRz(SF.nsw,zdum,2,waterType);
    
    case 9     % 9-band exponential decay of solar radiation
    Iz = get_SRz(SF.nsw,zdum,9,waterType);
end

% obeject for array-vector division
avd = dsp.ArrayVectorDivider('Dimension',2);

% obeject for array-vector multiplication
avm = dsp.ArrayVectorMultiplier('Dimension',2);

% obeject for array-vector subtraction
avs = dsp.ArrayVectorSubtractor('Dimension',2);

% obeject for array-vector addition
ava = dsp.ArrayVectorAdder('Dimension',2);

w_theta_r = -avd( -avs(Iz,SF.nsw'), cp')/rho0;
w_b_r     =  g*avm(w_theta_r,alpha');

%% MOST parameters

% buoyancy forcing
FS.Bf = -ava(w_b_r,w_b_0');

% surface layer scales
FS.Ustar =  sqrt(SF.tau/rho0);
FS.Sstar = -w_s_0 ./ FS.Ustar / kappa;
FS.Tstar = -avd( ava(w_theta_r,w_theta_0'), FS.Ustar') / kappa;
FS.Bstar = -avd( ava(w_b_r,    w_b_0'    ), FS.Ustar') / kappa;

if exist('sld','var')
    
    W     = load(Wname);
    index = isbetween(W.SD.datm,idatm(1),idatm(ntm));
    SD    = W.SD(index,:);
    
    % surface layer averaged Stokes drift
    FS.USt_sl = get_St_SL(SD.uSt',SD.vSt',W.zSt,sld);
    
    % Obukhov-Langmuir length [m]
%     Lo = FS.Ustar.^2 .* FS.USt_sl' ./ Bf / kappa;
end

% Obukhov length [m]
Lo = avm(1./FS.Bf, FS.Ustar'.^3) / kappa;

% stability parameter
zeta = abs(repmat(zdum,1,ntm)) ./ Lo;

end


