function [u_w, v_w, theta_w, s_w] = get_turb_flux(out, rotate_w)

% get_turb_flux
%==========================================================================
%
% USAGE:
%  [u_w, v_w, theta_w, s_w] = get_turb_flux(out, rotate_w)

%
% DESCRIPTION:
%  Compute the turbulent momentum flux and turbulent heat, salt flux using 
%  output from GOTM simulation
%
% INPUT:
%
%  out - A struct containing all the model output from GOTM
%  rotate_w - 1 or 0 (1 represents rotating current according to wind to 
%    get downwind and crosswind component)
% 
% OUTPUT:
%
%  u_w - turbulent x-momentum flux [m^2/s^2]
%  v_w - turbulent y-momentum flux [m^2/s^2]
%  theta_w - turbulent heat flux [K*m/s]
%  s_w - turbulent salt flux [PSU*m/s]
%
% AUTHOR:
%  Sep. 16 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%  Oct. 24 2018. Add turbulent salt flux

%% Note

% 1. cmue1 = sqrt(2)*Sm ... 
% 2. num = cmue1*sqrt(tke)*L = Sm*sqrt(2*tke)*L = Sm*q*L ...

%% Read relevant variables

% sacrifice data at both end of z-direction, staggered grid
nu_m = out.nu_m(2:end-1,:);
nu_cl = out.nu_cl(2:end-1,:);
nu_h = out.nu_h(2:end-1,:);
nu_s = out.nu_s(2:end-1,:);

%% Get vertical gradients
[u_z, v_z, uStokes_z, vStokes_z, temp_z, salt_z] = ...
    get_z_gradient(out,rotate_w);

%% Computation
u_w = -(nu_m.*u_z + nu_cl.*uStokes_z);
v_w = -(nu_m.*v_z + nu_cl.*vStokes_z);
theta_w = -(nu_h.*temp_z);
s_w = -(nu_s.*salt_z);

end
