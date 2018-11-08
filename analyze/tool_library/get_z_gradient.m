function [u_z, v_z, uStokes_z, vStokes_z, temp_z, salt_z] = ...
    get_z_gradient(out, rotate_w)

% get_z_gradient
%==========================================================================
%
% USAGE:
%  [u_z, v_z, uStokes_z, vStokes_z, temp_z, salt_z] = ...
%    get_z_gradient(out, rotate_w)

%
% DESCRIPTION:
%  Compute the vertical gradient of velocity, Stokes drift, temperature and
%  salinity using the output from GOTM simulation
%
% INPUT:
%
%  out - A struct containing all the model output from GOTM
%  rotate_w - 1 or 0 (1 represents rotating current according to wind to 
%    get downwind and crosswind component)
% 
% OUTPUT:
%
%  u_z - vertical gradient of x-direction velocity [1/s]
%  v_z - vertical gradient of y-direction velocity [1/s]
%  uStokes_z - vertical gradient of x-direction Stokes drift [1/s]
%  vStokes_z - vertical gradient of y-direction Stokes drift [1/s]
%  temp_z - vertical gradient of temperature [K/m]
%  salt_z - vertical gradient of salinity [PSU/m]
%
% AUTHOR:
%  Oct. 24 2018. Zhihua Zheng                       [ zhihua@uw.edu ]

%% TO-DO

% 1. check if the velocity and temprature output are averaged quantities
%    - Using these output as averaged quantities is ok!

%% Read relevant variables

z = mean(out.z,2);
temp = out.temp;
salt = out.salt;

%% Rotation of coordinate
if rotate_w
    
     new_vec = rotate_coor(out);
     
     u = new_vec.u;
     v = new_vec.v;
     u_stokes = new_vec.u_stokes;
     v_stokes = new_vec.v_stokes;
else
    u = out.u;
    v = out.v;
    u_stokes = out.u_stokes;
    v_stokes = out.v_stokes;
end

%% Eulerian shear
u_z = center_diff(u,z,1);
v_z = center_diff(v,z,1);

%% Stokes shear
uStokes_z = center_diff(u_stokes,z,1);
vStokes_z = center_diff(v_stokes,z,1);

%% Temperature & salinity gradient
temp_z = center_diff(temp,z,1);
salt_z = center_diff(salt,z,1);

end
