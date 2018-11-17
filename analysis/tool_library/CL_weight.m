function qS = CL_weight(out, q, rotate_w)

% CL_weight
%==========================================================================
%
% USAGE:
%  qS = CL_weight(out, q, rotate_w)
%
% DESCRIPTION:
%  Compute the average of a quantity weighted by positive CL vortex 
%  production in the upper ocean. (see Harcourt 2015)
%
% INPUT:
%
%  out - A struct containing all the model output from GOTM
%  q - A matrix (t,z) containing the quantity to be weighted-averaged
%  rotata_w - 1 or 0 (1 represents rotating current according to wind to 
%    get downwind and crosswind component)
%
% OUTPUT:
%
%  qS - 1-D vector in temporal dimension. Averages of quantity q, weighted 
%    by positive CL vortex production in the upper ocean
%
% AUTHOR:
%  October 22 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%
%% Read relevant variables
z = mean(out.z,2);

%% Rotation of coordinate
if rotate_w
    
     new_vec = rotate_coor(out);

     u_stokes = new_vec.u_stokes;
     v_stokes = new_vec.v_stokes;
else
    u_stokes = out.u_stokes;
    v_stokes = out.v_stokes;
end

%% Stokes shear
uStokes_z = center_diff(u_stokes,z,1);
vStokes_z = center_diff(v_stokes,z,1);

%% Turbulence fluxes
[u_w, v_w, ~] = get_turb_flux(out,rotate_w);

%% Positive CL vortex production
PS = -(u_w.*uStokes_z + v_w.*vStokes_z); % CL vortex production
PS_plus = max(zeros(size(PS)),PS);

%% Layer thickness
h_new = diff(out.zi(2:end,:));

%% Weighted averges
qS = sum(q.*PS_plus.*h_new)./sum(PS_plus.*h_new);

end