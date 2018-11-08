function tke_comps = get_tke_comp(model_par, out, rotate_w)

% get_tke_comp
%==========================================================================
%
% USAGE:
%  tke_comps = get_tke_comp(model_par, out, rotate_w)

%
% DESCRIPTION:
%  Compute different components of turbulent kinetic energy (TKE) based on
%  given information of the turbulent fluxes, length scale, TKE magnitude
%  and SMC model setting.
%
% INPUT:
%
%  model_par - A struct containing parameters used in the model
%  out - A struct containing all the model output from GOTM
%  rotata_w - 1 or 0 (1 represents rotating current according to wind to 
%    get downwind and crosswind component)
%
% OUTPUT:
%
%  tke_comps - 3-D matrix containning different components of TKE
%  tke_comps(:,:,1) - horizontal(x) velocity fluctuation variance [m^2/s^2]
%  tke_comps(:,:,2) - horizontal(y) velocity fluctuation variance [m^2/s^2]
%  tke_comps(:,:,3) - vertical(z) velocity fluctuation variance [m^2/s^2]
%
% AUTHOR:
%  September 16 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%

%% Note 

% 1. After applying rescale of l/q according to Ozmidov scale, some l/q is
%    still too large, causing negative vTKE, tracing back to very small
%    buoyancy frequency.
% 2. Not all the negative values of vTKE are due to l/q > 0.53/N. Found 3
%    points except initial column. index = (1223, 1224, 8447). the middle
%    one is the minimum of w_w.

%% Read relevant variables

% Z = out.z;
% T = repmat(out.time',size(Z,1),1);
% NN = out.NN; % buoyancy frequency

% sacrifice data at both end of z-direction, staggered grid

% Zi = out.zi(2:end-1,:);
% Ti = repmat(out.time',size(Zi,1),1);
L = out.L(2:end-1,:);
q2 = 2*out.tke(2:end-1,:); 
q = sqrt(q2); % turbulent velocity scale [m/s]

rescale_r = model_par.rescale_r;
A1 = model_par.A1;
B1 = model_par.B1;
g = 9.81;

% divided by (-rho_0) to get thermal expanison coefficient (positive)
alpha = - model_par.dtr0/model_par.rho_0; 

% divided by (-rho_0) to get saline expanison coefficient (positive)
belta =  model_par.dsr0/model_par.rho_0; 

%% Get vertical gradients
[u_z, v_z, uStokes_z, vStokes_z, temp_z, salt_z] = ...
    get_z_gradient(out,rotate_w);

%% Turbulence fluxes
[u_w, v_w, theta_w, s_w] = get_turb_flux(out,rotate_w);

%% Rescale l/q under stable stratification

% avoid using the NN output, the grid level is not aligned with L ...
% compute the NN from temperature and salinity gradient instead.

NN = g*(alpha*temp_z - belta*salt_z);

NN(NN<0) = NaN;
N = sqrt(NN);
% N = interp2(T,Z,N,Ti,Zi,'linear');

l_over_q = L./q; % time scale for turbulence
r_Ozm = 0.53./N; % time scale cooresponding to Ozmidov length scale 

% replace the large values of l/q with r_Ozm
if rescale_r
    l_over_q(l_over_q > r_Ozm) = r_Ozm(l_over_q > r_Ozm);
end

%% Surface proximity function
fS = spf(out,rotate_w);

%% Angle of horizontal unit vector
phi = angle_h15(out,rotate_w);
phi = repmat(phi,size(fS,1),1);

%% Computation

% TO-DO: incorporate turbulent salinity fluxes

% u_u
tke_comps(:,:,1) = q2*(1-6*A1/B1)/3 - (6*A1*l_over_q).*(u_w.*u_z + ...
    (sin(phi)).^2.*fS.*(u_w.*uStokes_z + v_w.*vStokes_z));

% v_v
tke_comps(:,:,2) = q2*(1-6*A1/B1)/3 - (6*A1*l_over_q).*(v_w.*v_z + ...
    (cos(phi)).^2.*fS.*(u_w.*uStokes_z + v_w.*vStokes_z));

% w_w
tke_comps(:,:,3) = q2*(1-6*A1/B1)/3 + (6*A1*l_over_q).*(alpha*g*theta_w...
    - belta*g*s_w + (fS-1).*(u_w.*uStokes_z + v_w.*vStokes_z));


end