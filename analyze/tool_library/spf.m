function fS = spf(out, rotate_w)

% spf
%==========================================================================
%
% USAGE:
%  fS = = spf(out)
%
% DESCRIPTION:
%  Compute the surface proximity function in SMClT closure model 
%  (see Harcourt 2015)
%
% INPUT:
%
%  out - A struct containing all the model output from GOTM
%  rotata_w - 1 or 0 (1 represents rotating current according to wind to 
%    get downwind and crosswind component)
%
% OUTPUT:
%
%  fS - 2-D matrix (t,d) containing values for surface proximity function
%    used in SMCLT 
%
% AUTHOR:
%  October 21 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%

%% Read relevant variables
zi = out.zi(2:end-1,:);
L = out.L(2:end-1,:);
CS_l = 0.25;

%% Special length scale lS

% The length scale lS is an average of the dissipation lenght scale
% weighted by positive CL vortex production.

lS = CL_weight(out,L,rotate_w);

%% Computation
fS = 1 + tanh(CS_l*zi./repmat(lS,size(zi,1),1));

end



