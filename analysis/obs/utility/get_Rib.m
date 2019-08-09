function Ribprof = get_Rib(Bprof,Ustar,z)
%
% get_Rib
%==========================================================================
%
% USAGE:
%  Ribprof = get_Rib(Bprof,Ustar,z)
%
% DESCRIPTION:
%  Compute nondimensional bulk Richardson number for different levels
%  Rib = N^2 / S^2 = [(Bs - Bh)/(zs - zh)]/[2*ustar/kappa/(zs + zh)]^2
%
% INPUT:
%
%  Bprof - [z,t] 2-D matrix of buoyancy profile time series
%  Ustar - [t] 1-D vector of waterside friction velocity [m/s]
%  z - [z] 1-D vector of vertical coordinates for buoyancy profile [-, m]
%     
% OUTPUT:
%
%  Ribprof - bulk Richardson number at different levels
%
% AUTHOR:
%  August 8 2019, Zhihua Zheng                            [ zhihua@uw.edu ]
%==========================================================================

kappa = 0.4;
nz    = length(z);
ntm   = length(Ustar);

%% law of wall approximation of velocity shear frequency squared [1/s^2]

mid_z = (z(1) + z(2:end))/2;
S2 = (repmat(Ustar,nz-1,1) ./ repmat(kappa*mid_z,1,ntm)).^2;

%% buoyancy frequency squared [1/s^2]

del_B = repmat(Bprof(1,:),nz-1,1) - Bprof(2:end,:);
del_z = repmat(z(1),nz-1,1) - z(2:end);

N2 = del_B ./ repmat(del_z,1,ntm);

%% Richardson number

Ribprof = N2 ./ S2;

end

