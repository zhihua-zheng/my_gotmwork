function Rig = get_Rig(Bzz,Ustar,zz)
%
% get_Rig
%==========================================================================
%
% USAGE:
%  Rig get_Rig(Bprof,Ustar,zz)
%
% DESCRIPTION:
%  Compute nondimensional gradient Richardson number at middle level
%  Rig = N^2 / S^2
%
% INPUT:
%
%  Bzz - 2-D matrix containing buoyancy values at 2 different levels
%  Ustar - 1-D vector containing waterside friction velocity [m/s]
%  zz - vertical coordinates for the choosen 2 levels [-, m]
%     
% OUTPUT:
%
%  Rig - gradient Richardson number at middle level
%
% AUTHOR:
%  July 25 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

kappa = 0.4;
zmid  = mean(zz);

% law of wall approximation of velocity shear frequency squared [1/s^2]
S2 = (Ustar ./ (kappa*abs(zmid))).^2;

% buoyancy frequency squared [1/s^2]
N2 = center_diff(Bzz,zz,1,'mid');

Rig = N2 ./ S2;

end

