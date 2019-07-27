function IntPhiX = get_MOST_Delta(Xstar,zeta,Lo,zz)
%
% get_MOST_Delta
%==========================================================================
%
% USAGE:
%  IntPhi = get_MOST_Delta(Fstar,)
%
% DESCRIPTION:
%  Use the integral approach to estimate the difference of mean quantities
%  at different levels based on traditional Monin-Obukhov scaling.
%
% INPUT:
%
%  Xstar - the friction scale for the chosen quantity
%  zeta - Monin-Obukhov stability parameter
%  Lo - Monin-Obukhov length for evaluated depth [m]
%  zz - vertical coordinates for the choosen 2 levels [-, m]
%         
% OUTPUT:
%
%  IntPhiX - difference of quantity at two different levels
%
% AUTHOR:
%  July 24 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

ntm  = size(Lo,2);

z1   = zz(1);
z2   = zz(2);
zdum = linspace(z2,z1)'; % dummy variable for integration
Zdum = abs(repmat(zdum,1,ntm));

% compute empirical similarity function (scalar)
[phiXdum,~] = get_theo_phi(zeta);

intgrd  = phiXdum .* Xstar ./ Zdum;
IntPhiX = trapz(zdum,intgrd);

end