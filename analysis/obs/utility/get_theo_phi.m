function [phiS,phiM] = get_theo_phi(zeta)
%
% get_theo_phi
%==========================================================================
%
% USAGE:
%  [phiS,phiM] = get_theo_phi(zeta)
%
% DESCRIPTION:
%  Compute empirical Monin-Obukhov similarity functions for scalars and
%  monetum according to the form in Lerge et al. 1994
%
% INPUT:
%
%  zeta - Monin-Obukhov stability parameter = |z|/Lo
%     
% OUTPUT:
%
%  phiS - Empirical similarity function for scalars
%  phiM - Empirical similarity function for momentum
%
% AUTHOR:
%  July 19 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Constants

zetaS = -1;
cS    =  98.96;
aS    = -28.86;

zetaM = -0.2;
cM    =  8.38;
aM    =  1.26;

%% Theoretical dimensionaless functions (Large et al. 1994)

phiS = nan(size(zeta));
phiM = nan(size(zeta));

iStable     = find(zeta >= 0);
iWUnstableS = find(zeta < 0 & zeta >= zetaS);
iWUnstableM = find(zeta < 0 & zeta >= zetaM);
iSUnstableS = find(zeta <= zetaS);
iSUnstableM = find(zeta <= zetaM);

phiS(iStable) = 1 + 5*zeta(iStable);
phiM(iStable) = 1 + 5*zeta(iStable);

phiS(iWUnstableS) = (1 - 16*zeta(iWUnstableS)).^(-1/2);
phiM(iWUnstableM) = (1 - 16*zeta(iWUnstableM)).^(-1/4);

phiS(iSUnstableS) = (aS - cS*zeta(iSUnstableS)).^(-1/3);
phiM(iSUnstableM) = (aM - cM*zeta(iSUnstableM)).^(-1/3);

end

