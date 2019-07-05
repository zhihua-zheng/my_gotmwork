function [phiS,phiM] = get_theo_MOST(zeta)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here


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

