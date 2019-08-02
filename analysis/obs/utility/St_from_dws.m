function [uSt,vSt] = St_from_dws(WQ,fctr,bw,zSt)
%
% St_from_dws
%==========================================================================
%
% USAGE:
%  [uSt,vSt] = St_from_dws(WQ,fctr,bw,zSt)
%
% DESCRIPTION:
%  Compute Stokes drift profile from directional wave spectrum
%
% INPUT:
%
%  WQ - 
%  fctr - 
%  bw - 
%  zSt - 
%     
% OUTPUT:
%
%  uSt - Stokes drift profile in x direction
%  vSt - Stokes drift profile in  ydirection
%
% AUTHOR:
%  July 31 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Loading

datm = WQ.datm;
spec = WQ.wave_spec';
a1   = WQ.wave_a1';
b1   = WQ.wave_b1';

%% Constants

g     = 9.81;
multi = 16*pi^3/g;
decoe = 8*pi^2/g;
fc    = fctr(end) + bw(end)/2; % right edge cutoff frequency [Hz]

%% Reshaping

df    = repmat(bw,  1,length(datm),length(zSt));
f     = repmat(fctr,1,length(datm),length(zSt));
z     = repmat(zSt, 1,length(datm),length(bw));
z     = permute(z,[3 2 1]);

Spec  = repmat(spec,1,1,length(zSt));
A1    = repmat(a1,  1,1,length(zSt));
B1    = repmat(b1,  1,1,length(zSt));

%% Resolved Stokes spectrum

decay = exp(decoe*z .* f.^2);

uSt_re = squeeze(-multi*sum(f.^3 .* Spec .* B1 .* decay .* df))';
vSt_re = squeeze(-multi*sum(f.^3 .* Spec .* A1 .* decay .* df))';

%% Append high frequency tail (Breivik et al. 2014)

mu    = -decoe * zSt;
tailC =  multi * fc^4 * (exp(-mu*fc^2) - ...
         fc * sqrt(mu*pi) .* (1 - erf(fc*sqrt(mu))) );

uSt_tail = -tailC * (spec(end,:) .* b1(end,:));
vSt_tail = -tailC * (spec(end,:) .* a1(end,:));

%% Add up

uSt = uSt_re + uSt_tail;
vSt = vSt_re + vSt_tail;

end

