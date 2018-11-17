function mld = get_mld(rho, z)

% get_mld
%==========================================================================
%
% USAGE:
%  mld = get_mld(rho, z)
%
%
% DESCRIPTION:
%  Compute the mixed layer depth from the density profile
%
% INPUT:
%
%  rho - 2-D matrix (t,z) containing density profiles
%  z - 1-D vector containing vertical coordinates for density profile
% 
% OUTPUT:
%
%  mld - mixed layer depth [m]
%
%
% AUTHOR:
%  Oct. 29 2018. Zhihua Zheng                       [ zhihua@uw.edu ]

sigma = rho - 1000;
d_sig = 1.01*sigma(end,:);  %% 0.01 criteria
mld = zeros(size(d_sig'));

where_mld = sigma - repmat(d_sig,length(z),1);

for j = 1:length(d_sig)
    
        % last positive value & first negative value
        mld_down = find(where_mld(:,j)>0,1,'last');
        mld_up = find(where_mld(:,j)<0,1,'first');
        
        % interpolate to find where the zero is
        mld(j) = -interp1(where_mld([mld_down,mld_up],j),...
            z([mld_down,mld_up]),0);
end

end
