function mld = get_mld(A, z, flag)
%
% get_mld
%==========================================================================
%
% USAGE:
%  mld = get_mld(A, z, flag)
%
% DESCRIPTION:
%  Compute the mixed layer depth from the density/temperature profile
%
% INPUT:
%
%  A - 2-D matrix (z,t) containing density/temperature profiles
%  z - 1-D column vector with vertical coordinates for density/temperature 
%      profile, bottom to top [-, m]
%  flag - 1 or 2 to indicate using density criteria or temperature criteria
% 
% OUTPUT:
%
%  mld - mixed layer depth [+, m]
%
% AUTHOR:
%  Oct. 29 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%==========================================================================

switch flag
    case 1
        
    sigma = A - 1000;

    % Get reference values
    sigma_10 = zeros(1,size(sigma,2));

    for j = 1:size(sigma,2)
        sigma_10(j) = interp1(z,sigma(:,j),-10);
    end

    % mld is where the density is 0.03 kg/m^3 higher than
    % a surface reference value at 10 m (see de Boyer Montégut et al. 2004)
    d_sig = sigma_10 + 0.03;

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
    
    case 2
    
    temp = A;
        
    % Get reference values
    temp_10 = zeros(1,size(temp,2));

    for j = 1:size(temp,2)
        temp_10(j) = interp1(z,temp(:,j),-10);
    end

    % mld as the shallowest depth where the temp. is 0.2 C colder than
    % a surface reference value at 10 m (see de Boyer Montégut et al. 2004)
    dt = temp_10 - 0.2;

    mld = zeros(size(dt'));

    where_mld = temp - repmat(dt,length(z),1);

    for j = 1:length(dt)

        % last negative value & first positive value
        mld_down = find(where_mld(:,j)<0,1,'last');
        mld_up = find(where_mld(:,j)>0,1,'first');

        % interpolate to find where the zero is
        mld(j) = -interp1(where_mld([mld_down,mld_up],j),...
            z([mld_down,mld_up]),0);
    end
end
        
end
