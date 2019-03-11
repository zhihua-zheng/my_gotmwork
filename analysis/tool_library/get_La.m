function La = get_La(u_star,uStokes,vStokes,method,mld,z,h)
%
% get_La
%==========================================================================
%
% USAGE:
%  get_La(u_star,uStokes,vStokes,method,mld,z,h)
%
% DESCRIPTION:
%  Function to compute Langmuir number
%
% INPUT:
%
%  u_star - 1-D vector [t], friction velocity [m/s]
%  uStokes - 2-D matrix [t,z], x-direction Stokes drift velocity [m/s]
%  vStokes - 2-D matrix [t,z], y-direction Stokes drift velocity [m/s]
%  mld - 1-D vector [t], mixed layer depth [+,m]
%  z - 1-D vector [z], depth level for bins in GOTM [m]
%  h - 1-D vector [z], thickness for bins in GOTM [m]
%  method - type of Langmuir number being computed
%    1 - turbulent Langmuir number (McWilliams et al. 1997)
%    2 - surface layer Langmuir number (Harcourt & D'Asaro 2008)
%
% OUTPUT:
%
%  La - nondimensional Langmuir number
%
% AUTHOR:
%  December 19 2018. Zhihua Zheng                       [ zhihua@uw.edu ]


uSt = sqrt(uStokes.^2 + vStokes.^2);

switch method
    case 1
        La = sqrt(u_star'./uSt(end,:));
        
    case 2
        z_sl = mld/5; % surface layer depth
        z_ref = 0.765*mld; % reference level
        
        % find index for surface layer depth in z
        [SL, Z] = meshgrid(z_sl,z);
        where_sl = Z + SL;
        sl_inx = ones(size(z_sl))*length(z); 
        for j = 1:length(z_sl)
            if where_sl(end,j)>=0
                % first positive value
                sl_inx(j) = find(where_sl(:,j)>=0,1,'first'); 
            end
        end
        
        % average the Stokes drift within surface layer
        uSt_sl = zeros(size(z_sl));
        for j = 1:length(z_sl)
            h_sl = h(sl_inx(j):end);
            uSt_sl(j) = nansum(uSt(sl_inx(j):end,j).*h_sl)./(sum(h_sl));
        end
        
        % find index for reference level in z
        [REF, Z] = meshgrid(z_ref,z);
        where_ref = Z + REF;
        ref_inx = ones(size(z_ref))*length(z); 
        for j = 1:length(z_ref)
            % first positive value
            if where_ref(end,j)>=0
                ref_inx(j) = find(where_ref(:,j)>=0,1,'first'); 
            end
        end
        
        uSt_ref = uSt(sub2ind(size(uSt),ref_inx,(1:length(z_ref))'));
        
        La = sqrt(u_star./(uSt_sl - uSt_ref));    
end

        
end