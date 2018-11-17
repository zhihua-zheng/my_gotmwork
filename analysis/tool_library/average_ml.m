function q_ml = average_ml(mld, q, zi, mld_smooth)

% average_ml
%==========================================================================
%
% USAGE:
%  q_ml = average_ml(mld, q, zi, mld_smooth)
%
% DESCRIPTION:
%  Compute the average of a quantity in the entire mixed layer 
%
% INPUT:
%
%  mld - 1-D vector in t dimension with values for mixed layer depth
%  q - A matrix (t,z) containing the quantity to be averaged
%  zi - 1-D vector in z dimension with values for interface depth in GOTM
%  mld_smooth - 1 or 0 (1 means mld is smoothed while 0 means original mld 
%    from GOTM output)
%
% OUTPUT:
%
%  q_ml - 1-D vector in t dimension. Mixed layer averages of quantity q
%
% AUTHOR:
%  October 23 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%

%% Identify where the mixed layer depth is

[MLD, Zi] = meshgrid(mld,zi);

where_mld = Zi + MLD;

% 1-D vector with depth index for mld 
mld_inx = ones(size(mld)); 

if mld_smooth 
    
    for j = 1:length(mld)
        % first positive value
        mld_inx(j) = find(where_mld(:,j)>0,1,'first'); 
    end
else
    % if mld isn't smoothed, we can use linear index to locate the zeros
    mld_inx = mld_inx*length(zi);
    [row, col] = ind2sub(size(where_mld),find(~where_mld));
    mld_inx(col) = row;    
end 

%% Average the qunatity in the mixed layer

q_ml = zeros(size(mld));

for j = 1:length(mld)
    
    h_ml = diff(zi(mld_inx(j):end)); % thickness of cells in mixed layer
    q_ml(j) = sum(q(mld_inx(j)-1:end,j).*h_ml)./(sum(h_ml));
end

end