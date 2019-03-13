function A_filled = vert_fill(A,z)
%
% vert_fill
%==========================================================================
%
% USAGE:
%  vert_fill(A,z)
%
% DESCRIPTION:
%  Function to vertically fill the missing data, using linear interpolation
%  as default method.
%
% INPUT:
%  
%  A - 2-D matrix with NaNs in columns
%  z - 1-D vector for the row coordinates of A 
%
% OUTPUT:
%
%  A_filled - matrix A with gaps filled
%
% AUTHOR:
%  March 12 2019. Zhihua Zheng                       [ zhihua@uw.edu ]

A_filled = A;

[~, cols] = find(isnan(A));
col_nan = unique(cols);

for j = 1:length(col_nan)
    
    tmp = A(:,col_nan(j));
    A_filled(:,col_nan(j)) = interp1(z(~isnan(tmp)),tmp(~isnan(tmp)),z);
end
    


end

