function grad_F = center_diff(F,coor,dim)

% center_diff
%==========================================================================
%
% USAGE:
%  grad_F = center_diff(F,coor,dim)
%
% DESCRIPTION:
%  Compute the first order derivative of a 2-D matrix along one dimension
%  using center difference method. The resulting gradient is evaluated at
%  the mid-point between adjacent grid points (which is input 'coor' here).
%
% INPUT:
%
%  F - the 2-D matrix whose gradient is being calculated
%  coor - an 1-D vector containing coordinates of F along one dimension 
%  dim - specify gradient in which dimension is being calculated (1,2)
%
% OUTPUT:
%
%  grad_F - gradient of F along one certain dimension
%
% AUTHOR:
%  September 19 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%

[m, n] = size(F); 

% determine differential dimension and repeating dimension
if dim == 1
    grad_F = zeros(m-1,n)*NaN;
    diff_n = m; 
    rep_n = n;
    
    for j=1:rep_n
        for i = 1:diff_n-1
            grad_F(i,j) = (F(i+1,j) - F(i,j))/(coor(i+1) - coor(i));
        end
    end

else
    grad_F = zeros(m,n-1)*NaN;
    diff_n = n;
    rep_n = m;
    
    for i=1:rep_n
        for j = 1:diff_n-1
            grad_F(i,j) = (F(i,j+1) - F(i,j))/(coor(j+1) - coor(j));
        end
    end
end


end

