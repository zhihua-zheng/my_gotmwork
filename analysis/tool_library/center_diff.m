function grad_F = center_diff(F,coor,dim)
%
% center_diff
%==========================================================================
%
% USAGE:
%  grad_F = center_diff(F,coor,dim)
%
% DESCRIPTION:
%  Compute the first order derivative of a 3-D matrix along one dimension
%  using center difference method. The resulting gradient is evaluated at
%  the mid-point between adjacent grid points (which is input 'coor' here).
%
% INPUT:
%
%  F - the matrix whose gradient is to be calculated
%  coor - an 1-D vector containing coordinates of F along one dimension 
%  dim - specify gradient in which dimension is being calculated (1,2,3)
%
% OUTPUT:
%
%  grad_F - gradient of F along one certain dimension
%
% AUTHOR:
%  September 19 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%  May       22 2019. Expand to 3-D 
%==========================================================================

[m, n, l] = size(F); 

% determine differential dimension and repeating dimension

switch dim 
    case 1 % first dimension
    grad_F = zeros(m-1,n,l)*NaN;
    
    for k = 1:l
        for j = 1:n
            for i = 1:m-1
                grad_F(i,j,k) = (F(i+1,j,k) - F(i,j,k))/(coor(i+1) - coor(i));
            end
        end
    end

    case 2 % second dimension
    grad_F = zeros(m,n-1,l)*NaN;
    
    for k = 1:l
        for i = 1:m
            for j = 1:n-1
                grad_F(i,j,k) = (F(i,j+1,k) - F(i,j,k))/(coor(j+1) - coor(j));
            end
        end
    end
    
    case 3 % third dimension
    grad_F = zeros(m,n-1,l)*NaN;
    
    for j = 1:n
        for i = 1:m
            for k = 1:l-1
                grad_F(i,j,k) = (F(i,j,k+1) - F(i,j,k))/(coor(k+1) - coor(k));
            end
        end
    end 
    
end


end

