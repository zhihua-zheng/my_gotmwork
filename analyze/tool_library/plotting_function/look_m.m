function look_m(mat,plot_method)

% test_m
%==========================================================================
%
% USAGE:
%  test_m(mat)

%
% DESCRIPTION:
%  Give a general look into the structure of matrix
%
% INPUT:
%
%  mat - input 2-D matrix
%  plot_method - scalar to specify which plotting method to use 
%    contourf - 1
%      pcolor - 2
%     imagesc - 3
%
% OUTPUT:
%
%  figure displaying the structure of the matrix
%
% AUTHOR:
%  September 17 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%

figure('position', [0, 0, 900, 300])

switch plot_method
    case 1 % contourf
        contourf(mat,'LineStyle','none')
    
    case 2 % pcolor
        [nr,nc] = size(mat);
        h = pcolor([mat nan(nr,1); nan(1,nc+1)]); % padding the last row and column 
        set(h, 'EdgeColor', 'none');
        
    case 3 % imagesc
        imAlpha=ones(size(mat));
        imAlpha(isnan(mat))=0; % set color of NaNs to background color
        imagesc(mat,'AlphaData',imAlpha)
        axis('xy') % flip y axis, imagesc mapps the matrix from lower-left corner
end

cmocean('balance')
c_max = max(max(abs(mat)));
caxis([-c_max c_max]);

colorbar

end

