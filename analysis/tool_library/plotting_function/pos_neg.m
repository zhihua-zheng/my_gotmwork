function pos_neg(x,A,Ac)
%
% pos_neg
%==========================================================================
%
% USAGE:
%  pos_neg(x,A,Ac)
%
% DESCRIPTION:
%  Function to generate a line plot with different shading aread 
%  above and below a critical value Ac
%
% INPUT:
%
%  x - 1-D vector of the horizontal coordinates
%  A - 1-D vector of the quantity to displayed in the plot
%  Ac - critical value of A
%
% OUTPUT:
%
%  figure of the line plot with color shading above and below Ac
%
% AUTHOR:
%  July 19 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

pos = A;
neg = A;

pos(pos-Ac<0) = nan;
neg(neg-Ac>0) = nan;

area(x,neg,'linewidth',.8,'FaceColor',rgb('cerulean'),'FaceAlpha',.6,'EdgeAlpha',.5)
hold on
area(x,pos,'linewidth',.8,'FaceColor',rgb('coral')   ,'FaceAlpha',.6,'EdgeAlpha',.5)
xlim(x([1,end]))

end

