function   plot_time_depth(t, d, scalar, plot_info)

% plot_time_depth
%==========================================================================
%
% USAGE:
%  plot_time_depth(t, d, scalar, plot_info)
%
% DESCRIPTION:
%  Function to plot temporal evolution of scalar profile
%
% INPUT:
%
%  t - one dimensional array of date number
%  d - one dimensional array of depth (deep to shallow)
%  scalar - two dimensional array of scalar value
%  plot_info - struct containing additional annotating information
%
% OUTPUT:
%
%  figure
%
% AUTHOR:
%  September 2 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%

%% determine color value limits -------------------------------------------
if isempty(plot_info.clim)
    
    CL = [min(min(scalar)) max(max(scalar))];  % default option
elseif strcmp(plot_info.clim,'symmetric')
    
    tmp1 = min(min(scalar));
    tmp2 = max(max(scalar));
    tmp = max(abs(tmp1), abs(tmp2)); % the one with largest magnitude
    CL = [-tmp tmp];
else
    CL = plot_info.clim;
end

%% Plot -------------------------------------------------------------------
cnum = 15;
conts = linspace(floor(CL(1)),ceil(CL(2)),cnum);
[T, Z] = meshgrid(t,d);

figure('position', [0, 0, 900, 300])

switch plot_info.plot_method
    case 1
        contourf(T,Z,scalar,conts,'LineWidth',0.01,'LineStyle','none')
    case 2
        [nr,nc] = size(scalar);
        % padding the last row and column, since pcolor ignores them 
        h = pcolor([T nan(nr,1); nan(1,nc+1)],...
            [Z nan(nr,1); nan(1,nc+1)],[scalar nan(nr,1); nan(1,nc+1)]); 
        set(h, 'EdgeColor', 'none');
    case 3
        imAlpha=ones(size(scalar));
        imAlpha(isnan(scalar))=0; % set color of NaNs to background color
        imagesc(t,d,scalar,'AlphaData',imAlpha)
        % flip y axis, imagesc mapps the matrix from lower-left corner
        axis('xy') 
end

cmocean(plot_info.color)

%% Annotation -------------------------------------------------------------
  caxis(CL);
  box on
  datetick('x',plot_info.timeformat)
  ylabel(plot_info.ylabel, 'fontname', 'computer modern', 'fontsize', ...
      18,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[t(1) t(end)],'YLim',plot_info.ylim,'fontsize',...
      18,'fontname','computer modern','TickLabelInterpreter','latex')
  
  h = colorbar('EastOutside');
  h.Label.String = plot_info.clabel;
  h.Label.Interpreter = 'latex';
  h.Label.FontName = 'computer modern';
  h.Label.FontSize = 20;
  set(h,'TickLabelInterpreter','latex','fontsize',16);

%% save or not ------------------------------------------------------------
if ~isempty(plot_info.save_path)
    
  set(gca,'LooseInset', get(gca,'TightInset')); % no blank edge
  saveas(gcf, plot_info.save_path, 'epsc');
  
  % clean the white lines in the patch
  epsclean([plot_info.save_path,'.eps'],'closeGaps',true) 
end

end