function   line_annotate(plot_info)

% line_annotate
%==========================================================================
%
% USAGE:
%  line_annotate(plot_info)
%
% DESCRIPTION:
%  Function to annotate a line plot
%
% INPUT:
%
%  plot_info - struct containing additional annotating information
%
% OUTPUT:
%
%  figure with annotation
%
% AUTHOR:
%  September 4 2018. Zhihua Zheng                       [ zhihua@uw.edu ]

box on

%% Grid option ------------------------------------------------------------
if plot_info.grid_on
    
    grid on
    set(gca,'gridlinestyle','--')
end

%% Label x axis -----------------------------------------------------------
if ~isempty(plot_info.xlabel)
    
    if strcmp(plot_info.xlabel, 'time')
        datetick('x',plot_info.timeformat)
    else 
        xlabel(plot_info.xlabel, 'fontname',...
        'computer modern', 'fontsize', 18,'Interpreter', 'latex')    
    end
end

%% Label y axis -----------------------------------------------------------
if ~isempty(plot_info.ylabel)
    ylabel(plot_info.ylabel, 'fontname',...
        'computer modern', 'fontsize', 18,'Interpreter', 'latex')
end

%% X-Y limit --------------------------------------------------------------
if ~isempty(plot_info.x_lim) 
    
    % use setDateAxes() instead of set() to be compatible with time axes
    setDateAxes(gca,'XLim',plot_info.x_lim,...
        'fontsize',18,'fontname','computer modern',...
        'XMinorTick','on','TickLabelInterpreter','latex')
else
    setDateAxes(gca,'fontsize',18,'fontname','computer modern',...
        'XMinorTick','on','TickLabelInterpreter','latex')
end


if ~isempty(plot_info.y_lim)   
    set(gca,'YLim',plot_info.y_lim,'YMinorTick','on')
else
    set(gca,'YMinorTick','on')
end
    
%% Legend -----------------------------------------------------------------
if ~isempty(plot_info.lgd)
    
    if ~isempty(plot_info.lgd_pos)
        pos = plot_info.lgd_pos;
    else
        pos = 'best';
    end
    
    lgd = legend(plot_info.lgd,'Location',pos);
    set(lgd,'Interpreter','latex','fontsize',20,'color','none')
end

%% Save option ------------------------------------------------------------
if ~isempty(plot_info.save_path)
    
    export_fig(plot_info.save_path,'-png','-transparent','-painters')
end

end