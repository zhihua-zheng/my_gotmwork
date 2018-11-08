function   line_annotate(spec_info)

% line_annotate
%==========================================================================
%
% USAGE:
%  line_annotate(spec_info)
%
% DESCRIPTION:
%  Function to annotate a line plot
%
% INPUT:
%
%  spec_info - struct containing additional annotating information
%
% OUTPUT:
%
%  figure with annotation
%
% AUTHOR:
%  September 4 2018. Zhihua Zheng                       [ zhihua@uw.edu ]

box on

%% Grid option ------------------------------------------------------------
if spec_info.grid_on
    
    grid on
    set(gca,'gridlinestyle','--')
end

%% Label x axis -----------------------------------------------------------
if ~isempty(spec_info.xlabel)
    
    if strcmp(spec_info.xlabel, 'time')
        datetick('x',spec_info.timeformat)
    else 
        xlabel(spec_info.xlabel, 'fontname',...
        'computer modern', 'fontsize', 18,'Interpreter', 'latex')    
    end
end

%% Label y axis -----------------------------------------------------------
if ~isempty(spec_info.ylabel)
    ylabel(spec_info.ylabel, 'fontname',...
        'computer modern', 'fontsize', 18,'Interpreter', 'latex')
end

%% X-Y limit --------------------------------------------------------------
if ~isempty(spec_info.x_lim) 
    
    % use setDateAxes() instead of set() to be compatible with time axes
    setDateAxes(gca,'XLim',spec_info.x_lim,...
        'fontsize',18,'fontname','computer modern',...
        'XMinorTick','on','TickLabelInterpreter','latex')
else
    setDateAxes(gca,'fontsize',18,'fontname','computer modern',...
        'XMinorTick','on','TickLabelInterpreter','latex')
end


if ~isempty(spec_info.y_lim)   
    set(gca,'YLim',spec_info.y_lim,'YMinorTick','on')
else
    set(gca,'YMinorTick','on')
end
    
%% Legend -----------------------------------------------------------------
if ~isempty(spec_info.lgd)
    
    if ~isempty(spec_info.lgd_pos)
        pos = spec_info.lgd_pos;
    else
        pos = 'best';
    end
    
    lgd = legend(spec_info.lgd,'Location',pos);
    set(lgd,'Interpreter','latex','fontsize', 20)
end

%% Save option ------------------------------------------------------------
if ~isempty(spec_info.save_path)
    
    export_fig(spec_info.save_path,'-eps','-transparent','-painters')
end

end