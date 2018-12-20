function scat_annotate(scat_info)
%
% scat_annotate
%==========================================================================
%
% USAGE:
%  scat_annotate(scat_info)
%
% DESCRIPTION:
%  Function to annotate a scatter plot
%
% INPUT:
%
%  scat_info - struct containing additional annotating information
%
% OUTPUT:
%
%  figure with annotation
%
% AUTHOR:
%  December 18 2018. Zhihua Zheng                       [ zhihua@uw.edu ]

box on
axis square

%% Grid option ------------------------------------------------------------
if scat_info.grid_on
    
    grid on
    set(gca,'gridlinestyle','--')
end

%% Label x axis -----------------------------------------------------------
if ~isempty(scat_info.xlabel) 
    xlabel(scat_info.xlabel, 'fontname',...
        'computer modern', 'fontsize', 18,'Interpreter', 'latex')    
end

%% Label y axis -----------------------------------------------------------
if ~isempty(scat_info.ylabel)
    ylabel(scat_info.ylabel, 'fontname',...
        'computer modern', 'fontsize', 18,'Interpreter', 'latex')
end

%% X-Y limit --------------------------------------------------------------
if ~isempty(scat_info.x_lim) 
    
    set(gca,'XLim',scat_info.x_lim,...
        'fontsize',18,'fontname','computer modern',...
        'XMinorTick','on','TickLabelInterpreter','latex')
else
    set(gca,'fontsize',18,'fontname','computer modern',...
        'XMinorTick','on','TickLabelInterpreter','latex')
end


if ~isempty(scat_info.y_lim)   
    set(gca,'YLim',scat_info.y_lim,'YMinorTick','on')
else
    set(gca,'YMinorTick','on')
end

%% X-Y ticks --------------------------------------------------------------
if ~isempty(scat_info.ticks)
    xticks(scat_info.ticks)
    yticks(scat_info.ticks)
end
%% Legend -----------------------------------------------------------------
if ~isempty(scat_info.lgd)
    
    if ~isempty(scat_info.lgd_pos)
        pos = scat_info.lgd_pos;
    else
        pos = 'best';
    end
    
    lgd = legend(scat_info.lgd_obj,scat_info.lgd,'Location',pos);
    set(lgd,'Interpreter','latex','fontsize',20,'color','none')
end

%% Save option ------------------------------------------------------------
if ~isempty(scat_info.save_path)
    
    export_fig(scat_info.save_path,'-eps','-transparent','-painters')
end

end