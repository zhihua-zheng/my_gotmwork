
%% do_sst_scatter_season

% TO-DO: 
%        3.look the time series of deviation
%        4.do a spectral analysis of deviation

season_rmse = zeros(4,2); % rmse for 4 seasons

figure('position', [0, 0, 690, 680])
   
for j = 1:4
        
    subplot(2,2,j)
    
    % index for season in cell array, update for every j
    all_season_inx = cellfun(@(A) ismember(A,season_list(j)),...
        all_month,'UniformOutput',false);
    
    % number of points in each seasons
    season_num = cellfun(@(A,c) numel(A(c)),all_sst,all_season_inx,'UniformOutput',false);
    
    % root mean square for SST
    rmse = cellfun(@(A,B,c) sqrt(sum(A(c))/B),sst_se,season_num,all_season_inx);
    
    % SST in one season
    sst_season = cellfun(@(A,c) A(c),all_sst,all_season_inx,...
        'UniformOutput',false);
    sst_obs_season = cellfun(@(A,c) A(c),all_sst_obs,all_season_inx,...
        'UniformOutput',false);
    
    % composite averaged rmse
    season_rmse(j,:) = mean(rmse);
    
    % ---- for turb_method{2}
    y(:,2) = vertcat(sst_season{:,2});
    x(:,2) = vertcat(sst_obs_season{:,2});
 
    dummy = -rand(size(y(:,2)));
    s2 = scatter3(x(:,2),...
        y(:,2),dummy,20,'filled',...
    'MarkerFaceColor',rgb('bright red'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
    s2.MarkerFaceAlpha = 0.6;
    hold on
    
    sl2 = x(:,2)\y(:,2);

    % ---- for turb_method{1}
    
    % index for season in matrix
    y(:,1) = vertcat(sst_season{:,1});
    x(:,1) = vertcat(sst_obs_season{:,1});

    dummy = -rand(size(y(:,1)));
    s2 = scatter3(x(:,1),...
        y(:,1),dummy,20,'filled',...
    'MarkerFaceColor',rgb('neon blue'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
    s2.MarkerFaceAlpha = 0.6;
    
    sl1 = x(:,1)\y(:,1);

    % Change viewpoint 
    view(2)
    
    hold off
    axis square
    box on
    
%     xlim([0.9*min(min(sst_min)) 1.01*max(max(sst_max))])
%     ylim([0.9*min(min(sst_min)) 1.01*max(max(sst_max))])
    
    xlim([4 21])
    ylim([4 21])
    xticks((4:2:20))
    yticks((4:2:20))
    
    hold on
    h_ref = plot3((4:1:21),(4:1:21),zeros(size((4:1:21))),'Color',[.2 .2 .2],'LineWidth',3,'LineStyle','--');
    hold on
    fit2 = plot3((4:1:21),sl2*(4:1:21),zeros(size((4:1:21))),'Color',rgb('light red'),'LineWidth',3);
    hold on
    fit1 = plot3((4:1:21),sl1*(4:1:21),zeros(size((4:1:21))),'Color',rgb('bright sky blue'),'LineWidth',3);
    hold off
    
    lgd = legend([fit1 fit2 h_ref],['slope = ',...
        num2str(round(sl1,2))],['slope = ',...
        num2str(round(sl2,2))],'slope = 1','Location','southeast');
    set(lgd,'Interpreter','latex','fontsize',16,'color','none')
    
    title([season_str{j}],'fontsize',20,'fontname',...
        'computer modern','Interpreter','latex')
    set(gca,'fontsize',15,'fontname','computer modern','gridlinestyle','--',...
    'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')
    
    clear x y
end

[ax1, h1] = suplabel('\textit{model SST} ($^{\circ}C$)','y');
[ax2, h2] = suplabel('\textit{obs. SST} ($^{\circ}C$)','x');
set(h1,'FontSize',20,'FontName','computer modern','Interpreter','latex');
set(h2,'FontSize',20,'FontName','computer modern','Interpreter','latex');

% lgd = legend([s1 s2],closure);
% newPosition = [0.38 0.49 0.3 0.05];
% newUnits = 'normalized';
% set(lgd,'Position',newPosition,'Units',newUnits,'Interpreter',...
%     'latex','fontsize',15,'orientation','horizontal');

set(gca,'LooseInset', get(gca,'TightInset')); % no blank edge
saveas(gcf,'./figs/season_sst_scatter', 'png');

% ------ sesonal RMSE line plot -------------------------------------------
figure('position', [0, 0, 390, 680])
line(season_rmse(:,1),(1:1:4),'Marker','s','MarkerSize',15,'MarkerFaceColor',...
    rgb('neon blue'),'MarkerEdgeColor',[.5 .5 .5],'Color',...
    rgb('bright sky blue'),'LineWidth',2)
line(season_rmse(:,2),(1:1:4),'Marker','s','MarkerSize',15,'MarkerFaceColor',...
    rgb('bright red'),'MarkerEdgeColor',[.5 .5 .5],'Color',...
    rgb('light red'),'LineWidth',2)
line([1 1],[0.5 4.5],'Color',[.3 .3 .3],'LineStyle','--','LineWidth',2)

ytickangle(45)
set(gca,'ytick',(1:1:4),'yticklabel',season_str)
set(gca,'XAxisLocation','top','ydir','reverse');

spec_info.grid_on = 1;
spec_info.lgd = turb_method;
spec_info.lgd_pos = 'east';
spec_info.x_lim = [0.9*min(min(season_rmse)) 1.05*max(max(season_rmse))];
spec_info.y_lim = [0.5 4.5];
spec_info.xlabel = 'SST rmse $(^{\circ}C)$';
spec_info.ylabel = [];
spec_info.save_path = './figs/season_rmse_sst';
%spec_info.save_path = [];

line_annotate(spec_info)
% -------------------------------------------------------------------------