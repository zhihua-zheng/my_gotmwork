
%% do_mld_scatter_season

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
    season_num = cellfun(@(A,c) numel(A(c)),all_mld,all_season_inx,'UniformOutput',false);
    
    % root mean square for MLD
    rmse = cellfun(@(A,B,c) sqrt(sum(A(c))/B),mld_se,season_num,all_season_inx);
    
    % mld for season
    mld_season = cellfun(@(A,c) A(c),all_mld,all_season_inx,...
        'UniformOutput',false);
    mld_obs_season = cellfun(@(A,c) A(c),all_mld_obs,all_season_inx,...
        'UniformOutput',false);
    
    % composite averaged rmse
    season_rmse(j,:) = mean(rmse(1:6,:));
    
    % ---- for turb_method{2}
    y(:,2) = vertcat(mld_season{1:6,2});
    x(:,2) = vertcat(mld_obs_season{1:6,2});
 
    dummy = -rand(size(y(:,2)));
    s2 = scatter3(x(:,2),...
        y(:,2),dummy,20,'d','filled',...
    'MarkerFaceColor',rgb('bright red'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
    s2.MarkerFaceAlpha = 0.6;
    hold on
    
    sl2 = x(:,2)\y(:,2);

    % ---- for turb_method{1}
    
    % index for season in matrix
    y(:,1) = vertcat(mld_season{1:6,1});
    x(:,1) = vertcat(mld_obs_season{1:6,1});

    dummy = -rand(size(y(:,1)));
    s2 = scatter3(x(:,1),...
        y(:,1),dummy,20,'d','filled',...
    'MarkerFaceColor',rgb('neon blue'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
    s2.MarkerFaceAlpha = 0.6;
    
    sl1 = x(:,1)\y(:,1);

    % Change viewpoint 
    view(2)
    
    hold off
    axis square
    box on
    
    % xlim([0.8*min(min([mld_min mld_obs_min])) 1.1*max(max([mld_max mld_obs_max]))])
    % ylim([0.8*min(min([mld_min mld_obs_min])) 1.1*max(max([mld_max mld_obs_max]))])
    xlim([0 160])
    ylim([0 160])
    xticks((0:20:160))
    yticks((0:20:160))

    hold on
    h_ref = plot3((0:5:160),(0:5:160),zeros(size((0:5:160))),'Color',[.2 .2 .2],'LineWidth',2,'LineStyle','--');
    hold on 
    fit2 = plot3((0:5:160),sl2*(0:5:160),zeros(size((0:5:160))),'Color',rgb('light red'),'LineWidth',3);
    hold on
    fit1 = plot3((0:5:160),sl1*(0:5:160),zeros(size((0:5:160))),'Color',rgb('bright sky blue'),'LineWidth',3);
    
    lgd = legend([fit1 fit2 h_ref],['slope = ',...
        num2str(round(sl1,2))],['slope = ',...
        num2str(round(sl2,2))],'slope = 1','Location','northwest');
    set(lgd,'Interpreter','latex','fontsize',15.5,'color','none')
    
    title([season_str{j}],'fontsize',20,'fontname',...
        'computer modern','Interpreter','latex')
    set(gca,'fontsize',15,'fontname','computer modern','gridlinestyle','--',...
    'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')
    
    clear x y tmp
end

[ax1, h1] = suplabel('\textit{model MLD} ($m$)','y');
[ax2, h2] = suplabel('\textit{obs. MLD} ($m$)','x');
set(h1,'FontSize',20,'FontName','computer modern','Interpreter','latex');
set(h2,'FontSize',20,'FontName','computer modern','Interpreter','latex');

% lgd = legend([s1 s2],closure);
% newPosition = [0.38 0.49 0.3 0.05];
% newUnits = 'normalized';
% set(lgd,'Position',newPosition,'Units',newUnits,'Interpreter',...
%     'latex','fontsize',15,'orientation','horizontal');

set(gca,'LooseInset', get(gca,'TightInset')); % no blank edge
saveas(gcf,'./figs/season_mld_scatter', 'png');

% ------ sesonal RMSE line plot -------------------------------------------
figure('position', [0, 0, 390, 680])
line(season_rmse(:,1),(1:1:4),'Marker','s','MarkerSize',15,'MarkerFaceColor',...
    rgb('neon blue'),'MarkerEdgeColor',[.5 .5 .5],'Color',...
    rgb('bright sky blue'),'LineWidth',2)
line(season_rmse(:,2),(1:1:4),'Marker','s','MarkerSize',15,'MarkerFaceColor',...
    rgb('bright red'),'MarkerEdgeColor',[.5 .5 .5],'Color',...
    rgb('light red'),'LineWidth',2)
line([10 10],[0.5 4.5],'Color',[.3 .3 .3],'LineStyle','--','LineWidth',2)

ytickangle(45)
set(gca,'ytick',(1:1:4),'yticklabel',season_str)
set(gca,'XAxisLocation','top','ydir','reverse');

spec_info.grid_on = 1;
spec_info.lgd = turb_method;
spec_info.lgd_pos = [];
spec_info.x_lim = [0.9*min(min(season_rmse)) 1.05*max(max(season_rmse))];
spec_info.y_lim = [0.5 4.5];
spec_info.xlabel = 'MLD rmse $(m)$';
spec_info.ylabel = [];
spec_info.save_path = './figs/season_rmse_mld';
%spec_info.save_path = [];

line_annotate(spec_info)
% -------------------------------------------------------------------------