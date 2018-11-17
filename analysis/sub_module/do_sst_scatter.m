

%% do_sst_scatter

figure('position', [0, 0, 690, 680])

%-------- for turb_method{2} ----------------------------------------------

sst_model(:,2) = vertcat(all_sst{:,2});
sst_obs(:,2) = vertcat(all_sst_obs{:,2});

dummy = -rand(size(sst_model(:,2))); % dummy variable for 3D plot
s2 = scatter3(sst_obs(:,2),sst_model(:,2),dummy,35,'filled',...
    'MarkerFaceColor',rgb('bright red'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
s2.MarkerFaceAlpha = 0.6;
hold on

sl2 = sst_obs(:,2)\sst_model(:,2);

%-------- for turb_method{1} ----------------------------------------------

sst_model(:,1) = vertcat(all_sst{:,1});
sst_obs(:,1) = vertcat(all_sst_obs{:,1});

dummy = -rand(size(sst_model(:,1))); % dummy variable for 3D plot
s1 = scatter3(sst_obs(:,1),sst_model(:,1),dummy,35,'filled',...
    'MarkerFaceColor',rgb('neon blue'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
s1.MarkerFaceAlpha = 0.6;

sl1 = sst_obs(:,1)\sst_model(:,1);

% Change viewpoint 
view(2)

hold off 
axis square
box on

% xlim([0.8*min(min([sst_min sst_obs_min])) 1.3*max(max([sst_max sst_obs_max]))])
% ylim([0.8*min(min([sst_min sst_obs_min])) 1.3*max(max([sst_max sst_obs_max]))])
xlim([4 21])
ylim([4 21])
xticks((4:2:20))
yticks((4:2:20))

hold on
h_ref = plot3((4:1:21),(4:1:21),zeros(size((4:1:21))),'Color',[.2 .2 .2],'LineWidth',3,'LineStyle','--');
hold on
fit2 = plot3((4:1:21),sl2*(4:1:21),zeros(size((4:1:21))),'Color',rgb('light red'),'LineWidth',4);
hold on
fit1 = plot3((4:1:21),sl1*(4:1:21),zeros(size((4:1:21))),'Color',rgb('bright sky blue'),'LineWidth',4);
hold off

lgd = legend([s1 s2 fit1 fit2 h_ref],[turb_method{1},', rmse $\sim$ ',...
    num2str(round(all_sst_rmse(1),2))],[turb_method{2},', rmse $\sim$ ',...
    num2str(round(all_sst_rmse(2),2))],['slope = ',num2str(round(sl1,2))],...
    ['slope = ',num2str(round(sl2,2))],'slope = 1','Location','best');
set(lgd,'Interpreter','latex','fontsize', 25)

xlabel('\textit{obs. SST} ($$^{\circ}C$$)', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
ylabel('\textit{model SST} ($$^{\circ}C$$)', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
set(gca,'fontsize',20,'fontname','computer modern','gridlinestyle','--',...
    'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')

set(gca,'LooseInset', get(gca,'TightInset')); % no blank edge
saveas(gcf, './figs/all_sst_scatter', 'png');

