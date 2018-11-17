

%% do_mld_scatter

figure('position', [0, 0, 690, 680])

%-------- for turb_method{2} --------------------------------------------------

mld_model(:,2) = vertcat(all_mld{1:6,2});
mld_obs(:,2) = vertcat(all_mld_obs{1:6,2});

dummy = -rand(size(mld_model(:,2))); % dummy variable for 3D plot
s2 = scatter3(mld_obs(:,2),mld_model(:,2),dummy,40,'d','filled',...
    'MarkerFaceColor',rgb('bright red'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
s2.MarkerFaceAlpha = 0.6;
hold on

sl2 = mld_obs(:,2)\mld_model(:,2);

%-------- for turb_method{1} --------------------------------------------------

mld_model(:,1) = vertcat(all_mld{1:6,1});
mld_obs(:,1) = vertcat(all_mld_obs{1:6,1});

dummy = -rand(size(mld_model(:,1))); % dummy variable for 3D plot
s1 = scatter3(mld_obs(:,1),mld_model(:,1),dummy,40,'d','filled',...
    'MarkerFaceColor',rgb('neon blue'),'LineWidth',1,...
    'MarkerEdgeColor',[.65 .65 .65]);
s1.MarkerFaceAlpha = 0.6;

sl1 = mld_obs(:,1)\mld_model(:,1);

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
h_ref = plot3((0:5:160),(0:5:160),zeros(size((0:5:160))),'Color',[.2 .2 .2],'LineWidth',3,'LineStyle','--');
hold on
fit2 = plot3((0:5:160),sl2*(0:5:160),zeros(size((0:5:160))),'Color',rgb('light red'),'LineWidth',4);
hold on
fit1 = plot3((0:5:160),sl1*(0:5:160),zeros(size((0:5:160))),'Color',rgb('bright sky blue'),'LineWidth',4);
hold off

lgd = legend([s1 s2 fit1 fit2 h_ref],[turb_method{1},', rmse $\sim$ ',...
    num2str(round(all_mld_rmse(1),2))],[turb_method{2},', rmse $\sim$ ',...
    num2str(round(all_mld_rmse(2),2))],['slope = ',num2str(round(sl1,2))],...
    ['slope = ',num2str(round(sl2,2))],'slope = 1','Location','best');
set(lgd,'Interpreter','latex','fontsize', 20)

xlabel('\textit{obs. MLD} ($$m$$)', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
ylabel('\textit{model MLD} ($$m$$)', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
set(gca,'fontsize',20,'fontname','computer modern','gridlinestyle','--',...
    'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')

set(gca,'LooseInset', get(gca,'TightInset')); % no blank edge
saveas(gcf, './figs/all_mld_scatter', 'png');

