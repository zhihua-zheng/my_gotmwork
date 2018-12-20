%% do_sst_analysis

% Subroutine to ananlyze sst result from GOTM simulation output

% Zhihua Zheng, UW-APL, Oct. 18 2018

%% --------- read variables -----------------------------------------------
sst = out.sst;
sst_obs = out.sst_obs;
sst_from_prof = temp(128,:)';

%% --------- line plot with time ------------------------------------------ 

figure('position', [0, 0, 900, 300])
line(time,sst_from_prof,'LineWidth',.8,'Color',[.6 .4 .2])
line(time,sst_obs,'LineWidth',.8,'Color',[.3 .6 .4])

plot_info.grid_on = 1;
plot_info.y_lim = [0 25];
plot_info.xlabel = 'time';
plot_info.ylabel = 'sea surface temperature ($$^{\circ}C$$)';
plot_info.lgd = {'SMCLT','observation'};
%plot_info.save_path = [];
plot_info.save_path = './figs/sst_line';

line_annotate(plot_info)

% %%
% 
% line(time(1:10),sst_from_prof(1:10),'LineWidth',.8,'Color',[.6 .4 .2])
% line(time(1:10),sst_obs(1:10),'LineWidth',.8,'Color',[.3 .6 .4])
% 
% datetick('x','hh')

%% --------- scatter plot -------------------------------------------------

figure('position', [0, 0, 450, 450])
scatter(sst_obs,sst_from_prof,9,'filled')

axis equal
box on

tmp1 = 0.9*min(min([sst_obs sst_from_prof]));
tmp2 = 1.1*max(max([sst_obs sst_from_prof]));

xlim([tmp1 tmp2])
ylim([tmp1 tmp2])

h_ref = refline(1,0);
h_ref.Color = [.5 .5 .5];
h_ref.LineWidth = 1.5;

xlabel('obs. SST ($$^{\circ}C$$)', 'fontname',...
    'computer modern', 'fontsize', 14,'Interpreter', 'latex')
ylabel('Model SST ($$^{\circ}C$$)', 'fontname',...
    'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'fontsize',11,'fontname','computer modern',...
    'TickLabelInterpreter','latex')

export_fig('./figs/sst_scatter','-eps','-transparent','-painters')
