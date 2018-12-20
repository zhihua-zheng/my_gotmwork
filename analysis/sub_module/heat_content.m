%% heat_content

% Subroutine to ananlyze the heat content of the water column simualted

% Zhihua Zheng, UW-APL, Sep. 5 2018

%% --------- read relevant variables --------------------------------------

int_total = out.int_total;
int_heat = out.int_heat;
int_swr = out.int_swr;
temp = out.temp;
temp_obs = out.temp_obs;
salt = out.salt;
salt_obs = out.salt_obs;

% decide to use yearday or time string for x lable
if lable_as_yd 
    plot_info.xlabel = []; 
    t = yd;
else
    plot_info.xlabel = 'time';
    t = time;
end
    
%% Plot integrated heat info from observation ----------------------------

figure('position', [0, 0, 900, 300])

line(t,int_total./10^(6),'LineWidth',3,'Color',[.2 .6 .7])
line(t,int_heat./10^(6),'LineWidth',.8,'Color',[.8 .2 .2])
line(t,int_swr./10^(6),'LineWidth',.8,'Color',[.3 .2 .5])
line(t,zeros(size(t)),'LineWidth',.6,'Color',[.3 .3 .3],'LineStyle','--')

% mark the net heat taken by the water column from surface measurement
line([t(end) t(end)],[0 int_total(end)/10^(6)],'Color',...
    [.2 .2 .2],'LineStyle','-','LineWidth',4);
% text(t(end-8500),-int_total(end)/10^6,...
%     ['net heat into the ocean $\sim$ ',num2str(round(int_total(end)/10^(6))),...
%     ' $MJ/m^{2}$'],'fontname','computer modern','Interpreter','latex','fontsize',15)

% figure specification
plot_info.grid_on = 1;
plot_info.ylabel = 'integrated heat ($MJ/m^{2}$)';
plot_info.x_lim = [t(1) t(end)];
plot_info.y_lim = [];
plot_info.lgd = {'total heat exchange','surface heat flux',...
    'short wave radiation'};
plot_info.lgd_pos = [];
%plot_info.save_path = './figs/int_heat';
plot_info.save_path = [];

line_annotate(plot_info)
%--------------------------------------------------------------------------


%% Calculation of heat content and rate of change with time ---------------

% HC for the whole domain

cp = 3985; % specific heat capacity from GOTM setup

rho = sw_dens0(salt,temp);
HC = sum(rho.*(temp+273.15).*h*cp); %[J/m^2]
HC_delta = HC(end) - HC(1);
HC_t = gradient(HC,dt); %[W/m^2]

% from observation
rho_obs = sw_dens0(salt_obs,temp_obs);
HC_obs = sum(rho_obs.*(temp_obs+273.15).*h*cp);
HC_delta_obs = HC_obs(end) - HC_obs(1);
HC_t_obs = gradient(HC_obs,dt);

% from surface heat info
HC_surface = int_total+HC(1);
HC_t_surface = gradient(int_total,dt);

HC_error = 100*(int_total(end) - HC_delta)/int_total(end); % in percent
%--------------------------------------------------------------------------


%% Heat content time series -----------------------------------------------

figure('position', [0, 0, 900, 300])

line(t,HC./10^(6),'LineWidth',.3,'Color',[.7 .4 .6])
line(t,HC_surface./10^(6),'LineWidth',.3,'Color',[.1 .6 .7])
line(t,HC_obs./10^(6),'LineWidth',.05,'Color',[.4 .8 .6])
line(t,ones(size(t))*HC(1)/10^(6),'LineWidth',.9,'Color',[.3 .3 .3],'LineStyle','--')

% mark the deviation from surface heat input
line([t(end) t(end)],[HC(1)/10^(6) HC_surface(end)./10^(6)],...
    'Color',[.6 .1 .3],'LineStyle','-','LineWidth',3);
% text(t(end-2200),HC(end)/10^(6),...
%     ['$\sim$ ',num2str(round(HC_error,2)),'$\%$'],'color',[.6 .1 .3],...
%     'fontname','computer modern','Interpreter','latex','fontsize',13)

% mark the heat content change relative to observation
line([t(end) t(end)],[HC(1)/10^(6) HC(end)/10^(6)],'Color',...
    [.1 .1 .1],'LineStyle','-','LineWidth',3);
% text(t(end-8500),HC(8800)/10^(6),[turb_method,...
%     ' heat content change $\sim$ ',num2str(round((HC_delta-HC_delta_obs)/10^(6))),...
%     ' $MJ/m^{2}$'],'Color',[.1 .1 .1],'fontname','computer modern',...
%     'Interpreter','latex','fontsize',13)
  
% figure specification
plot_info.lgd = {[turb_method,' HC'],'surface heat exchange','obs. HC'};
plot_info.ylabel = 'heat content ($MJ/m^{2}$)';
% plot_info.save_path = './figs/HC';
plot_info.save_path = [];

line_annotate(plot_info)  
%--------------------------------------------------------------------------

%% Plot temporal variation of heat content ---------------------------------

figure('position', [0, 0, 900, 300])

line(t,HC_t_surface,'LineWidth',1,'Color','k')
line(t,HC_t,'LineWidth',1,'Color',[.4 .9 .7])
line(t,zeros(size(t)),'LineWidth',.6,'Color',[.5 .5 .5],'LineStyle',':')
 
% figure specification
plot_info.lgd = {'surface heat exchange rate',...
    ['$\partial_{t}HC$ in ',turb_method]};
plot_info.ylabel = 'temporal heat variation ($W/m^{2}$)';
% plot_info.save_path = './figs/HC_t_surf';
plot_info.save_path = [];

line_annotate(plot_info)

%--------
figure('position', [0, 0, 900, 300])

line(t,HC_t_obs,'LineWidth',1,'Color','k')
line(t,HC_t,'LineWidth',1,'Color',[.4 .9 .7])
line(t,zeros(size(t)),'LineWidth',.6,'Color',[.5 .5 .5],'LineStyle',':')

plot_info.lgd = {'$\partial_{t}HC$ in obs.', ...
    ['$\partial_{t}HC$ in ',turb_method]};
plot_info.ylabel = 'temporal heat variation ($W/m^{2}$)';
% plot_info.save_path = './figs/HC_t_obs';
plot_info.save_path = [];

line_annotate(plot_info)
