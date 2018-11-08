%% heat_content

% Subroutine to ananlyze the heat content of the water column simualted

% Zhihua Zheng, UW-APL, Sep. 5 2018

%% --------- read relevant variables --------------------------------------

int_total = out.int_total;
int_heat = out.int_heat;
int_swr = out.int_swr;
rho = out.rho;
temp = out.temp;
temp_obs = out.temp_obs;

%% Plot integrated heat infor from observation ----------------------------

figure('position', [0, 0, 900, 300])

line(time,int_total./10^(6),'LineWidth',3,'Color',[.2 .6 .7])
line(time,int_heat./10^(6),'LineWidth',.8,'Color',[.8 .2 .2])
line(time,int_swr./10^(6),'LineWidth',.8,'Color',[.3 .2 .5])
line(time,zeros(size(time)),'LineWidth',.6,'Color',[.3 .3 .3],'LineStyle','--')

% mark the net heat taken by the water column from surface measurement
line([time(end) time(end)],[0 int_total(end)/10^(6)],'Color',...
    [.2 .2 .2],'LineStyle','-','LineWidth',4);
% text(time(end-8500),-int_total(end)/10^6,...
%     ['net heat into the ocean $\sim$ ',num2str(round(int_total(end)/10^(6))),...
%     ' $MJ/m^{2}$'],'fontname','computer modern','Interpreter','latex','fontsize',15)

% figure specification
spec_info.grid_on = 1;
spec_info.xlabel = 'time';
spec_info.ylabel = 'integrated heat ($MJ/m^{2}$)';
spec_info.x_lim = [time(1) time(end)];
spec_info.y_lim = [];
spec_info.lgd = {'total heat exchange','surface heat flux',...
    'short wave radiation'};
spec_info.lgd_pos = [];
spec_info.save_path = './figs/int_heat';

line_annotate(spec_info)
%--------------------------------------------------------------------------


%% Calculation of heat content and temporal derivative --------------------

% observation only goes down to 200m, so we will focus on the upper 200m
% specific heat capacity - gsw_cp0, from GSW toolbox

% depth z index 22 ~ upper 201.4027m

% HC = sum(rho(22:end,:).*(temp(22:end,:)+273.15).*h(22:end,:)*gsw_cp0); %[J/m^2]

cp = 3985; % specific heat capacity from GOTM setup
HC = sum(rho(22:end,:).*(temp(22:end,:)+273.15).*h(22:end,:)*cp); %[J/m^2]
HC_delta = HC(end) - HC(1);
HC_t = gradient(HC,dt*nsave); %[W/m^2]

% from observation
HC_obs = sum(rho(22:end,:).*(temp_obs(22:end,:)+273.15).*h(22:end,:)*gsw_cp0);
HC_delta_obs = HC_obs(end) - HC_obs(1);
HC_t_obs = gradient(HC_obs,dt*nsave);

% from surface heat info
HC_surface = int_total+HC(1);
HC_t_surface = gradient(int_total,dt*nsave);

HC_error = 100*(int_total(end) - HC_delta)/int_total(end); % in percent
%--------------------------------------------------------------------------


%% Heat content time series -----------------------------------------------

figure('position', [0, 0, 900, 300])

line(time,HC./10^(6),'LineWidth',.3,'Color',[.7 .4 .6])
line(time,HC_surface./10^(6),'LineWidth',.3,'Color',[.1 .6 .7])
line(time,HC_obs./10^(6),'LineWidth',.05,'Color',[.4 .8 .6])
line(time,ones(size(time))*HC(1)/10^(6),'LineWidth',.9,'Color',[.3 .3 .3],'LineStyle','--')

% mark the deviation from surface heat input
line([time(end) time(end)],[HC(1)/10^(6) HC_surface(end)./10^(6)],...
    'Color',[.6 .1 .3],'LineStyle','-','LineWidth',3);
% text(time(end-2200),HC(end)/10^(6),...
%     ['$\sim$ ',num2str(round(HC_error,2)),'$\%$'],'color',[.6 .1 .3],...
%     'fontname','computer modern','Interpreter','latex','fontsize',13)

% mark the heat content change relative to observation
line([time(end) time(end)],[HC(1)/10^(6) HC(end)/10^(6)],'Color',...
    [.1 .1 .1],'LineStyle','-','LineWidth',3);
% text(time(end-8500),HC(8800)/10^(6),[turb_method,...
%     ' heat content change $\sim$ ',num2str(round((HC_delta-HC_delta_obs)/10^(6))),...
%     ' $MJ/m^{2}$'],'Color',[.1 .1 .1],'fontname','computer modern',...
%     'Interpreter','latex','fontsize',13)
  
% figure specification
spec_info.lgd = {[turb_method,' HC'],'surface heat exchange','obs. HC'};
spec_info.ylabel = 'heat content ($MJ/m^{2}$)';
spec_info.save_path = './figs/HC';

line_annotate(spec_info)  
%--------------------------------------------------------------------------

%% Plot temporal variation of heat content ---------------------------------

figure('position', [0, 0, 900, 300])

line(time,HC_t_surface,'LineWidth',1,'Color','k')
line(time,HC_t,'LineWidth',1,'Color',[.4 .9 .7])
line(time,zeros(size(time)),'LineWidth',.6,'Color',[.5 .5 .5],'LineStyle',':')
 
% figure specification
spec_info.lgd = {'surface heat exchange rate',...
    ['$\partial_{t}HC$ in ',turb_method]};
spec_info.ylabel = 'temporal heat variation ($W/m^{2}$)';
spec_info.save_path = './figs/HC_t_surf';

line_annotate(spec_info)

%--------
figure('position', [0, 0, 900, 300])

line(time,HC_t_obs,'LineWidth',1,'Color','k')
line(time,HC_t,'LineWidth',1,'Color',[.4 .9 .7])
line(time,zeros(size(time)),'LineWidth',.6,'Color',[.5 .5 .5],'LineStyle',':')

spec_info.lgd = {'$\partial_{t}HC$ in obs.', ...
    ['$\partial_{t}HC$ in ',turb_method]};
spec_info.ylabel = 'temporal heat variation ($W/m^{2}$)';
spec_info.save_path = './figs/HC_t_obs';

line_annotate(spec_info)
