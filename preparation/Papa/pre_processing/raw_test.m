%% raw_test

% Test the the raw data quality in 'met_forcing_p2007.mat' before any
% data manipulation.


%% wind gustiness

bad_w_gust = find(w_gust>100);

% set the upper and bottom limit for plot
v_max = max(w_gust(w_gust<100));
v_min = min(w_gust(w_gust<100));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

% set the value for gap bars in the plot
w_gust_nan = ones(size(w_gust))*NaN;
w_gust_nan(bad_w_gust) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 200])

scatter(time(w_gust<100),w_gust(w_gust<100),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,w_gust_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('zonal wind speed ($$m/s$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/w_gust_sereis','-pdf','-transparent','-painters')


%% zonal wind speed

bad_w_u = find(w_u>100);

v_max = max(w_u(w_u<100));
v_min = min(w_u(w_u<100));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

w_u_nan = ones(size(w_u))*NaN;
w_u_nan(bad_w_u) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 200])

scatter(time(w_u<100),w_u(w_u<100),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,w_u_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('zonal wind speed ($$m/s$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/w_u_sereis','-pdf','-transparent','-painters')

%% meridional wind speed


bad_w_v = find(w_v>100);

v_max = max(w_v(w_v<100));
v_min = min(w_v(w_v<100));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

w_v_nan = ones(size(w_v))*NaN;
w_v_nan(bad_w_v) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 200])

scatter(time(w_v<100),w_v(w_v<100),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,w_v_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('meridional wind speed ($$m/s$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/w_v_sereis','-pdf','-transparent','-painters')

%% air temperature
 
bad_t_air = find(t_air>100);

v_max = max(t_air(t_air<100));
v_min = min(t_air(t_air<100));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

t_air_nan = ones(size(t_air))*NaN;
t_air_nan(bad_t_air) = v_max + 0.05*(v_max-v_min);


figure('position', [0, 0, 1000, 200])

scatter(time(t_air<100),t_air(t_air<100),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,t_air_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('air temperature ($$^{o}C$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/tair_sereis','-pdf','-transparent','-painters')


%% relative humidity


bad_rh = find(rh>1000);

v_max = max(rh(rh<1000));
v_min = min(rh(rh<1000));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

rh_nan = ones(size(rh))*NaN;
rh_nan(bad_rh) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 200])

scatter(time(rh<1000),rh(rh<1000),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,rh_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('relative humidity ($$\%$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/rh_sereis','-pdf','-transparent','-painters')


%% SST
 
bad_sst = find(sst>100);

v_max = max(sst(sst<100));
v_min = min(sst(sst<100));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

sst_nan = ones(size(sst))*NaN;
sst_nan(bad_sst) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 200])

scatter(time(sst<100),sst(sst<100),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,sst_nan,'LineWidth',6,'Color',[.9 .1 .2])

% hold on
% scatter(time(sst_2017),ts_papa,.9,'MarkerEdgeColor',rand(1,3),...
%               'MarkerFaceColor',rand(1,3),'LineWidth',1)
box on
datetick('x','yyyy')
ylabel('SST ($$^{o}C$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')


%export_fig ('./test/sst_sereis','-pdf','-transparent','-painters')


%% downgoing shortwave radiation

bad_Rs = find(Rs>1000);

v_max = max(Rs(Rs<1000));
v_min = min(Rs(Rs<1000));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

Rs_nan = ones(size(Rs))*NaN;
Rs_nan(bad_Rs) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 250])

scatter(time(Rs<1000),Rs(Rs<1000),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,Rs_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('downgoing shortwave radiation ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/Rs_sereis','-pdf','-transparent','-painters')


%% downgoing longwave radiation

bad_Rl = find(Rl>1000);

v_max = max(Rl(Rl<1000));
v_min = min(Rl(Rl<1000));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

Rl_nan = ones(size(Rl))*NaN;
Rl_nan(bad_Rl) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 250])

scatter(time(Rl<1000),Rl(Rl<1000),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time,Rl_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('downgoing longwave radiation ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/Rl_sereis','-pdf','-transparent','-painters')

%% barometric pressure


bad_P = find(P>10000);

v_max = max(P(P<10000));
v_min = min(P(P<10000));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

P_nan = ones(size(P))*NaN;
P_nan(bad_P) = v_max + 0.05*(v_max-v_min);

figure('position', [0, 0, 1000, 200])

%line(time(P<10000),P(P<10000),'LineWidth',.4,'Color',rand(1,3))
scatter(time(P<10000),P(P<10000),.9,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .5 .5],'LineWidth',1)
line(time,P_nan,'LineWidth',7,'Color',[.9 .1 .2])

% hold on
% scatter(time(P_2010),slp_papa_2010,.9,'MarkerEdgeColor',rand(1,3),...
%               'MarkerFaceColor',rand(1,3),'LineWidth',1)
% hold on
% scatter(time(P_2013),slp_papa_2013,.9,'MarkerEdgeColor',rand(1,3),...
%               'MarkerFaceColor',rand(1,3),'LineWidth',1)
box on
datetick('x','yyyy')
ylabel('barometric pressure ($$hPa$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('April 21, 2018')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/p_sereis','-pdf','-transparent','-painters')


%% rain
 
bad_rain = find(rain>100);

v_max = max(rain(rain<100));
v_min = min(rain(rain<100));
y_up = v_max + 0.1*(v_max-v_min);
y_bot = v_min - 0.1*(v_max-v_min);

rain_nan = ones(size(rain))*NaN;
rain_nan(bad_rain) = v_max + 0.05*(v_max-v_min);


figure('position', [0, 0, 1000, 200])

scatter(time_rain(rain<100),rain(rain<100),.9,'MarkerEdgeColor',rand(1,3),...
              'MarkerFaceColor',rand(1,3),'LineWidth',1)
line(time_rain,rain_nan,'LineWidth',6,'Color',[.9 .1 .2])
box on
datetick('x','yyyy')
ylabel('precipitation rate ($$mm/hr$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[datenum('June 08, 2007') datenum('June 15, 2017')],...
'YLim',[y_bot y_up],...
'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

%export_fig ('./test/rain_sereis','-pdf','-transparent','-painters')

%% temperature profile

