%% visualization of flux time series for old osw_p case

% heat flux
figure('position', [0, 0, 1000, 200])
line(time_old,cell2mat(heatflux_old(:,2)),'LineWidth',.4,'Color',[.2 .7 .4])
line(time_old,zeros(size(time_old)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('surface heat flux ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 1975') datenum('June 16, 1981')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')
  
  export_fig ('sur_hf_old','-pdf','-transparent','-painters')


  
% momentum flux
figure('position', [0, 0, 1000, 200])
line(time_old,cell2mat(momentumflux_old(:,2)),'LineWidth',.4,'Color',[.6 .7 .4])
line(time_old,zeros(size(time_old)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('x - momentum flux ($$N/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 1975') datenum('June 16, 1981')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')
  
  export_fig ('sur_mfx_old','-pdf','-transparent','-painters')
%-----------
figure('position', [0, 0, 1000, 200])
line(time_old,cell2mat(momentumflux_old(:,3)),'LineWidth',.4,'Color',[.6 .7 .4])
line(time_old,zeros(size(time_old)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('y - momentum flux ($$N/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 1975') datenum('June 16, 1981')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')
  
  export_fig ('sur_mfy_old','-pdf','-transparent','-painters')

  
 
  
  
% net short wave radiation
figure('position', [0, 0, 1000, 200])
line(time_old,cell2mat(swr_old(:,2)),'LineWidth',.4,'Color',[.6 .2 .4])
%line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('net shortwave radiation ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 1975') datenum('June 16, 1981')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')
  
  export_fig ('sur_nsw_old','-pdf','-transparent','-painters')
  
  
  
%%


figure('position', [0, 0, 1000, 200])

subplot(4,1,1)
line(time_r,-hlb,'LineWidth',.4,'Color',rand(1,3))
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('latent ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 2010') datenum('June 16, 2016')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')


subplot(4,1,2)
line(time_r,-hsb,'LineWidth',.4,'Color',rand(1,3))
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('sensible ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 2010') datenum('June 16, 2016')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

subplot(4,1,3)
line(time_r,nlw,'LineWidth',.4,'Color',rand(1,3))
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('longwave ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 2010') datenum('June 16, 2016')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')
  
  
subplot(4,1,4)
line(time_r,nsw,'LineWidth',.4,'Color',rand(1,3))
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('shortwave ($$W/m^2$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[datenum('June 15, 2010') datenum('June 16, 2016')],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')
 
export_fig ('hf_separate','-pdf','-transparent','-painters')
