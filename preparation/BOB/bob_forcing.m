%% bob_forcing.m
%
% make forcing files for Bay of Bengal simulation in GOTM
%
% Zhihua Zheng, UW-APL, Nov. 7 2018
%
% ------ load meteological forcing data -----------------------------------
%
%  FWFlux: freshwater flux (evaporation is negative), unit? assume [m/s]
%  rhf: short-wave heat flux
%  thf: sensible + latent + long-wave heat flux
%  taux: x-direction wind stress [708x1 double]
%  tauy: y-direction wind stress[708x1 double]
%  time: datetime [708x1 double]
%  yd: year day [708x1 double]
%  sst: sea surface temperature [36x1 double] size?
%
% ------ load T-S profile data -----------------------------------
%
%  SG: salinity profiles [z,t]
%  TG: potential temperature profiles [z,t]
%  SigG: density anomaly profiles [z,t]
%  ZG: 2-D z coordinates [z,t]
%  YG: 2-D time coordinates as year day starting from 1 [z,t]
%  SxG: x-direction shear profile? [z,t]
%  SyG: y-direction shear profile? [z,t]

%% Load the Met. forcing and T-S profiles

clear

data_dir = '~/Documents/Study/Graduate_research/data_raw/BOB/';

load([data_dir 'Forcing.mat'])
forc.date = string(datestr(forc.time, 'yyyy/mm/dd HH:MM:SS'));

load([data_dir 'GriddedData.mat']);

basecase = '../../data/BOB/';

%% interpolate the T-S profiles

[nz, ntime] = size(SG);

% fill in deep holes
for iz = 1:nz
    b = find(isnan(SG(iz,:))); % index for NaN
    if ~isempty(b)
        g = find(~isnan(SG(iz,:))); % index for good data
        SG(iz,b) = interp1(YG(iz,g),SG(iz,g),YG(iz,b));
        TG(iz,b) = interp1(YG(iz,g),TG(iz,g),YG(iz,b));
        SigG(iz,b) = interp1(YG(iz,g),SigG(iz,g),YG(iz,b));
    end
end

MG = YG-1+datenum('01-Jan-2016');  % matlab time
MG_date = string(datestr(MG(1,:), 'yyyy/mm/dd HH:MM:SS'));
TTG = sw_temp(SG,TG,-ZG,0);  % get temp from potemp

%% salinity profile

if exist([basecase,'sprof_file.dat'],'file')
    delete([basecase,'sprof_file.dat'])
end

fileID = fopen([basecase,'sprof_file.dat'],'w');

for i = 1:ntime

    fprintf(fileID,'%s  81  2\n',MG_date(i));
    fprintf(fileID,'%6.1f   %9.6f\n',[flip(ZG(:,i))'; flip(SG(:,i))']);
end

fclose(fileID);

%% temperature profile

if exist([basecase,'tprof_file.dat'],'file')
    delete([basecase,'tprof_file.dat'])
end

fileID = fopen([basecase,'tprof_file.dat'],'w');

for i = 1:ntime

    fprintf(fileID,'%s  81  2\n',MG_date(i));
    fprintf(fileID,'%6.1f   %9.6f\n',[flip(ZG(:,i))'; flip(TTG(:,i))']);
end

fclose(fileID);

%% momentum flux file

fileID = fopen([basecase,'tau_file.dat'],'w');
H = [cellstr(forc.date) num2cell(forc.taux') num2cell(forc.tauy')];
formatSpec = '%s  % 8.6e % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

% one dot (.) - present folder, two dots (..) - parent of current folder

%% heat flux file

fileID = fopen([basecase,'heatflux_file.dat'],'w');
H = [cellstr(forc.date) num2cell(forc.thf')];
formatSpec = '%s   % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% precipitation file

% assume the precip file is the same as freshwater flux, and unit is m/s

forc.FWFlux = forc.FWFlux*1000/(1/3600); % [mm/hr]

fileID = fopen([basecase,'precip_file.dat'],'w');
H = [cellstr(forc.date) num2cell(forc.FWFlux')];
formatSpec = '%s   % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);


%% net short wave radiation (swr) file

fileID = fopen([basecase,'swr_file.dat'],'w');
H = [cellstr(forc.date) num2cell(forc.rhf')];
formatSpec = '%s  % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% visualization of flux time series for comparison

% heat flux
figure('position', [0, 0, 1000, 200])
line(time_r,hf,'LineWidth',.4,'Color',[.2 .7 .4])
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('surface heat flux ($$W/m^2$$)', 'fontname',...
      'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize',...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'fontsize',11,...
      'fontname','computer modern','TickLabelInterpreter', 'latex')

  export_fig ('./figs/sur_hf','-pdf','-transparent','-painters')


% precipitation
figure('position', [0, 0, 1000, 200])
line(time_r,rain_r,'LineWidth',.4,'Color',[.6 .2 .5])
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('precipitation rate ($$mm/hr$$)', 'fontname', ...
      'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', ...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'fontsize',11,...
      'fontname','computer modern','TickLabelInterpreter', 'latex')

  export_fig ('./figs/precip','-pdf','-transparent','-painters')


% momentum flux
figure('position', [0, 0, 1000, 200])
line(time_r,tau_x,'LineWidth',.4,'Color',[.6 .7 .4])
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('x - momentum flux ($$N/m^2$$)', 'fontname', ...
      'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', ...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'YLim',[-3 3],...
      'fontsize',11,'fontname','computer modern',...
      'TickLabelInterpreter', 'latex')

  export_fig ('./figs/sur_mfx','-pdf','-transparent','-painters')
%-----------
figure('position', [0, 0, 1000, 200])
line(time_r,tau_y,'LineWidth',.4,'Color',[.6 .7 .4])
line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('y - momentum flux ($$N/m^2$$)', 'fontname', ...
      'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', ...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'YLim',[-3 3],...
      'fontsize',11,'fontname','computer modern',...
      'TickLabelInterpreter', 'latex')

  export_fig ('./figs/sur_mfy','-pdf','-transparent','-painters')



% sst
figure('position', [0, 0, 1000, 200])
line(time_r,sst_r,'LineWidth',.4,'Color',[.8 .4 .2])
%line(time_r,zeros(size(time_r)),'LineWidth',.6,'Color',[.3 .4 .3])
  box on
  datetick('x','yyyy')
  ylabel('sea surface temperature ($$^{\circ}C$$)', 'fontname',...
  'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', ...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'fontsize',11,...
      'fontname','computer modern','TickLabelInterpreter', 'latex')

  export_fig ('./figs/sur_t','-pdf','-transparent','-painters')




% net short wave radiation
figure('position', [0, 0, 1000, 200])
line(time_r,nsw,'LineWidth',.4,'Color',[.6 .2 .4])
  box on
  datetick('x','yyyy')
  ylabel('net shortwave radiation ($$W/m^2$$)', 'fontname', ...
      'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', ...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'fontsize',11,...
      'fontname','computer modern','TickLabelInterpreter', 'latex')

  export_fig ('./figs/sur_nsw','-pdf','-transparent','-painters')
