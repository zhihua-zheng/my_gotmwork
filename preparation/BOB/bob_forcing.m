%% make_forcing.m
%
% make forcing files for GOTM Bay of Bengal simulation
%
% Zhihua Zheng, UW-APL, Nov. 7 2018
% 
% load meteological observation data from PMEL OWS P station mooring
% -------------------------------------------------------------------------
% %        w_u - x component wind velocity (m/s)
% %        w_v - y component wind velocity (m/s)
% %      w_spd - total wind speed (m/s)
% %      w_dir - wind direction (clockwise to the North), 
% %              in oceanographic sense (degree)
% %     z_wind - wind velocity measurement height (m)
% %      t_air - air temperature (degree centigrade)
% %       z_ta - air temperature measurement height (m)
% %         rh - relative humidity (%)
% %       z_rh - relative humidity measurement height (m)
% %         Rs - downward shortwave radiation (W/m^2)
% %         Rl - downward longwave radiation (W/m^2)
% %       rain - precipitation rate (mm/hr)
% %          P - sea level barometric pressure (hPa)
% %        sst - sea surface temperature (degree centigrade)
% %        sss - sea surface salinity (PSU)
% %        ssd - sea surface potential density anomaly (sigma-theta, kg/m^3)
% %      sprof - salinity profile (PSU)
% %      tprof - temperature profile (degree centigrade)
% %    rhoprof - potential density anomaly profile (kg/m^3)
% %    depth_t - depth for T profile (m)
% %    depth_s - depth for S profile (m), same as depth_rho
% %    cur_spd - total current speed (cm/s)
% %      cur_u - x component current velocity profile (cm/s)
% %      cur_v - y component current velocity profile (cm/s)
% %    cur_dir - current direction (clockwise to the North), (degree)
% %  depth_cur - depth for current profile (m)
% %        lat - mooring latitude (degree)
% %        lon - mooring longitude (degree)
% %       time - datenumbers for measurements (UTC)
% %       date - date strings for measurements (UTC)
% %   time_cur - datenumbers for current profile measurements (UTC)
% %   date_cur - date strings for current profile measurements (UTC)
% %  time_rain - datenumbers for precipitation measurements (UTC)
% %  date_rain - date strings for precipitation measurements (UTC)
% %     dn2007 - datenumber for start time of time series
%   dn2007_cur - datenumber for start time of current profile time series
%  dn2007_rain - datenumber for start time of precipitation time series

% load wave data from OWS P station CDIP buoy
% -------------------------------------------------------------------------
% %   wave_bw - frequency bandwidth (hertz)
% % wave_freq - band center frequency (hertz)
% % wave_mdir - band mean direction (counter-clockwise to the East), 
% %             in oceanographic sense (degree)
% % wave_spec - band energy density (m^2/hertz)
% % wave_date - date strings for wave measurements (UTC)
% % wave_time - datenumbers for wave measurements (UTC)

%% Load the Met. forcing 

clear

forcing_dir = '~/Documents/Study/Graduate_research/data_raw/BOB/';

% load the whole workspace as a struct
% met_data = load('met_forcing_p2007.mat'); 

load([forcing_dir 'Forcing.mat'])
% load('wave_p2010.mat')

% set the negative values of rain data to zeros
rain(rain<0) = 0;

%% pick out good data

raw_test; % initial examination of data quality

%% organize data 

% After preliminary analysis, noticed that most data is ruined for a
% long period between Nov. 2008 and Jun. 2009. We decided to truncate the 
% data and pick out the data after Jun. 2009

% Wanted to include rain data in the computation of flux, but the
% precipitation measurement is intermittent, hence rain data is ignored,
% and the sensible heat flux due to rain (small) is not included.

inx2009 = find(date=='2009/01/01 00:00:00');
inx2010 = find(date=='2010/01/01 00:00:00');

% loop through variables to find the latest good data beginning time 
% (pre_day)
vars = {'w_u','w_v','sst','P','t_air','rh','Rs','Rl','w_dir','rain',...
    'sprof','tprof'};

for i = 1:length(vars)
   
    % large gap from '2009/01/01 01:00:00' to '2010/01/01 00:00:00',
    % hence pick out data from that period for an examination
    if i < 11 
        var_test = met_data.(vars{i})(inx2009+1:inx2010);
    else
        var_test = met_data.(vars{i})(1,inx2009+1:inx2010);
    end

    % preferred starting time index
    if i < 2
        
        pre_day = find(var_test>10000,1,'last') + inx2009 + 1; 
    else
        tmp = find(var_test>10000,1,'last') + inx2009 + 1;

        if tmp > pre_day
            pre_day = tmp;
        end
    end   
end

clear tmp inx2009 inx2010

% Abandon all the data before pre_day, and linearly interpolate the small
% gaps

% Note that barometric pressure, rain rate and sea surface temperature have 
% relative long periods of missing data, hence the large gaps in P and sst 
% time series are filled by interpolating the reanalysis data from MERRA2 
% in space and time to the ows_papa mooring location and data time.

merra_fill; % subroutine to fill the other large gaps

% truncation (and interpolation)
sst_r = sst_r(pre_day:end);
P_r = P_r(pre_day:end);
rain_r = rain_r(pre_day:end);

% w_spd_r = interp1(time(w_spd<100),w_spd(w_spd<100),time(pre_day:end));
% extra large gap at the end, compared to w_u and w_v data...

w_u_r = interp1(time(w_u<100),w_u(w_u<100),time(pre_day:end));
w_v_r = interp1(time(w_v<100),w_v(w_v<100),time(pre_day:end));
w_spd_r = sqrt(w_u_r.^2 + w_v_r.^2);
sss_r = interp1(time(sss<100),sss(sss<100),time(pre_day:end));
t_air_r = interp1(time(t_air<50),t_air(t_air<50),time(pre_day:end));
rh_r = interp1(time(rh<150),rh(rh<150),time(pre_day:end));
Rs_r = interp1(time(Rs<1000),Rs(Rs<1000),time(pre_day:end));
Rl_r = interp1(time(Rl<1000),Rl(Rl<1000),time(pre_day:end));
w_dir_r = interp1(time(w_dir<=360),w_dir(w_dir<=360),time(pre_day:end));
time_r = time(pre_day:end); % truncated datenumbers for measurements (UTC)
date_r = date(pre_day:end); % truncated strings for measurements (UTC)


% TS profile data truncation and interpolation
sprof(sprof>1000) = NaN;
tprof(tprof>1000) = NaN;

sprof_tr = sprof(:,pre_day:end);
tprof_tr = tprof(:,pre_day:end);

% Absolute salinity and potential temperature conversion (TEOS-10)

[T_s, Z_s] = meshgrid(time_r,depth_s);
F = scatteredInterpolant(T_s(~isnan(sprof_tr)),Z_s(~isnan(sprof_tr)),...
    sprof_tr(~isnan(sprof_tr)),'linear'); % interpolation object
sprof_tr_i = F(T_s,Z_s);

sprof_r = gsw_SA_from_SP(sprof_tr_i,Z_s,lon,lat);
% henceforth absolute salinity


[T_t, Z_t] = meshgrid(time_r,depth_t);

% temperature at S-profile depth level
F = scatteredInterpolant(T_t(~isnan(tprof_tr)),Z_t(~isnan(tprof_tr)),...
    tprof_tr(~isnan(tprof_tr)),'linear'); 
tprof_depth_s = F(T_s,Z_s);

tprof_r = gsw_pt0_from_t(sprof_r,tprof_depth_s,Z_s); 
% henceforth potential temperature

%----- velocity data is extremely bad
% [T, Z] = meshgrid(time_prof_r,depth_cur);
% cur_u_r = griddata(T(cur_u<100),Z(cur_u<100),cur_u(cur_u<100),T,Z,'linear');


%% compute surface flux

% momentum flux = - surface wind stress
%
% net heatflux (positive in) = - latent heat flux (positive out) - sensible 
% heat flux (positive out) + net longwave radiation (positive in)


% use wind-speed dependent formulation, without wave information
A = coare35vn(w_spd_r,z_wind,t_air_r,z_ta,rh_r,z_rh,P_r,sst_r,Rs_r,...
    Rl_r,lat,NaN,rain_r,NaN,NaN);
% note we should use relative velocity when I have time to do that, but
% does it mattet? I don't have good surface current measurement record.

tau = A(:,2);
hsb = A(:,3); % sensible heat flux
hlb = A(:,4); % latent heat flux
u10 = A(:,29); % wind speed adjusted to 10 m

w_cos = w_u_r./w_spd_r;
w_sin = w_v_r./w_spd_r;

% surface momentum flux in x, y direction
tau_x = tau.*w_cos;
tau_y = tau.*w_sin;

% u10 in x, y direction
u10_x = u10.*w_cos;
u10_y = u10.*w_sin;

%------ get the date vector and decimal yearday, adjusted for leap year

date_vec = datevec(char(date_r)); 
% lp = leapyear(date_vec(:,1));  % leap year index
% date_vec = [date_vec lp];
% yd = yearday(date_vec(:,2),date_vec(:,3),date_vec(:,7)); 
% yearday - vectorization issue with function yearday
% date_vec = [date_vec yd];

yd = date2doy(time_r)-1; % using function date2doy from File Exchange
%-----------------------

% compute net short wave heat flux
nsw = swhf(yd,date_vec(:,1),(360-lon)*ones(size(time_r)),...
    lat*ones(size(time_r)),Rs_r);
%---- Need to change the longitude into format of West positive degree
%---- Note the operating time in the subroutine of swhf, soradna1 (no-sky
% solar radiation) is in UTC format, therefore the time zone shifting must 
% be perfomed after here.

% compute net long wave heat flux
nlw = lwhf(sst_r,Rl_r,Rs_r);

% compute net surfac heat flux (shortwave heat flux and advective heat
% flux are omitted)
hf = -hlb - hsb + nlw;


%% Compute Wave Spectrum (frequency bin width weighted)

% observed band directional component-x and component-y
xcmp = cos(wave_mdir*pi/180); 
ycmp = sin(wave_mdir*pi/180);

% wave spectrum
spec_h2 = wave_spec.*repmat(wave_bw,1,length(wave_time));

%% Time Conversion

% USE UTC TIME!!
date_r = string(datestr(time_r, 'yyyy-mm-dd HH:MM:SS'));

% % correction for time zone shift
% UTC_time = datetime(time_r,'ConvertFrom','datenum','TimeZone','UTC');
% Papa_time = datetime(UTC_time,'TimeZone','-10:00');
% 
% % date number, date string and date vector for papa station local time
% time_r = datenum(Papa_time);
% date_r = datestr(Papa_time,'yyyy/mm/dd HH:MM:SS'); % leave it as 'character'
% date_vec = datevec(date_r); 
% 
% % do it again for time_prof_r
% UTC_time = datetime(time_prof_r,'ConvertFrom','datenum','TimeZone','UTC');
% Papa_time = datetime(UTC_time,'TimeZone','-10:00');
% time_prof_r = datenum(Papa_time);
% date_prof_r = string(datestr(Papa_time,'yyyy/mm/dd HH:MM:SS'));

% Hereafter all the truncated time related variables are in local time zone.

% An updated workspace 'gotm_input_p2007.mat' is hence saved.
% save('gotm_input_p2007.mat');

%% Diurnal SWR check

diurnal = ones(3000,24)*NaN; % rearrange for every hour (0-23)

for i = 1:24
    
    inx = find(date_vec(:,4) == i-1);
    diurnal(1:size(inx),i) = Rs_r(inx);
    
end

one_day = nanmean(diurnal);

plot(one_day)
% plot(circshift(one_day,-10))


%% momentum flux file

fileID = fopen('../setup_files/tau_file.dat','w');
H = [cellstr(date_r) num2cell(tau_x) num2cell(tau_y)];
formatSpec = '%s  % 8.6e % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

% one dot (.) - present folder, two dots (..) - parent of current folder

%% U10 file

fileID = fopen('../setup_files/u10_file.dat','w');
H = [cellstr(date_r) num2cell(u10_x) num2cell(u10_y)];
formatSpec = '%s  % 8.6e % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% heat flux file

fileID = fopen('../setup_files/heatflux_file.dat','w');
H = [cellstr(date_r) num2cell(hf)];
formatSpec = '%s   % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% precipitation file

fileID = fopen('../setup_files/precip_file.dat','w');
H = [cellstr(date_r) num2cell(rain_r)];
formatSpec = '%s   % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% sea surface temparature (sst) file

fileID = fopen('../setup_files/sst_file.dat','w');
H = [cellstr(date_r) num2cell(sst_r)];
formatSpec = '%s %6.3f\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% sea surface salinity (sss) file

fileID = fopen('../setup_files/sss_file.dat','w');
H = [cellstr(date_r) num2cell(sss_r)];
formatSpec = '%s %6.3f\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% net short wave radiation (swr) file

fileID = fopen('../setup_files/swr_file.dat','w');
H = [cellstr(date_r) num2cell(nsw)];
formatSpec = '%s  % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% salinity profile

fileID = fopen('../setup_files/sprof_file.dat','w');

for i = 1:size(time_r,1)
    
    fprintf(fileID,'%s  18  2\n',date_r(i));
    fprintf(fileID,'%6.1f   %9.6f\n',[-depth_s'; sprof_r(:,i)']);
end

fclose(fileID);

%% temperature profile

fileID = fopen('../setup_files/tprof_file.dat','w');

for i = 1:size(time_r,1)
    
    fprintf(fileID,'%s  18  2\n',date_r(i));
    fprintf(fileID,'%6.1f   %9.6f\n',[-depth_s'; tprof_r(:,i)']);
end

fclose(fileID);

%% velocity profile

% fileID = fopen('../setup_files/cur_prof.dat','w');
% 
% for i = 1:size(time_prof_r,1)
%     
%     fprintf(fileID,'%s  18 2\n',date_prof_r(i));
%     fprintf(fileID,'% 9.4f   % 8.4f\n',[depth_t'; tprof_r(:,i)']);
% end
% 
% fclose(fileID);
% 
% copyfile('../setup_files/cur_prof.dat', '../setup_files/cur_prof_file.dat');

%% wave spectrum file


fileID = fopen('../setup_files/spec_file.dat','w');

for i = 1:size(wave_time,1)
    
    fprintf(fileID,'%s  64  1\n',wave_date(i));
    fprintf(fileID,'%10.5f   %8.6e   %8.6e   %8.6e\n',...
        [wave_freq'; spec_h2(:,i)'; xcmp(:,i)'; ycmp(:,i)']);
end

fclose(fileID);

gzip('../setup_files/spec_file.dat')
delete('../setup_files/spec_file.dat')

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
  
  
% salinity profiles
figure('position', [0, 0, 900, 250])

cnum = 15;
CL = [min(min(sprof_r)) max(max(sprof_r))];
conts = linspace(CL(1),CL(2),cnum);
cmocean('-haline')
[T, Z] = meshgrid(time_r,depth_s);
contourf(T,Z,sprof_r,conts,'LineWidth',0.01,'LineStyle','none')
  caxis(CL);
  box on
  axis ij
  datetick('x','yyyy')
  ylabel('depth ($$m$$)', 'fontname', 'computer modern', ...
      'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', ...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'fontsize',...
      11,'fontname','computer modern','TickLabelInterpreter', 'latex')
  
  h = colorbar('EastOutside');
  h.Label.String = 'salinity ($$psu$$)';
  h.Label.Interpreter = 'latex';
  h.Label.FontName = 'computer modern';
  h.Label.FontSize = 14;
  set(h,'TickLabelInterpreter','latex','fontsize',9);

  export_fig ('./figs/prof_sal','-pdf','-transparent','-painters')
  
  
  
% temperature profiles
figure('position', [0, 0, 900, 200])

cnum = 15;
CL = [min(min(tprof_r)) max(max(tprof_r))];
conts = linspace(CL(1),CL(2),cnum);
cmocean('matter')
[T, Z] = meshgrid(time_r,depth_t);
contourf(T,Z,tprof_r,conts,'LineWidth',0.01,'LineStyle','none')
  caxis(CL);
  box on
  axis ij
  datetick('x','yyyy')
  ylabel('depth ($$m$$)', 'fontname', 'computer modern', ...
      'fontsize', 14,'Interpreter', 'latex')
  xlabel('time', 'fontname', 'computer modern', 'fontsize', ...
      14,'Interpreter', 'latex')
  setDateAxes(gca,'XLim',[time_r(1) time_r(end)],'fontsize',...
      11,'fontname','computer modern','TickLabelInterpreter', 'latex')
  h = colorbar('EastOutside');
  h.Label.String = 'temperature ($$^{\circ}C$$)';
  h.Label.Interpreter = 'latex';
  h.Label.FontName = 'computer modern';
  h.Label.FontSize = 14;
  set(h,'TickLabelInterpreter','latex','fontsize',9);
  
  export_fig ('./figs/prof_temp','-pdf','-transparent','-painters')
  
 