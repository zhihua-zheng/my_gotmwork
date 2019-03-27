%% papa_forcing.m
%
% make forcing files for new GOTM OCSPapa test case
%
% Zhihua Zheng, UW-APL, updated on Mar. 25 2019
% 
% load meteological observation data from PMEL OCSP station mooring
% -------------------------------------------------------------------------
% %        w_u - x component wind velocity (m/s)
% %        w_v - y component wind velocity (m/s)
% %      w_spd - total wind speed (m/s)
% %      w_dir - wind direction (clockwise to the North), 
% %              in oceanographic sense (degree)
% %     z_wind - wind velocity measurement height (m)
% %         Ta - air temperature (degree centigrade)
% %       z_Ta - air temperature measurement height (m)
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
% %      dprof - potential density anomaly profile (kg/m^3)
% %    depth_t - depth for T profile (m)
% %   depth_sd - depth for S profile (m), same as depth_p
% %      uprof - x component current velocity profile (cm/s)
% %      vprof - y component current velocity profile (cm/s)
% %  depth_cur - depth for current profile (m)
% %        lat - mooring latitude (degree)
% %        lon - mooring longitude (degree)
% %       time - datenumbers for Met. measurements (UTC)

% load wave data from OCSP station CDIP buoy
% -------------------------------------------------------------------------
% %   wave_bw - frequency bandwidth (hertz)
% % wave_freq - band center frequency (hertz)
% % wave_mdir - band mean direction (counter-clockwise to the East), 
% %             in oceanographic sense (degree)
% % wave_spec - band energy density (m^2/hertz)
% % wave_date - date strings for wave measurements (UTC)
% % wave_time - datenumbers for wave measurements (UTC)

%% Load the Met. forcing and wave data

clear
data_dir = '~/Documents/Study/Grad_research/data/OCSP/';

% load the whole workspace as a struct
met_data = load([data_dir,'Mooring/Met_IC_ocsp.mat']); 

load([data_dir,'Mooring/Met_IC_ocsp.mat'])
load([data_dir,'CDIP_Wave/wave_2010.mat'])

% set the negative values of rain data to outlier
rain(rain<0) = 1e35;

%% initial look

raw_test; % initial examination of data quality

%% organize data 

% After preliminary analysis, noticed that most data is ruined for a
% long period before sometime in 2009. We decided to pcik a good starting 
% time in 2009 and only pick out the data after that

inx2009 = find(time == datenum(2009,1,1,0,0,0));
inx2010 = find(time == datenum(2010,1,1,0,0,0));

% loop through variables to find the latest date (pre_day), after whihc all
% the data have good quality 
vars = {'w_u','w_v','sst','P','Ta','rh','Rs','Rl','w_dir','rain',...
    'sprof','tprof'};

for i = 1:length(vars)
   
    % large gap from '2009/01/01 01:00:00' to '2010/01/01 00:00:00',
    % hence pick out data from that period for an examination
    if i < 11 
        var_test = met_data.(vars{i})(inx2009+1:inx2010);
    else
        var_test = met_data.(vars{i})(1,inx2009+1:inx2010);
    end

    % preferred truncating time index (after last bad flag)
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

% Abandon the data before pre_day, and linearly interpolate the small gaps

% Note that barometric pressure, rain rate and sea surface temperature have 
% relative long periods of missing data, hence the large gaps in P and sst 
% time series are filled by interpolating the reanalysis data MERRA2 in 
% space and time to the OCSP mooring location and time.

% function to fill other large gaps
[P_r,sst_r,rain_r] = merra_fill(P,sst,rain,lat,lon,time); 

% truncation for COARE inputs
sst_r  = sst_r(pre_day:end);
P_r    = P_r(pre_day:end);
rain_r = rain_r(pre_day:end);

time_r = time(pre_day:end);

% truncation and linear interpolation for COARE inputs
w_u_r  = interp1(time(w_u<100),w_u(w_u<100),time_r);
w_v_r  = interp1(time(w_v<100),w_v(w_v<100),time_r);
Ta_r   = interp1(time(Ta<50),Ta(Ta<50),time_r);
rh_r   = interp1(time(rh<150),rh(rh<150),time_r);
Rs_r   = interp1(time(Rs<1000),Rs(Rs<1000),time_r);
Rl_r   = interp1(time(Rl<1000),Rl(Rl<1000),time_r);

w_spd_r = sqrt(w_u_r.^2 + w_v_r.^2);

% truncation and linear interpolation for other variables
sss(sss>100) = NaN;
ssd(ssd>100) = NaN;
sss_r = interp1(time(~isnan(sss)),sss(~isnan(sss)),time_r);
ssd_r = interp1(time(~isnan(ssd)),ssd(~isnan(ssd)),time_r);

% T-S-D profile truncation and interpolation
tprof(tprof>100) = NaN;
sprof(sprof>100) = NaN;
dprof(dprof>100) = NaN; % density anomaly

tprof_r = tprof(:,pre_day:end);
sprof_r = sprof(:,pre_day:end);
dprof_r = dprof(:,pre_day:end);

tprof_r = vert_fill(tprof_r,-depth_t, time_r);
sprof_r = vert_fill(sprof_r,-depth_sd,time_r);
dprof_r = vert_fill(dprof_r,-depth_sd,time_r);

% Absolute salinity and potential temperature conversion (TEOS-10)
% GSW toolbox uses pressure, hence positive depth

SAprof = gsw_SA_from_SP(sprof_r,depth_sd,lon,lat); % absolute salinity

% Project tprof to depth_sd
tprof_sd = nan(size(SAprof));
for j = 1:length(time_r)
    
    tprof_sd(:,j) = interp1(-depth_t,tprof_r(:,j),-depth_sd);
end

PTprof = gsw_pt0_from_t(SAprof,tprof_sd,depth_sd); % potential temperature

z_tsd = -depth_sd;

% save T-S-D profiles
save('~/Documents/Study/Grad_research/data/OCSP/Mooring/profs.mat',...
    'tprof','tprof_r','PTprof','sprof','SAprof','dprof','dprof_r',...
    'depth_t','z_tsd','time','time_r')

clear depth_sd depth_t

%% Time Conversion

% USE UTC TIME!!
date_r = string(datestr(time_r, 'yyyy-mm-dd HH:MM:SS'));

%% compute surface flux

% net heatflux (positive in) = latent heat flux (positive in) + sensible 
% heat flux (positive in) + net longwave radiation (positive in)

% call COARE (3.5) bulk formula
A = coare35vn(w_spd_r,z_wind,Ta_r,z_Ta,rh_r,z_rh,P_r,sst_r,Rs_r,...
    Rl_r,lat,NaN,rain_r,NaN,NaN);
% note we should use relative velocity when I have time to do that, but
% does it matter? I don't have good surface current measurement record.

tau = A(:,2); % wind stress [N/m^2]
hsb = A(:,3); % sensible heat     flux into the ocean [W/m^2]
hlb = A(:,4); % latent   heat     flux into the ocean [W/m^2]
hbb = A(:,5); % surface  buoyancy flux into the ocean [W/m^2]

w_cos = w_u_r./w_spd_r;
w_sin = w_v_r./w_spd_r;

% surface momentum flux in x, y direction
tau_x = tau.*w_cos;
tau_y = tau.*w_sin;

date_vec = datevec(char(date_r)); 
yd = date2doy(time_r)-1; % yearday (starts at 1) from 'date2doy'
%-----------------------

% compute net short wave heat flux, longitude: West positive degree
nsw = swhf(yd,date_vec(:,1),(360-lon)*ones(time_r),...
    lat*ones(time_r),Rs_r);

% compute net longwave heat flux into the ocean
nlw = lwhf(sst_r,Rl_r,Rs_r);

% compute net surfac heat flux into the ocean
hf = hlb + hsb + nlw;

% save Meteological variables
save('~/Documents/Study/Grad_research/data/OCSP/Mooring/MF.mat',...
    'tau_x','tau_y','nsw','hf','hbb',...
    'sst_r','sss_r','ssd_r','rain_r','time_r')

%% Compute Wave Spectrum (frequency bin width weighted)

% observed band directional component-x and component-y
xcmp = cos(wave_mdir*pi/180); 
ycmp = sin(wave_mdir*pi/180);

% wave spectrum
spec_h2 = wave_spec.*repmat(wave_bw,1,length(wave_time));

%% Diurnal SWR check

diurnal = nan(3000,24); % rearrange for every hour (0-23)

for i = 1:24
    
    inx = find(date_vec(:,4) == i-1);
    diurnal(1:size(inx),i) = Rs_r(inx);
    
end

one_day = nanmean(diurnal);

plot(one_day)
% plot(circshift(one_day,-10))

%% save

basecase = '../../data/OCSPapa_20070608-20190325/';

%% momentum flux file

fileID = fopen([basecase,'tau_file.dat'],'w');
H = [cellstr(date_r) num2cell(tau_x) num2cell(tau_y)];
formatSpec = '%s  % 8.6e % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

% one dot (.) - present folder, two dots (..) - parent of current folder

%% U10 file

% fileID = fopen([basecase,'u10_file.dat'],'w');
% H = [cellstr(date_r) num2cell(u10_x) num2cell(u10_y)];
% formatSpec = '%s  % 8.6e % 8.6e\n';
% 
% for i = 1:size(H,1)
%     fprintf(fileID,formatSpec,H{i,:});
% end
% 
% fclose(fileID);

%% heat flux file

fileID = fopen([basecase,'heatflux_file.dat'],'w');
H = [cellstr(date_r) num2cell(hf)];
formatSpec = '%s   % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% precipitation file

fileID = fopen([basecase,'precip_file.dat'],'w');
H = [cellstr(date_r) num2cell(rain_r)];
formatSpec = '%s   % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% sea surface temparature (sst) file

fileID = fopen([basecase,'sst_file.dat'],'w');
H = [cellstr(date_r) num2cell(sst_r)];
formatSpec = '%s %6.3f\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% sea surface salinity (sss) file

fileID = fopen([basecase,'sss_file.dat'],'w');
H = [cellstr(date_r) num2cell(sss_r)];
formatSpec = '%s %6.3f\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% net short wave radiation (swr) file

fileID = fopen([basecase,'swr_file.dat'],'w');
H = [cellstr(date_r) num2cell(nsw)];
formatSpec = '%s  % 8.6e\n';

for i = 1:size(H,1)
    fprintf(fileID,formatSpec,H{i,:});
end

fclose(fileID);

%% salinity profile

fileID = fopen([basecase,'sprof_file.dat'],'w');

for i = 1:size(time_r,1)
    
    fprintf(fileID,'%s  18  2\n',date_r(i));
    fprintf(fileID,'%6.1f   %9.6f\n',[z_tsd'; SAprof(:,i)']);
end

fclose(fileID);

%% temperature profile

fileID = fopen([basecase,'tprof_file.dat'],'w');

for i = 1:size(time_r,1)
    
    fprintf(fileID,'%s  18  2\n',date_r(i));
    fprintf(fileID,'%6.1f   %9.6f\n',[z_tsd'; PTprof(:,i)']);
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


fileID = fopen([basecase,'spec_file.dat'],'w');

for i = 1:size(wave_time,1)
    
    fprintf(fileID,'%s  64  1\n',wave_date(i));
    fprintf(fileID,'%10.5f   %8.6e   %8.6e   %8.6e\n',...
        [wave_freq'; spec_h2(:,i)'; xcmp(:,i)'; ycmp(:,i)']);
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
  
 