%% papa_forcing.m
%
% make forcing files for modern GOTM OCSPapa test case
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
% %          P - sea level barometric pressure (hPa = mbar)
% %        sst - sea surface temperature (C), measured at -1 m
% %        sss - sea surface salinity (PSU), measured at -1 m
% %        ssd - sea surface potential density anomaly (sigma-theta, kg/m^3), measured at -1 m
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
% %    wave_bw - frequency bandwidth (hertz)
% %  wave_freq - band center frequency (hertz)
% %  wave_mdir - band mean direction (counter-clockwise to the East), 
% %              in oceanographic sense (degree)
% %  wave_spec - band energy density (m^2/hertz)
% %    wave_a1 - band a1 Fourier coefficients
% %    wave_b1 - band b1 Fourier coefficients
% %  wave_date - date strings for wave measurements (UTC)
% %  wave_time - datenumbers for wave measurements (UTC)

%% Load the Met. forcing and wave data

clear
ocsp_dir = '~/GDrive/UW/Research/Data/OCSP/';

% load the whole workspace as a struct
met_data = load([ocsp_dir,'Mooring/Met_2007_2019_raw.mat']); 

load([ocsp_dir,'Mooring/Met_2007_2019_raw.mat'])
load([ocsp_dir,'CDIP/Wave_2010.mat'])

% set the negative values of rain data to outlier
rain(rain<0) = 1e35;

%% initial look

raw_test; % initial examination of data quality

%% organize data 

% After preliminary analysis, noticed that most data is ruined for a
% long period before 2009. I decided to pick a good starting point in 2009
% and only use the data after that

inx2009 = find(time == datenum(2009,1,1,0,0,0));
inx2010 = find(time == datenum(2010,1,1,0,0,0));

% Loop through variables to find the latest date (pre_day), after which all
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
time_r = time(pre_day:end);

% truncation and linear interpolation for COARE inputs
w_u_r   = interp1(time(w_u<100),w_u(w_u<100),time_r);
w_v_r   = interp1(time(w_v<100),w_v(w_v<100),time_r);
Ta_r    = interp1(time(Ta<50),  Ta(Ta<50),   time_r);
rh_r    = interp1(time(rh<150), rh(rh<150),  time_r);
Rs_r    = interp1(time(Rs<1000),Rs(Rs<1000), time_r);
Rl_r    = interp1(time(Rl<1000),Rl(Rl<1000), time_r);
w_spd_r = sqrt(w_u_r.^2 + w_v_r.^2);

% Note that barometric pressure, rain rate and sea surface temperature have 
% relative long periods of missing data, hence the large gaps in P, sst and
% rain are filled by interpolating the reanalysis data MERRA2 in space and 
% time to the OCSP mooring location and time.
% function to fill large gaps after 2010
[P_r,sst_r,rain_r] = merra_fill(P(pre_day:end),sst(pre_day:end),...
                     rain(pre_day:end),lat,lon,time_r); 
                
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
save([ocsp_dir,'Mooring/profs.mat'],...
    'tprof','tprof_r','PTprof','sprof','sprof_r','SAprof','dprof',...
    'dprof_r','depth_t','z_tsd','time','time_r')

clear depth_sd depth_t

% save meteorological variables
MET.w_u  = w_u_r;
MET.w_v  = w_v_r;
MET.P    = P_r;
MET.Ta   = Ta_r;
MET.rh   = rh_r;
MET.Rs   = Rs_r;
MET.Rl   = Rl_r;
MET.sst  = sst_r;
MET.sss  = sss_r;
MET.ssd  = ssd_r;
MET.rain = rain_r;
MET.zw   = z_wind;
MET.zTa  = z_Ta;
MET.zrh  = z_rh;
MET.time = time_r;

save([ocsp_dir,'Mooring/Met_2009_2019.mat'],'MET')

%% Time Conversion

% USE UTC TIME!!
date_r = string(datestr(time_r, 'yyyy-mm-dd HH:MM:SS'));

%% compute surface flux

% net trubulent heat flux = latent heat flux (positive in) + sensible 
% heat flux (positive in) + net longwave radiation (positive in)

% call COARE (3.6) bulk formula
A = coare36vn_zrf(w_spd_r,z_wind,Ta_r,z_Ta,rh_r,z_rh,P_r,sst_r,...
                  Rs_r,Rl_r,lat,NaN,rain_r,sss_r,NaN,NaN,10,10,10);
% TO-DO: relative wind velocity, wave effect

% A_as = hfbulktc(w_spd_r,z_wind,Ta_r,z_Ta,rh_r,z_rh,P_r,sst_r,sss_r,Rl_r,Rs_r,nsw);
% complex number in the output?

tau  =  A(:,2);  % wind stress [N/m^2]
hsb  = -A(:,3);  % sensible heat     flux Into the ocean [W/m^2]
hlb  = -A(:,4);  % latent   heat     flux Into the ocean [W/m^2]
hbb  = -A(:,5);  % surface  buoyancy flux Into the ocean [W/m^2]
Le   =  A(:,27); % latent   heat of  vaporization [J/kg]
Evap =  A(:,37); % evaporation rate [mm/hr]

w_cos = w_u_r./w_spd_r;
w_sin = w_v_r./w_spd_r;

% surface momentum flux in x, y direction
tau_x = tau.*w_cos;
tau_y = tau.*w_sin;

date_vec = datevec(char(date_r)); 
yd       = date2doy(time_r)-1; % yearday (starts at 1) from 'date2doy'
%-----------------------

% net short wave heat flux Into the ocean, longitude: West positive degree
nsw = swhf(yd,date_vec(:,1),(360-lon)*ones(size(time_r)),...
           lat*ones(size(time_r)),Rs_r);

% net longwave heat flux Into the ocean
nlw = lwhf(sst_r,Rl_r,Rs_r);

% net trubulent heat flux Into the ocean
hf = hlb + hsb + nlw;

% save Meteological variables
save([ocsp_dir,'/Mooring/mFRC.mat'],...
    'tau_x','tau_y','nsw','nlw','hlb','hsb','hbb','Le','Evap',...
    'sst_r','sss_r','ssd_r','rain_r','time_r')

%% Compute Wave Spectrum (frequency bin width weighted)

% observed band directional component-x and component-y
% alpha = atan2(wave_b1,wave_a1);
% xcmp  = -sin(alpha); 
% ycmp  = -cos(alpha);
xcmp  = cosd(wave_mdir); 
ycmp  = sind(wave_mdir); 

df      = repmat(wave_bw,1,length(wave_time));
wave_r1 = sqrt(wave_a1.^2 + wave_b1.^2);

% wave spectrum input for GOTM
spec_h2 = wave_spec .* wave_r1 .* df;

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

gotmdata_root = '~/Documents/GitHub/GOTM/gotmwork/data/';
basecase      = [gotmdata_root,'OCSPapa_20070608-20190616/'];

TAUname  = [basecase,'tau_file.dat'];
HFname   = [basecase,'heatflux_file.dat'];
RAINname = [basecase,'precip_file.dat'];
SSSname  = [basecase,'sss_file.dat'];
SSTname  = [basecase,'sst_file.dat'];
SWRname  = [basecase,'swr_file.dat'];
SPname   = [basecase,'sprof_file.dat'];
TPname   = [basecase,'tprof_file.dat'];
WSPname  = [basecase,'spec_file.dat']; % wave spectrum file

write_gotm_flux(TAUname, [tau_x tau_y],date_r)
write_gotm_flux(HFname,  hf,           date_r)
write_gotm_flux(RAINname,rain_r,       date_r)
write_gotm_flux(SSSname, sss_r,        date_r)
write_gotm_flux(SSTname, sst_r,        date_r)
write_gotm_flux(SWRname, nsw,          date_r)

SAprof3 = reshape(SAprof,length(z_tsd),1,length(time_r));
PTprof3 = reshape(PTprof,length(z_tsd),1,length(time_r));

write_gotm_ini(SPname,SAprof3,date_r,z_tsd)
write_gotm_ini(TPname,PTprof3,date_r,z_tsd)

waveInput  = cat(3,spec_h2,xcmp,ycmp);
waveInput3 = permute(waveInput,[1 3 2]);

write_gotm_ini(WSPname,waveInput3,wave_date,wave_freq)
gzip(WSPname)
delete(WSPname)

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
  
 