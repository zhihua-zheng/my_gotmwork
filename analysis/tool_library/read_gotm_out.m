function  out = read_gotm_out(fname,order)

% read_gotm_out
%==========================================================================
%
% USAGE:
%  out = read_gotm_out(fname,order)
%
% DESCRIPTION:
%  Read the output variables from GOTM simulation netCDF file (fname)
%
% INPUT:
%
%  fname - name of output file
%  order - order of turbulence closure scheme
%
% OUTPUT:
%
%  out - the struct contains all variables in output file
%
% AUTHOR:
%  Sep.  2 2018. Zhihua Zheng                           [ zhihua@uw.edu ]
%  Oct. 20 2018. Neglect first value in time dimesion.


%% general_info
out.time = double(ncread(fname,'time',2,Inf));
% seconds since the initialized time in the simulation

out.lon = double(ncread(fname,'lon')); % longitude [degrees_east]
out.lat = double(ncread(fname,'lat'));% latitude [degrees_north]

%% column_structure
out.ga = double(squeeze(squeeze(ncread(fname,'ga',[1 1 1 2],[Inf Inf Inf Inf])))); % coordinate scaling
out.z = double(squeeze(squeeze(ncread(fname,'z',[1 1 1 2],[Inf Inf Inf Inf])))); % depth coordinate[m]
out.zi = double(squeeze(squeeze(ncread(fname,'zi',[1 1 1 2],[Inf Inf Inf Inf])))); % interface depth [m]
out.h = double(squeeze(squeeze(ncread(fname,'h',[1 1 1 2],[Inf Inf Inf Inf])))); % layer thickness [m]

%% temperature_and_salinity_and_density
out.temp = double(squeeze(squeeze(ncread(fname,'temp',[1 1 1 2],[Inf Inf Inf Inf])))); % potential temperature [Celsius]
out.temp_obs = double(squeeze(squeeze(ncread(fname,'temp_obs',[1 1 1 2],[Inf Inf Inf Inf])))); % observed potential temperature [Celsius]
out.salt = double(squeeze(squeeze(ncread(fname,'salt',[1 1 1 2],[Inf Inf Inf Inf])))); % absolute salinity [g/kg]
out.salt_obs = double(squeeze(squeeze(ncread(fname,'salt_obs',[1 1 1 2],[Inf Inf Inf Inf])))); % observed practical salinity [PSU]
out.rho = double(squeeze(squeeze(ncread(fname,'rho',[1 1 1 2],[Inf Inf Inf Inf])))); % potential density [kg/m^3]

%% surface
% out.KPP_OSBL = double(squeeze(squeeze(ncread(fname,'KPP_OSBL')))); % KPP boundary layer depth
out.mld_surf = double(squeeze(ncread(fname,'mld_surf',[1 1 2],[Inf Inf Inf]))); % surface mixed layer depth [m]
out.zeta = double(squeeze(ncread(fname,'zeta',[1 1 2],[Inf Inf Inf]))); % sea surface elevation [m]
out.sss = double(squeeze(ncread(fname,'sss',[1 1 2],[Inf Inf Inf]))); % sea surface salinity [PSU]
out.sst_obs = double(squeeze(ncread(fname,'sst_obs',[1 1 2],[Inf Inf Inf]))); % observed sea surface temperature [Celsius]
out.sst = double(squeeze(ncread(fname,'sst',[1 1 2],[Inf Inf Inf]))); % sea surface temperature [Celsius]
out.tx = double(squeeze(ncread(fname,'tx',[1 1 2],[Inf Inf Inf]))); % wind stress (tau_x/rho) [m^2/s^2]
out.ty = double(squeeze(ncread(fname,'ty',[1 1 2],[Inf Inf Inf]))); % wind stress (tau_y/rho) [m^2/s^2]
out.u_taus = double(squeeze(ncread(fname,'u_taus',[1 1 2],[Inf Inf Inf]))); % water-side surface friction velocity [m/s]
out.u10 = double(squeeze(ncread(fname,'u10',[1 1 2],[Inf Inf Inf]))); % 10m wind x-velocity [m/s]
out.v10 = double(squeeze(ncread(fname,'v10',[1 1 2],[Inf Inf Inf]))); % 10m wind y-velocity [m/s]
out.v0_stokes = double(squeeze(squeeze(ncread(fname,'v0_stokes',[1 1 2],[Inf Inf Inf])))); % surface Stokes drift y-component
out.u0_stokes = double(squeeze(squeeze(ncread(fname,'u0_stokes',[1 1 2],[Inf Inf Inf])))); % surface Stokes drift x-component

out.int_precip = double(squeeze(ncread(fname,'int_precip',[1 1 2],[Inf Inf Inf]))); % integrated precipitation [m/s]
out.int_evap = double(squeeze(ncread(fname,'int_evap',[1 1 2],[Inf Inf Inf]))); % integrated evaporation [m/s]
out.int_swr = double(squeeze(ncread(fname,'int_swr',[1 1 2],[Inf Inf Inf]))); % integrated short wave radiation [J/m^2]
out.int_heat = double(squeeze(ncread(fname,'int_heat',[1 1 2],[Inf Inf Inf]))); % integrated surface heat fluxes [J/m^2]
out.int_total = double(squeeze(ncread(fname,'int_total',[1 1 2],[Inf Inf Inf]))); % integrated total surface heat exchange [J/m^2]

out.rhoa = double(squeeze(ncread(fname,'rhoa',[1 1 2],[Inf Inf Inf]))); % air density [kg/m^3]
out.cloud = double(squeeze(ncread(fname,'cloud',[1 1 2],[Inf Inf Inf]))); % cloud cover
out.albedo = double(squeeze(ncread(fname,'albedo',[1 1 2],[Inf Inf Inf]))); % albedo
out.precip = double(squeeze(ncread(fname,'precip',[1 1 2],[Inf Inf Inf]))); % precipitation [m/s]
out.evap = double(squeeze(ncread(fname,'evap',[1 1 2],[Inf Inf Inf]))); % evaporation [m/s]
out.airt = double(squeeze(ncread(fname,'airt',[1 1 2],[Inf Inf Inf]))); % 2m air temperature [Celsius]
out.airp = double(squeeze(ncread(fname,'airp',[1 1 2],[Inf Inf Inf]))); % air pressure [Pa]
out.rh = double(squeeze(ncread(fname,'hum',[1 1 2],[Inf Inf Inf]))); % relative humidity [%]
out.es = double(squeeze(ncread(fname,'es',[1 1 2],[Inf Inf Inf]))); % saturation water vapor pressure [Pa]
out.ea = double(squeeze(ncread(fname,'ea',[1 1 2],[Inf Inf Inf]))); % actual water vapor pressure [Pa]
out.qs = double(squeeze(ncread(fname,'qs',[1 1 2],[Inf Inf Inf]))); % saturation specific humidity
out.qa = double(squeeze(ncread(fname,'qa',[1 1 2],[Inf Inf Inf]))); % specific humidity

% heat_fluxes
out.heat = double(squeeze(ncread(fname,'heat',[1 1 2],[Inf Inf Inf]))); % net surface heat flux [W/m^2]
out.qe = double(squeeze(ncread(fname,'qe',[1 1 2],[Inf Inf Inf]))); % sensible heat flux [W/m^2]
out.qh = double(squeeze(ncread(fname,'qh',[1 1 2],[Inf Inf Inf]))); % latent heat flux [W/m^2]
out.qb = double(squeeze(ncread(fname,'qb',[1 1 2],[Inf Inf Inf]))); % long-wave back radiation [W/m^2]
out.I_0 = double(squeeze(ncread(fname,'I_0',[1 1 2],[Inf Inf Inf]))); % incoming short wave radiation [W/m^2]

%% bottom
out.mld_bott = double(squeeze(ncread(fname,'mld_bott',[1 1 2],[Inf Inf Inf]))); % bottom mixed layer depth [m]
out.taub = double(squeeze(ncread(fname,'taub',[1 1 2],[Inf Inf Inf]))); % bottom stress [Pa]
out.u_taub = double(squeeze(ncread(fname,'u_taub',[1 1 2],[Inf Inf Inf]))); % bottom friction velocity [m/s]

%% velocities
out.u = double(squeeze(squeeze(ncread(fname,'u',[1 1 1 2],[Inf Inf Inf Inf])))); % x-mean flow [m/s]
out.v = double(squeeze(squeeze(ncread(fname,'v',[1 1 1 2],[Inf Inf Inf Inf])))); % y-mean flow [m/s]
out.u_obs = double(squeeze(squeeze(ncread(fname,'u_obs',[1 1 1 2],[Inf Inf Inf Inf])))); % observed x-velocity [m/s]
out.v_obs = double(squeeze(squeeze(ncread(fname,'v_obs',[1 1 1 2],[Inf Inf Inf Inf])))); % observed y-velocity [m/s]
out.v_stokes = double(squeeze(squeeze(ncread(fname,'v_stokes',[1 1 1 2],[Inf Inf Inf Inf])))); % Stokes drift y-component
out.u_stokes = double(squeeze(squeeze(ncread(fname,'u_stokes',[1 1 1 2],[Inf Inf Inf Inf])))); % Stokes drift x-component

%% turbulence
% out.taux = double(squeeze(squeeze(ncread(fname,'taux')))); % turbulent flux of x-momentum [m^2/s^2]
% out.tauy = double(squeeze(squeeze(ncread(fname,'tauy')))); % turbulent flux of y-momentum [m^2/s^2]

out.Px = double(squeeze(squeeze(ncread(fname,'xP',[1 1 1 2],[Inf Inf Inf Inf])))); % extra turbulence production [m^2/s^3]
out.tke = double(squeeze(squeeze(ncread(fname,'tke',[1 1 1 2],[Inf Inf Inf Inf])))); % turbulent kinetic energy [m^2/s^2]
out.Rig = double(squeeze(squeeze(ncread(fname,'Rig',[1 1 1 2],[Inf Inf Inf Inf])))); % gradient Richardson number

out.D_e = double(squeeze(squeeze(ncread(fname,'avh',[1 1 1 2],[Inf Inf Inf Inf])))); % eddy diffusivity [m^2/s]
out.nu_m = double(squeeze(squeeze(ncread(fname,'num',[1 1 1 2],[Inf Inf Inf Inf])))); % turbulent diffusivity of momentum - down the Eulerian shear [m^2/s]
out.nu_cl = double(squeeze(squeeze(ncread(fname,'nucl',[1 1 1 2],[Inf Inf Inf Inf])))); % Craik-Leibovich turbulent diffusivity of momentum - down the Stokes shear [m^2/s]
out.nu_h = double(squeeze(squeeze(ncread(fname,'nuh',[1 1 1 2],[Inf Inf Inf Inf])))); % turbulent diffusivity of heat [m^2/s]
out.nu_s = double(squeeze(squeeze(ncread(fname,'nus',[1 1 1 2],[Inf Inf Inf Inf])))); % turbulent diffusivity of salt [m^2/s]

out.gamu = double(squeeze(squeeze(ncread(fname,'gamu',[1 1 1 2],[Inf Inf Inf Inf])))); % non-local flux of x-momentum [m^2/s^2]
out.gamv = double(squeeze(squeeze(ncread(fname,'gamv',[1 1 1 2],[Inf Inf Inf Inf])))); % non-local flux of y-momentum [m^2/s^2]
out.gamh = double(squeeze(squeeze(ncread(fname,'gamh',[1 1 1 2],[Inf Inf Inf Inf])))); % non-local heat flux [K*m/s]
out.gams = double(squeeze(squeeze(ncread(fname,'gams',[1 1 1 2],[Inf Inf Inf Inf])))); % non-local salinity flux [(g/kg)*(m/s)]
out.gam = double(squeeze(squeeze(ncread(fname,'gam',[1 1 1 2],[Inf Inf Inf Inf])))); % non-dimensional non-local buoyancy flux

if order == 2

 out.eps = double(squeeze(squeeze(ncread(fname,'eps',[1 1 1 2],[Inf Inf Inf Inf])))); % energy dissipation rate [m^2/s^3]

 out.Rif = double(squeeze(squeeze(ncread(fname,'xRf',[1 1 1 2],[Inf Inf Inf Inf])))); % flux Richardson number

 out.gamb = double(squeeze(squeeze(ncread(fname,'gamb',[1 1 1 2],[Inf Inf Inf Inf])))); % non-local buoyancy flux [m^2/s^3]

 out.L = double(squeeze(squeeze(ncread(fname,'L',[1 1 1 2],[Inf Inf Inf Inf])))); % turbulence length scale [m]
 out.r = double(squeeze(squeeze(ncread(fname,'r',[1 1 1 2],[Inf Inf Inf Inf])))); % turbulent time scale ratio

 out.ab = double(squeeze(squeeze(ncread(fname,'an',[1 1 1 2],[Inf Inf Inf Inf])))); % non-dimensional buoyancy time scale
 out.as = double(squeeze(squeeze(ncread(fname,'as',[1 1 1 2],[Inf Inf Inf Inf])))); % non-dimensional shear time scale
 out.abk = double(squeeze(squeeze(ncread(fname,'at',[1 1 1 2],[Inf Inf Inf Inf])))); % non-dimensional (half) buoyancy variance

 out.cmue1 = double(squeeze(squeeze(ncread(fname,'cmue1',[1 1 1 2],[Inf Inf Inf Inf])))); % stability function for momentum diffusivity
 out.cmue2 = double(squeeze(squeeze(ncread(fname,'cmue2',[1 1 1 2],[Inf Inf Inf Inf])))); % stability function for scalar diffusivity

% buoyancy
 out.buoy_P = double(squeeze(squeeze(ncread(fname,'Pb',[1 1 1 2],[Inf Inf Inf Inf])))); % production of (half) buoyancy variance [m^2/s^5]
 out.Pb = double(squeeze(squeeze(ncread(fname,'G',[1 1 1 2],[Inf Inf Inf Inf])))); % buoyancy production [m^2/s^3]
 out.epsb = double(squeeze(squeeze(ncread(fname,'epsb',[1 1 1 2],[Inf Inf Inf Inf])))); % destruction of (half) buoyancy variance [m^2/s^5]
 out.kb = double(squeeze(squeeze(ncread(fname,'kb',[1 1 1 2],[Inf Inf Inf Inf])))); % (half) buoyancy variance [m^2/s^4]
end

out.buoy = double(squeeze(squeeze(ncread(fname,'buoy',[1 1 1 2],[Inf Inf Inf Inf])))); % buoyancy [m/s^2]
out.NNT = double(squeeze(squeeze(ncread(fname,'NNT',[1 1 1 2],[Inf Inf Inf Inf])))); % contribution of T-gradient to buoyancy frequency squared [1/s^2]
out.NNS = double(squeeze(squeeze(ncread(fname,'NNS',[1 1 1 2],[Inf Inf Inf Inf])))); % contribution of S-gradient to buoyancy frequency squared [1/s^2]
out.NN = double(squeeze(squeeze(ncread(fname,'NN',[1 1 1 2],[Inf Inf Inf Inf])))); % buoyancy frequency squared [1/s^2]

% shear
if order == 2

 out.uu = double(squeeze(squeeze(ncread(fname,'uu',[1 1 1 2],[Inf Inf Inf Inf])))); % variance of turbulent x-velocity (fluctuation) [m^2/s^2]
 out.vv = double(squeeze(squeeze(ncread(fname,'vv',[1 1 1 2],[Inf Inf Inf Inf])))); % variance of turbulent y-velocity (fluctuation) [m^2/s^2]
 out.ww = double(squeeze(squeeze(ncread(fname,'vv',[1 1 1 2],[Inf Inf Inf Inf])))); % variance of turbulent z-velocity (fluctuation) [m^2/s^2]
 out.Ps = double(squeeze(squeeze(ncread(fname,'P',[1 1 1 2],[Inf Inf Inf Inf])))); % shear production [m^2/s^3]
 out.PS = double(squeeze(squeeze(ncread(fname,'PS',[1 1 1 2],[Inf Inf Inf Inf])))); % Stokes production [m^2/s^3]

end

out.SS = double(squeeze(squeeze(ncread(fname,'SS',[1 1 1 2],[Inf Inf Inf Inf])))); % shear frequency squared [1/s^2]
out.drag = double(squeeze(squeeze(ncread(fname,'drag',[1 1 1 2],[Inf Inf Inf Inf])))); % drag coefficient in water column
out.fric = double(squeeze(squeeze(ncread(fname,'fric',[1 1 1 2],[Inf Inf Inf Inf])))); % extra friction coefficient in water column


% Langmuir
out.theta_WW = double(squeeze(squeeze(ncread(fname,'theta_WW',[1 1 2],[Inf Inf Inf])))); % angle between wind and waves
out.theta_WL = double(squeeze(squeeze(ncread(fname,'theta_WL',[1 1 2],[Inf Inf Inf])))); % angle between wind and Langmuir cells
out.La_SLP2 = double(squeeze(squeeze(ncread(fname,'La_SLP2',[1 1 2],[Inf Inf Inf])))); % surface layer averaged, projected Langmuir number (RWHGK16)
out.La_SLP1 = double(squeeze(squeeze(ncread(fname,'La_SLP1',[1 1 2],[Inf Inf Inf])))); % surface layer averaged, projected Langmuir number (VFSHH12)
out.La_SL = double(squeeze(squeeze(ncread(fname,'La_SL',[1 1 2],[Inf Inf Inf])))); % surface layer averaged Langmuir number
out.La_Turb = double(squeeze(squeeze(ncread(fname,'La_Turb',[1 1 2],[Inf Inf Inf])))); % turbulent Langmuir number
out.delta = double(squeeze(squeeze(ncread(fname,'delta',[1 1 2],[Inf Inf Inf])))); % Stokes drift penetration depth

%% column_integrals
out.Eturb = double(squeeze(ncread(fname,'Eturb',[1 1 2],[Inf Inf Inf]))); % turbulent kinetic energy [J]
out.Epot = double(squeeze(ncread(fname,'Epot',[1 1 2],[Inf Inf Inf]))); % potential energy [J]
out.Ekin = double(squeeze(ncread(fname,'Ekin',[1 1 2],[Inf Inf Inf]))); % kinetic energy [J]

%% light
out.bioshade = double(squeeze(squeeze(ncread(fname,'bioshade',[1 1 1 2],[Inf Inf Inf Inf]))));
% fraction of visible light that is not shaded by overlying biogeochemistry

out.swr = double(squeeze(squeeze(ncread(fname,'rad',[1 1 1 2],[Inf Inf Inf Inf])))); % short-wave radiation [W/m^2]

%% -- Get date strings

  % get the reference time
  t_ref = ncreadatt(fname,'time','units'); % get the attribute 'units' for 'time'
  t_ref = t_ref(15:end); % truncate to get the time string
  t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS'); % datenumber for initialized time of simulation

out.time = t_ref + out.time./3600/24; % seconds since t_ref
out.date = string(datestr(out.time, 'yyyy/mm/dd HH:MM:SS'));

end
