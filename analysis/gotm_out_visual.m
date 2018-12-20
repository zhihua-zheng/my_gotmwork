%% gotm_out_visual

% Main program to analyze output data from GOTM simulations

% Zhihua Zheng, UW-APL, Sep. 5 2018


%% General Configuration before Analysis

init_analy;

%% Heat Content

lable_as_yd = 0; % 1 to use yearday, 0 to use time string
heat_content;

%% Mixed Layer Depth (diagnosed from Ri criteria, mld_method = 2)

mixed_layer_d;

plot_info.clim = [];
plot_info.clabel = '$log_{10}(TKE)$ $m^2/s^2$';
plot_info.ylabel = 'depth (m)';
plot_info.ylim = [-300 0];
plot_info.color = 'tempo';
plot_info.plot_method = 3;
plot_info.save_path = [];
plot_time_depth(time,zi,log10(out.tke),plot_info)
hold on
plot(time,-mld,'Color',rgb('pinkish'),'LineWidth',.1)
% hold on 
% contour(log10(out.tke),[-5 -5],'LineWidth',0.01,'LineColor','k')

saveas(gcf,'./figs/mld_on_tke','fig');


%% Potential Energy (PE)

% PE per unit volume for the whole water column
g = 9.81; %[m/s^2]
pe = sum(rho.*h*g); %[J/m^2]
pe_obs = sum(rho_obs.*h*g); %[kg/(m*s^2)]

figure('position', [0, 0, 900, 300])
line(time,pe_obs./10^6,'LineWidth',.1,'Color',rgb('light teal'))
line(time,pe./10^6,'LineWidth',1,'Color',[.7 .4 .6])

plot_info.lgd = {'obs. PE',[turb_method,' PE']};
plot_info.ylabel = 'potential energy ($MJ/m^{2}$)';
plot_info.save_path = './figs/PE';
line_annotate(plot_info)

windowSize = 8; % average in one day
b = (1/windowSize)*ones(1,windowSize);
a = 1;
pe_diff = filter(b,a,pe-pe_obs);

figure('position', [0, 0, 900, 300])
line(time,pe_diff,'LineWidth',.1,'Color',rgb('ocean blue'))
line(time,zeros(size(pe)),'LineWidth',.1,'Color',rgb('royal purple'),'LineStyle','--')
plot_info.lgd = turb_method;
plot_info.ylabel = 'model PE - obs. PE ($J/m^{2}$)';
plot_info.save_path = './figs/PE_diff';
line_annotate(plot_info)

%% Currents

u = out.u;
v = out.v;
u_stokes = out.u_stokes;
v_stokes = out.v_stokes;

% find where the mixed layer depth is in variable z
ml_mask = out.z >= -repmat(mld',length(z),1);  

% u(u>1 | u<-1) = 0;
% v(v>1 | v<-1) = 0;
u_surf = u(end,:);
v_surf = v(end,:);
cur_surf = complex(u_surf,v_surf);
cur = complex(u,v);

% average current in mixed layer
cur_a = zeros(size(time));
for j = 1:length(time)
     
     tmp = cur(:,j);
     cur_a(j) = mean(tmp(ml_mask(:,j)));
     
     % fill the outlier with average value
%      tmp(abs(tmp)>=1) = cur_a(j);
%      cur(:,j) = tmp;
end

cur_a(ml_mask(end,:) == 0) = 0; % avoid NaN when mld is 0

% temporal evolution of current vector vertex
hodogram_t(time,cur_a,0)
  
%------- FFT power spectral density estimate ------------------------------

% f = (1/(nsave*dt))*(0:(n/2))/n; % frequency domain
% 
% p_u_surf = abs((fft(u_surf,n))/n);
% loglog(f,p_u_surf(1:n/2+1))
%--------------------------------------------------------------------------


%------- Welch's power spectral density estimate --------------------------

cur_a = cur_a - mean(cur_a);
 
% the next power of 2 greater than signal length - FFT transfrom length
n = 2^nextpow2(length(cur_a));

[p_cur, f] = pwelch(cur_a,[],[],n,24*3600/(nsave*dt)); % f in cycle/day

rotary_spec(f,p_cur,24/t_Coriolis,1)
%--------------------------------------------------------------------------


%--- use jLab functions ---------------------------------------------------

% psi = sleptap(length(cur_a),16); 
% [f,spp,snn] = mspec(3/24,cur_a,psi);
% plot(f/2/pi,[snn spp]),xlim([f(2) f(end)]/2/pi),xlog,ylog

%--------------------------------------------------------------------------

% plot_info.save = 0;
% plot_info.plot_method = 1;
% plot_time_depth(time,z,abs(cur),plot_info)

%% Turbulence Statistics - turbulent fluxes & TKE components

tke_sub;

%% Eddy Diffusivity

% index = find(dateVec(:,1)==2011 & dateVec(:,2)==8 & dateVec(:,3)==5);

D_e = out.D_e;
nu_m = out.nu_m;
nu_h = out.nu_h;
nu_s = out.nu_s;

figure('position', [0, 0, 400, 600])
semilogx(De_keps_pick,z,'LineWidth',.4,'Color',[.8 .7 .2])

% hold on
% semilogx(nu_h_keps_pick,zi,'LineWidth',.4,'Color',[.1 .7 .2])
% same the the turbulent diffusivity for salt

hold on
semilogx(nu_m_keps_pick,zi,'LineWidth',.4,'Color',[.9 .4 .8])
% very different than eddy diffusivity and the turbulent diffusivity for salt

hold on
semilogx(nu_s_keps_pick,zi,'LineWidth',.4,'Color',[.4 .3 .5])
% similar to eddy diffusivity

% hold on
% semilogx(De_kpp_pick,z,'LineWidth',.4,'Color',[.4 .2 .8])

  hold off
  box on
  lgd = legend('eddy diffusivity','$\nu_{m}$','$\nu_{s}$','Location','best');
  % lgd = legend('$k-\varepsilon$','KPP','Location','best');
  set(lgd,'Interpreter','latex','fontsize', 14)
  ylabel('depth (m)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  xlabel('eddy diffusivity ($$m^{2}/s$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
  setDateAxes(gca,'YLim',[-100 0],'XLim',[10^(-8) .1],...
      'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

  export_fig ('./figs/keps_diffusivity_comparison_2011','-pdf','-transparent','-painters')
  % export_fig ('./figs/eddy_diffusivity_comparison_2011','-pdf','-transparent','-painters')

%% Evolution of Temperature Profile (prediction and observation)

% model prediction
plot_info.ylabel = 'depth ($$m$$)';
plot_info.clim = [];
plot_info.clabel = 'potential temperature ($$^{\circ}C$$)';
plot_info.color = 'haline';
plot_info.plot_method = 1;
plot_info.ylim = [-150, 0];
plot_info.save_path = './figs/temp';

plot_time_depth(time,z,temp,plot_info)

% observation
plot_info.clim = [min(min(temp)) max(max(temp))];
plot_info.clabel = 'obs. potential temperature ($$^{\circ}C$$)';
plot_info.color = 'haline';
plot_info.save_path = './figs/temp_obs';

plot_time_depth(time,z,temp_obs,plot_info)

% prediction - observation
plot_info.clim = 'symmetric';
plot_info.clabel = 'temperature diffence ($$^{\circ}C$$)';
plot_info.color = 'curl';
plot_info.save_path = [];
plot_info.save_path = './figs/temp_diff';

plot_time_depth(time,z,temp-temp_obs,plot_info)


%% SST

do_sst_analysis;

%% length scale

plot_info.clim = [];
plot_info.color = 'tempo';
plot_info.save_path = './figs/length';
plot_info.clabel = 'length scale ($$m$$)';
plot_info.ylabel = 'depth (m)';
plot_info.ylim = [zi(1) 0];
plot_info.plot_method = 1;

plot_time_depth(time,zi,out.L,plot_info)

hold on 
line(time,-mld,'LineWidth',.1,'Color',[.3 .2 .1],'LineStyle','--')

set(gca(),'LooseInset', get(gca(),'TightInset')); % no blank edge
saveas(gcf, plot_info.save_path, 'epsc');

% clean the white lines in the patch
epsclean([plot_info.save_path,'.eps'],'closeGaps',true) 

