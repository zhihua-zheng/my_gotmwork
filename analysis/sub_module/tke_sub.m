%% tke_sub

% Subroutine to ananlyze TKE components from GOTM simulation output

% Zhihua Zheng, UW-APL, Sep. 22 2018

%% ----- specify model parameters -----------------------------------------
model_par.dtr0 = -0.2; % derivative of density w.r.t. temperature
model_par.dsr0 = 0.78; % derivative of density w.r.t. salinity
model_par.A1 = 0.92;
model_par.B1 = 16.6;
model_par.rho_0 = 1027; % reference density of seawater
model_par.rescale_r = 1;

%% ----- computation ------------------------------------------------------

% Note - the TKE components result is messy for OCSPapa simulation!!

% compute components of TKE
tke_comps = get_tke_comp(model_par,out,0);

% water-side friction velocity square
u_star2 = out.u_taus.^2;
u_star = out.u_taus;

% normalize TKE components by water-side friction velocity
tke_comps_n = tke_comps./(repmat(u_star2',length(zi(2:end-1)),1,3));

%% ----- plot evolution of column -----------------------------------------
plot_info.ylabel = 'depth ($$m$$)';
plot_info.clim = [0 nanmax(nanmax(nanmax(tke_comps_n)))];
%plot_info.clim = [];
plot_info.clabel = '$$\overline{w^{\prime}w^{\prime}}/u_{*}^{2}$$';
plot_info.color = 'tempo';
plot_info.plot_method = 3;
plot_info.ylim = [zi(1), 0];
plot_info.save_path = './figs/ww_norm';

% z-direction TKE
plot_time_depth(time,zi,tke_comps_n(:,:,3),plot_info)

% downwind-direction TKE
plot_info.clabel = '$$\overline{u^{\prime}u^{\prime}}/u_{*}^{2}$$';
plot_info.save_path = './figs/uu_norm';
plot_time_depth(time,zi,tke_comps_n(:,:,1),plot_info)

% crosswind-direction TKE
plot_info.clabel = '$$\overline{v^{\prime}v^{\prime}}/u_{*}^{2}$$';
plot_info.save_path = './figs/vv_norm';
plot_time_depth(time,zi,tke_comps_n(:,:,2),plot_info)

% 2*TKE 
tke_n = 2*out.tke./(repmat(u_star2',length(zi),1));

plot_info.clabel = '$$q^{2}/u_{*}^{2}$$';
plot_info.save_path = './figs/qq_norm';
%plot_info.clim = [];
plot_time_depth(time,zi,tke_n,plot_info)

%% ---- plot time averaged profiles ---------------------------------------

% TO-DO: how to choos the averaging length?

% del_t = 3; % box length ~ 3 hours
% new_t_length = floor(length(time)/del_t);
% tke_comps_n_av = nan(length(zi),new_t_length,3);
% 
% for i = 1:new_t_length
%     
%     tke_comps_n_av(:,i,:) = nanmean(tke_comps_n(:,(i-1)*del_t+1:i*del_t,:),2);
% end
% 
% plot(tke_comps_n(:,45,1),zi/mld)
% ylim([-2 0])

%% ---- averaged vertical rms velocity .vs. friction velocity -------------

ww_ml = average_ml(mld,tke_comps(:,:,3),zi,mld_smooth);
w_rms = sqrt(ww_ml);

figure('position', [0, 0, 500, 480])

scatter(u_star,w_rms,50,'s','MarkerFaceColor',...
    rgb('light turquoise'),'MarkerEdgeColor',[.5 .5 .5]);
hold on
v_lim = max(max([u_star w_rms]));

inx_good = ~isnan(w_rms);
sl = u_star(inx_good)\w_rms(inx_good);
h_fit = plot((0:.01:.15),sl*(0:.01:.15),'Color',...
    rgb('turquoise'),'LineWidth',1.5);

h_ref = refline(1.07,0);
set(h_ref,'Color',[.2 .2 .2],'LineWidth',1.5,'LineStyle','--');


xlim([0 .15])
ylim([0 .15])
xticks([0 .05 .1 .15])
yticks([0 .05 .1 .15])
% boundedline((0:23),hr_tr,hr_tr_error,...
%     'orientation','vert','alpha','transparency',0.4,'cmap',rand(1,3));
    
axis square
box on
grid on

lgd = legend([h_fit h_ref],{['slope = ',num2str(round(sl,2))],...
    'slope = 1.07'},'Location','best');
set(lgd,'Interpreter','latex','fontsize',22,'color','none')

xlabel('friction velocity $$u_*$$', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
ylabel('vertical velocity $$w_{rms}$', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
set(gca,'fontsize',20,'fontname','computer modern','gridlinestyle','--',...
    'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')

export_fig('./figs/w_ustar','-eps','-transparent','-painters')
