%% ideal_hurricane_sec

% Main program to analyze and compare output data from a set of GOTM 
% simulations for a transect across an idealized hurricane track

% Zhihua Zheng, UW-APL, Nov. 30 2018

%% Build data in mat form

% root_dir='~/Documents/GitLab/GOTM_dev/run/Idealized_Hurricane_Experiment/';
% ST_dir='Idealized_Hurricane_SMC_20110401-20110404/STORAGE';
% LT_dir='Idealized_Hurricane_SMCLT_20110401-20110404/STORAGE';
% 
% cd([root_dir,LT_dir])
% [outLT, ~] = load_gotm_out();
% 
% cd([root_dir,ST_dir])
% [outST, ~] = load_gotm_out();
% 
% out = [outST, outLT];
% 
% cd(root_dir)
% save('IH_all','out');
% 
% clear outLT outST

%% Load data

% GOTM output
load IH_all

% LES output
% load TC021_PROF
% load TC031_PROF

%% General variabels
time = out{1}.time;
dt = (time(2) - time(1))*24*3600; 
yd = date2doy(time)-1; % yearday output from date2doy starts as 1

z = mean(out{1}.z,2);
zi = mean(out{1}.zi,2);
h = mean(out{1}.h,2); % layer thickness
mld_smooth = 0;

%% Wind Stress

figure('position', [0, 0, 950, 450])
set(gcf,'defaultAxesColorOrder',[rgb('black'); rgb('red')]);

plot_info.grid_on = 0;
plot_info.xlabel = [];
plot_info.ylabel = [];
plot_info.x_lim = [90 93];
plot_info.y_lim = [-0.45 0.45];
plot_info.lgd = [];
plot_info.lgd_pos = [];
plot_info.save_path = [];

subplot(2,1,1)
yyaxis right
quiver(yd',zeros(size(yd')),out{46,2}.u_stokes(end,:),...
    out{46,2}.v_stokes(end,:),1.5,'Marker','.','Linewidth',1.5,...
    'Color',rgb('red'),'ShowArrowHead','off')
line_annotate(plot_info)
text(92.1,-0.2,'$u^{St}_0$ at RHS maximum [m/s]','Interpreter',...
    'latex','fontsize',15,'Color',rgb('red'));

yyaxis left
quiver(yd,zeros(size(yd)),out{46,2}.tx,out{46,2}.ty,1.5,'Marker','.',...
    'Color','k','ShowArrowHead','off','Linewidth',3);
text(90.1,-0.2,'$u^2_*$ at RHS maximum [$m^2/s^2$]','Interpreter','latex','fontsize',15);
line_annotate(plot_info)

subplot(2,1,2)
yyaxis right
quiver(yd',zeros(size(yd')),out{56,2}.u_stokes(end,:),...
    out{56,2}.v_stokes(end,:),1.5,'Marker','.','Linewidth',1.5,...
    'Color',rgb('red'),'ShowArrowHead','off')
line_annotate(plot_info)
text(92.1,0.2,'$u^{St}_0$ at LHS maximum [$m/s$]','Interpreter',...
    'latex','fontsize',15,'Color',rgb('red'));

yyaxis left
quiver(yd,zeros(size(yd)),out{56,2}.tx,out{56,2}.ty,1.5,'Marker','.',...
    'Color','k','ShowArrowHead','off','Linewidth',3);
text(90.1,0.2,'$u^2_*$ at LHS maximum [$m^2/s^2$]','Interpreter','latex','fontsize',15);
line_annotate(plot_info)

export_fig('./figs/S_wind_wave','-eps','-transparent','-painters');

%% Langmuir number

figure('position', [0, 0, 900, 450])
plot_info.grid_on = 1;
plot_info.xlabel = [];
plot_info.ylabel = [];
plot_info.x_lim = [90 93];
plot_info.y_lim = [0 2.8];
plot_info.lgd = {'$La_t$','$La_{SL}$'};
plot_info.lgd_pos = [];
plot_info.save_path = [];
plot_inx = 1;
loc_text = {'Langmuir number at RHS maximum',...
    'Langmuir number at LHS maximum'};

for j = [46 56]
    subplot(2,1,plot_inx)
    La_t = get_La(u_star{j,2},out{j,2}.u_stokes,...
        out{j,2}.v_stokes,1,[],[],[]);
    La_sl = get_La(u_star{j,2},out{j,2}.u_stokes,...
        out{j,2}.v_stokes,2,out{j,2}.mld_surf,z,h);

    line(yd,La_t,'LineWidth',3,'Color',rgb('azure'));
    line(yd,La_sl,'LineWidth',3,'Color',rgb('aqua'));
    line(yd,ones(size(yd))*0.3,'LineWidth',2,'Color',[.3 .3 .3],'LineStyle',':')
    line_annotate(plot_info)
    text(91,2,loc_text{plot_inx},'Interpreter','latex','fontsize',15.5);

    plot_inx = plot_inx+1;
end
export_fig('./figs/La_number','-eps','-transparent','-painters');

%% ----- TKE components anaylysis -----------------------------------------

model_par.dtr0 = -0.2; % derivative of density w.r.t. temperature
model_par.dsr0 = 0.78; % derivative of density w.r.t. salinity
model_par.A1 = 0.92;
model_par.B1 = 16.6;
model_par.rho_0 = 1027; % reference density of seawater
model_par.rescale_r = 1;

% compute components of TKE
tke_comps = cellfun(@(A) get_tke_comp(model_par,A,0),out,'UniformOutput',false);

% water-side friction velocity
u_star = cellfun(@(A) A.u_taus,out,'UniformOutput',false);
u_star2 = cellfun(@(A) A.^2,u_star,'UniformOutput',false);

% normalize TKE components by water-side friction velocity
tke_comps_n = cellfun(@(A,B) A./(repmat(B',length(zi(2:end-1)),1,3)),...
    tke_comps,u_star2,'UniformOutput',false);

% average vertical TKE within the mixed layer 
ww_ml = cellfun(@(A,B) average_ml(B.mld_surf,A(:,:,3),zi,mld_smooth),...
    tke_comps,out,'UniformOutput',false);
w_rms = cellfun(@sqrt,ww_ml,'UniformOutput',false);

%% Scaling of vertical velocity

v_lim = 0.11;

figure('position', [0, 0, 800, 780])
scat_info.grid_on = 1;
scat_info.ticks = (0:0.02:v_lim);
scat_info.xlabel = [];
scat_info.ylabel = [];
scat_info.x_lim = [0 v_lim];
scat_info.y_lim = [0 v_lim];
scat_info.lgd = [];    
scat_info.lgd_obj = [];
scat_info.lgd_pos = [];
scat_info.save_path = [];
% scat_info.save_path = './figs/w_ustar';
plot_inx = 1;
s_mark = ['s','o'];
loc_text = {'RHS maximum','LHS maximum';'RHS maximum','LHS maximum'};

for j = [46 56]
    
    % maxmium wind time
    [~,m] = max(u_star{j,i});
    
    % start of high wind forcing
    s = find(u_star{j,i}>=0.015,1,'first');
    
    % end of high wind forcing
    e = find(u_star{j,i}>=0.015,1,'last');
    
    for i = 1:2 
        
    subplot(2,2,plot_inx)
    s1 = scatter(u_star{j,i}(s:m-1),w_rms{j,i}(s:m-1),200,s_mark(i),...
        'MarkerFaceColor',rgb('neon blue'),'MarkerEdgeColor',[.5 .5 .5]);
    hold on 
    s2 = scatter(u_star{j,i}(m:e),w_rms{j,i}(m:e),200,s_mark(i),...
        'MarkerFaceColor',rgb('neon red'),'MarkerEdgeColor',[.5 .5 .5]);
    hold on
    
    text(.01,.09,loc_text{plot_inx},'Interpreter','latex','fontsize',18)
    text(.07,.02,out{j,i}.turb_method,'Interpreter','latex','fontsize',18)
    sl = u_star{j,i}(s:m-1)\w_rms{j,i}(s:m-1);
    plot((0:.01:v_lim),sl*(0:.01:v_lim),'Color',...
        rgb('neon blue'),'LineWidth',2);

    sl = u_star{j,i}(m:e)\w_rms{j,i}(m:e);
    plot((0:.01:v_lim),sl*(0:.01:v_lim),'Color',...
        rgb('neon red'),'LineWidth',2);

    h_ref = refline(1.07,0);
    set(h_ref,'Color',[.2 .2 .2],'LineWidth',1.5,'LineStyle','--');

    % boundedline((0:23),hr_tr,hr_tr_error,...
    %    'orientation','vert','alpha','transparency',0.4,'cmap',rand(1,3));
    
    scat_annotate(scat_info)
    plot_inx = plot_inx + 1;
    end
end 

[~,h1] = suplabel('vertical velocity $w_{rms}$','y');
[~,h2] = suplabel('friction velocity $u_*$','x');
set([h1 h2],'FontSize',20,'FontName','computer modern','Interpreter','latex');

export_fig('./figs/w_ustar','-eps','-transparent','-painters')

% {'before maximum','after maximum'}
%% 

