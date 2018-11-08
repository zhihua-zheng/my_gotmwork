%% multi_year_ocsp

% Main program to analyze and compare output data from GOTM simulations for
% Ocean Climate Station Papa in multiple years

% Zhihua Zheng, UW-APL, Oct. 18 2018

%% read outputs
clear
main_dir = '~/Documents/GitLab/GOTM_dev/run/OCSPapa_multi_year';
cd(main_dir)
dir_str = genpath(main_dir);

% regular expression for any substring doesn't contain path seperator
expression = ['[^',pathsep,']*'];
sep_idx = regexp(dir_str,pathsep); % index of path seperator

% truncate the string to ignore the main folder
sub_folders = regexp(dir_str(sep_idx(1):end),expression,'match');

% eliminate figs folder, delete the entry satisfying condition
sub_folders(endsWith(sub_folders,'/figs')) = [];

turb_method = cell(size(sub_folders));
start_date = cell(7,size(sub_folders,2));
out_all = cell(7,size(sub_folders,2));

% loop through subfolders
for j = 1:length(sub_folders)
    
    path = sub_folders{j};
    tmp = strsplit(path,'/');
    tmp = strsplit(tmp{end},'_');
    turb_method(j) = tmp(2); % get turbulence closure method name
    
    cd(path)
    dinfo = dir(fullfile('./*.nc'));
    fname = fullfile('./',{dinfo.name});
    
    % read outputs
    for i = 1:length(fname)
        out_all{i,j} = read_gotm_out(fname{i},2);
        % get start date
        tmp = strsplit(fname{i},{'_','.','-'});
        start_date(i,j) = strsplit(tmp{end-2});
    end
end

clear tmp spe_inx expression fname dinfo
cd(main_dir)

%% read obs. density profiles 
cd ~/Documents/Study/Graduate_research/data_raw/OCS_P/Mooring/original_data_2007/profiles/
rho_obs = ncread('d50n145w_hr.cdf','STH_71');
rho_obs = double(squeeze(rho_obs));
rho_obs(rho_obs>100) = NaN;

rho_obs_time = double(ncread('d50n145w_hr.cdf','time'));
t_ref = datenum('2007-06-08 00:00:00', 'yyyy-mm-dd HH:MM:SS');
rho_obs_time = rho_obs_time/24 + t_ref;
rho_obs_date = string(datestr(rho_obs_time, 'yyyy/mm/dd HH:MM:SS'));
rho_obs_z = -double(ncread('d50n145w_hr.cdf','depth'));

cd(main_dir)

%% save data

save('ocsp_years','out_all','start_date','turb_method','rho_obs',...
    'rho_obs_z','rho_obs_time','-v7.3');

%% read variables
all_sst = cell(size(out_all));
all_sst_obs = cell(size(out_all));
all_time = cell(size(out_all));
all_date = cell(size(out_all));

% compute mld based on density profile (0.1*sigma)
all_rho = cell(size(out_all));

for i = 1:size(out_all,1)
    for j = 1:size(out_all,2)
        
         all_sst{i,j} = out_all{i,j}.temp(128,:)';
         all_sst_obs{i,j} = out_all{i,j}.sst_obs; 
         all_time{i,j} = out_all{i,j}.time;
         all_date{i,j} = datevec(char(out_all{i,j}.date));
         all_rho{i,j} = out_all{i,j}.rho;
    end
end

all_rho_obs = cell(size(out_all));
z_model = cell(size(out_all));
t_model = cell(size(out_all));

for i = 1:size(out_all,1)
    for j = 1:size(out_all,2)
        
        first_point = all_time{i,j}(1);
        last_point = all_time{i,j}(end);
        first_inx = find(rho_obs_time == first_point);
        last_inx = find(rho_obs_time == last_point);

        tmp = rho_obs(:,first_inx:last_inx);
        
        z_model{i,j} = out_all{i,j}.z(22:end-3,:); % 1-200m
        t_model{i,j} = repmat(all_time{i,j}',size(z_model{i,j},1),1);
        [T,Z] = meshgrid(rho_obs_time(first_inx:last_inx),rho_obs_z);

        % interpolate to the model grid
        F = scatteredInterpolant(T(~isnan(tmp)),Z(~isnan(tmp)),...
            tmp(~isnan(tmp)),'linear'); % interpolation object
        all_rho_obs{i,j} = F(t_model{i,j},z_model{i,j});
    end
end

clear tmp first_point first_inx last_point last_inx F T Z

%% Mixed layer depth
all_mld = cell(size(out_all));
all_mld_obs = cell(size(out_all));

for i = 1:size(out_all,1)
    for j = 1:size(out_all,2)
    
        % 48-hr running mean
        all_mld{i,j} = movmean(get_mld(all_rho{i,j},mean(out_all{i,j}.z,2)),16);
        all_mld_obs{i,j} = movmean(get_mld(1000+all_rho_obs{i,j},...
            mean(z_model{i,j},2)),16);
    end
end

%% SST & MLD statistics

% find bounding values in SST & MLD
sst_min = cellfun(@min,all_sst);
sst_max = cellfun(@max,all_sst);
sst_obs_min = cellfun(@min,all_sst_obs);
sst_obs_max = cellfun(@max,all_sst_obs);

mld_min = cellfun(@min,all_mld);
mld_max = cellfun(@max,all_mld);
mld_obs_min = cellfun(@min,all_mld_obs);
mld_obs_max = cellfun(@max,all_mld_obs);

num = cellfun(@numel,all_sst,'UniformOutput',false); % number of points

% squared model-obs. error for SST & MLD
sst_se = cellfun(@(A,B) (A-B).^2,all_sst,all_sst_obs,'UniformOutput',false);
mld_se = cellfun(@(A,B) (A-B).^2,all_mld,all_mld_obs,'UniformOutput',false); 

% RMSE for each simulation 
sst_rmse = cellfun(@(A,B) sqrt(sum(A)/B),sst_se,num);
mld_rmse = cellfun(@(A,B) sqrt(sum(A)/B),mld_se,num);

% composite mean of RMSE
all_sst_rmse = mean(sst_rmse); 
all_mld_rmse = mean(mld_rmse(1:6,:)); % avoid using the simulation for last year

% ----- seasonal parameters -----------------------------------------------

season_list = [6  7  8;
               9  10 11;
               12 1  2;
               3  4  5];
season_str = {'J-J-A','S-O-N','D-J-F','M-A-M'};

all_month = cellfun(@(A) A(:,2),all_date,'UniformOutput',false);

% ----- seasonal parameters -----------------------------------------------


% sst_diff_detr = cell(size(out_all)); % detrended all_sst_diff
% [s, int] = trend(all_sst_diff{i},all_time{i});
% sst_diff_detr{i} = all_sst_diff{i} - (int + all_time{i}.*s);      


% for i = 1:7
%     for j=1:2
%     
%     figure
%     contourf(all_rho_obs{i,j},'LineStyle','none');colorbar
%     end
% end

%% SST error spectral ananlysis

for i = 1:length(out_all)
    
    plot(all_sst_diff{i},'Color',rand(1,3))
    hold on   
end

for i = 1:length(out_all)
    
    plot(all_mld_diff{14},'Color',rand(1,3))
    hold on   
end

figure('position', [0, 0, 600, 200])

% the next power of 2 greater than signal length - FFT transfrom length
n = 2^nextpow2(length(sst_diff_detr{i}));

[p_sst, f] = pwelch(sst_diff_detr{i},[],[],n,24*3600/(1800*6)); % f in cycle/day

% counter-clockwise, negative
loglog(f(1:n/2),p_sst(1:n/2),'Color',[.1 .6 .7],'LineStyle','-','LineWidth',.3) 

% plot(f,p_sst)
% set(gca, 'TickLabelInterpreter', 'latex','yscale','log')

%% --------- overall SST scatter plot -------------------------------------
do_sst_scatter

%% --------- seasonal SST scatter plot ----------------------------------------

do_sst_scatter_season

%% --------- overall MLD scatter plot -------------------------------------
do_mld_scatter

%% --------- seasonal MLD scatter plot ------------------------------------

do_mld_scatter_season

%% vertical velocity

% TO-DO: 1.Unrealistic large value for ww at i=2,3 simulation, but TkE
%          value is normal. The bug is due to the computation of TKE
%          components.
%       

model_par.dtr0 = -0.2; % derivative of density w.r.t. temperature
model_par.dsr0 = 0.78; % derivative of density w.r.t. salinity
model_par.A1 = 0.92;
model_par.B1 = 16.6;
model_par.rho_0 = 1027; % reference density of seawater
model_par.rescale_r = 1;
kapa = 0.4; % Von Karman constant;

figure('position', [0, 0, 500, 480])

all_u_star = cellfun(@(A) A.u_taus,out_all,'UniformOutput',false);
all_heat = cellfun(@(A) A.heat,out_all,'UniformOutput',false);
all_swr = cellfun(@(A) A.I_0,out_all,'UniformOutput',false);
all_hf = cellfun(@(A,B) A+B,all_heat,all_swr,'UniformOutput',false); % positive into the ocean
all_MOL = cellfun(@(A,B) A.^3./(kapa*B),all_u_star,all_hf,'UniformOutput',false); % Monin?Obukhov length scale
all_stable_inx = cellfun(@(A) A>0,all_hf,'UniformOutput',false);

all_tke_comps = cellfun(@(A) get_tke_comp(model_par,A,0),out_all,'UniformOutput',false);
zi = mean(out_all{1}.zi,2);

% average the ww in the mixed layer
all_ww_ml = cellfun(@(A,B) average_ml(B,A(:,:,3),zi,1),all_tke_comps,...
    all_mld,'UniformOutput',false);

% apply stable criteria
all_ww_ml_st = cellfun(@(A,B) A(B),all_ww_ml,all_stable_inx,'UniformOutput',false);
all_u_star_st = cellfun(@(A,B) A(B),all_u_star,all_stable_inx,'UniformOutput',false);

ww_ml(:,2) = vertcat(all_ww_ml_st{:,2});
ww_ml(:,1) = vertcat(all_ww_ml_st{:,1});
ww_ml(ww_ml<0) = NaN;

u_star(:,2) = vertcat(all_u_star_st{:,2});
u_star(:,1) = vertcat(all_u_star_st{:,1});

w_rms = sqrt(ww_ml);
v_lim = max(max([w_rms u_star]));

% for turb_method{2}
dummy = randn(size(u_star(:,2)));        
s2 = scatter3(u_star(:,2),w_rms(:,2),dummy,12,rgb('coral'),'filled');
s2.MarkerFaceAlpha = 0.9;
hold on

% for turb_method{1}
dummy = randn(size(u_star(:,1)));        
s1 = scatter3(u_star(:,1),w_rms(:,1),dummy,12,rgb('aqua'),'filled');
s1.MarkerFaceAlpha = 0.9;
hold on

% Change viewpoint 
view(2)

hold off 
axis square
box on
grid on

xlim([0 1.01*v_lim])
ylim([0 1.01*v_lim])

h_ref = refline(1.07,0);
h_ref.Color = [.5 .5 .5];
h_ref.LineWidth = 1.5;

lgd = legend([s1 s2],[turb_method{1}],...
    [turb_method{2}],'Location','best');
set(lgd,'Interpreter','latex','fontsize', 22)

xlabel('friction velocity $$u_*$$', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
ylabel('vertical velocity $$w_{rms}$', 'fontname',...
    'computer modern', 'fontsize', 28,'Interpreter', 'latex')
set(gca,'fontsize',20,'fontname','computer modern','gridlinestyle','--',...
    'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')

export_fig('./figs/w_ustar','-eps','-transparent','-painters')


