%% cdip_spec_check

% Zhihua Zheng, June 21 2019

%% Loading

ana_root = '~/Documents/GitHub/GOTM/gotmwork/analysis/';
moor_dir = '~/GDrive/UW/Research/Data/OCSP/Mooring/';
cdip_dir = '~/GDrive/UW/Research/Data/OCSP/CDIP/';

METname  = fullfile(moor_dir,'Met_2007_2019.mat');
CDIPname = fullfile(cdip_dir,'Wave_2010.mat');

load(METname,'w_spd','w_dir','time','z_wind','Ta');
load(CDIPname)

w_spd(w_spd>100)  = NaN;
Ta(Ta>100)        = NaN;
w_dir(w_dir>1000) = NaN;

w_dir = deg2rad(90 - w_dir); % counter-clockwise from East (radian)

%%

g  = 9.81;
df = repmat(wave_bw,  1,length(wave_time));
f  = repmat(wave_freq,1,length(wave_time));

[u10,ustar] = spshfttc(w_spd,z_wind,10,Ta); % following TOGA/COARE, Smith (1988)

% project to wave time
u10_wt   = interp1(time,u_10,wave_time);
ustar_wt = interp1(time,ustar,wave_time);
w_dir_wt = interp1(time,w_dir,wave_time);

% nondimensional wave energy: epsilon = eta2 * g^2 / (u10^4)
eta2    = sum(wave_spec .* df); % [m^2]
epsilon = eta2*g^2 ./ (u10_wt.^4)';

eps_s  = epsilon(epsilon<1);
s_frac = round(length(eps_s)/length(epsilon),2)*100;
eps_b  = epsilon(epsilon>=1);
b_frac = round(length(eps_b)/length(epsilon),2)*100;

edges = [linspace(0,1e0,1000) 8e4];
figure
histogram(epsilon,edges,'Normalization','countdensity')
set(gca,'Xscale','log')
xlim([1e-4 1e4])
xlabel('\eta^2 g^2/u_{10}^4','fontsize',14)
ylabel('Count Density','fontsize',14)
title([num2str(b_frac),'% of total is larger than 1'],'fontsize',14)
saveas(gcf,[ana_root,'obs/figs/nondimensional_waveE'], 'png');

% subplot(2,1,1)
% histogram(eps_s,'Normalization','probability')
% title(['\eta^2 g^2/u_{10}^4 < 1,  ',num2str(s_frac),'% of total'],'fontsize',14)
% % 62% ~ 1e-3
% subplot(2,1,2)
% histogram(eps_b,'Normalization','probability')
% title(['\eta^2 g^2/u_{10}^4 >= 1,  ',num2str(b_frac),'% of total'],'fontsize',14)

% surface Stokes drift
const    = 16*pi^3/g;
ux_St0   = -const*sum(f.^3.*wave_spec.*wave_b1.*df);
uy_St0   = -const*sum(f.^3.*wave_spec.*wave_a1.*df);
u_St0    = sqrt(ux_St0.^2 + uy_St0.^2);
uSt0_dir = atan2(uy_St0,ux_St0); % counter-clockwise from East

% turbulent Langmuir number
La_t = sqrt(u_St0 ./ ustar_wt');
figure
histogram(La_t,'Normalization','probability')
xlabel('La_t','fontsize',14)
ylabel('Probability','fontsize',14)
saveas(gcf,[ana_root,'obs/figs/La_t'], 'png');

% direction check
wd_we_diff = w_dir_wt - uSt0_dir'; % direction difference
figure
polarhistogram(wd_we_diff,'Normalization','probability','FaceColor','red')
title('probability of direction difference: wind - surface Stokes drift','fontsize',14)
saveas(gcf,[ana_root,'obs/figs/wd_we_dir'], 'png');

