%% spursi_mooring

% Script to annalyze mooring data from SPURS-I site

% Zhihua Zheng, UW-APL, July 18 2019

clear

%% Loading

spursi_root = '~/GDrive/UW/Research/Data/SPURSI/';
Mname       = fullfile(spursi_root,'Mooring/spursi_prof_hrBox.mat');

load(Mname)

iPROF   = PROF(1:end-2,:);
iProf   = Prof(1:end-2,:);
PTprof  = iPROF.PTprof';
SAprof  = iPROF.SAprof';
PDprof  = iPROF.PDprof';
Bprof   = iPROF.Bprof';
NSQprof = iPROF.NSQprof';
idatm   = iPROF.datm;
itime   = datenum(idatm);

DEPname   = 'SPURS-I';
dateBound = [datetime('27-May-2013') datetime('05-Jun-2013')];
iD = 1;

% tshift  = timezone(lon);

%% Density profile

mld = get_mld(flipud(PDprof),-flipud(depth_t),1,-1);

sld = mld/5;

Hov = 0;
if Hov
    [TM,Z] = meshgrid(itime,-depth_t);
    figure('position',[0 0 600 300])
    contourf(TM,Z,PDprof,'linestyle','none')
    hold on
    plot(itime,-mld,'--k')
    ylim([-150 0])
    datetick('x','mmm');
end

%% Temperature difference

nz  = length(depth_t);
ntm = length(idatm);

% SPURSI has very strong negative gradients not predicted by MOST
kappa = 0.4;
dz    = diff(depth_t);
ze    = -(depth_t(2:end) + depth_t(1:end-1))/2; % evaluation depth

dPTprof = -diff(PTprof)';

%%

[zl,zu] = meshgrid(-depth_t,-depth_t);
zMat    = cat(3,zu,zl);
bn      = nan(size(zMat)); % coefficients for night time fitting line 
bicn    = nan([size(zMat) 2]); % confidence interval for night time fitting line 
ntfb    = nan(size(zMat)); % info for night time figure bounds
zeta_c  = 2;
Rig_c   = 0.25;

night_title = ['Night time ',DEPname(iD,:)];

ntfb(1,2,:) = [-.008 .002];
ntfb(1,3,:) = [-.013 .001];
ntfb(1,4,:) = [-.016 .001];
ntfb(1,5,:) = [-.020 .001];
ntfb(1,6,:) = [-.021 .001];

ntfb(2,3,:) = [-.008 .002];
ntfb(2,4,:) = [-.011 .002];
ntfb(2,5,:) = [-.015 .002];
ntfb(2,6,:) = [-.016 .002];

ntfb(3,4,:) = [-.007 .004];
ntfb(3,5,:) = [-.009 .004];
ntfb(3,6,:) = [-.009 .004];

ntfb(4,5,:) = [-.005 .002];
ntfb(4,6,:) = [-.006 .002];

ntfb(5,6,:) = [-.004 .002];

%% MOST parameters and obs-theory comparison

for i = 1:2         % index for upper level
    for j = (i+1):3 % index for lower level
        
zz  = squeeze(zMat(i,j,:));
Bzz = Bprof([i,j],:);
        
% -------------------------------------------------------------------------
[Lo,zeta,FS] = MOSTpar_from_flux(idatm,zz,'SPURSI');

% integral approach to estimate difference at two levels
IntPhiT = get_MOST_Delta(FS.Tstar,zeta,Lo,zz);

% gradient Richardson number, shear as law of wall
Rig = get_Rig(Bzz,FS.Ustar',zz);

% subset
zetaM    = max(abs(zeta)); % maximum magnitude of zeta for each time point
Ishallow = find(sld < zz(2));
Ibig     = find(zetaM > zeta_c); % strongly stable/unstable conidtion
Ilaminar = find(Rig > Rig_c); % laminar case
Iexc     = union(union(Ishallow,Ibig),Ilaminar);
% -------------------------------------------------------------------------

dT_MOi = IntPhiT';
dT_obs = (PTprof(i,:) - PTprof(j,:))';

dT_MOi(Iexc) = NaN;
dT_obs(Iexc) = NaN;

% day time

% night time
dT_MOin = dT_MOi(FS.Inighti);
dT_obsn = dT_obs(FS.Inighti);
tboundn = squeeze(ntfb(i,j,:));
xref    = linspace(tboundn(1),tboundn(2));
texti   = tboundn(1) + 8e-4;
textj   = tboundn(1) + 6e-4;
xltext   = ['\Delta \Theta (',num2str(-zz(1)),' - ', ...
            num2str(-zz(2)),' m) from MOST [C]'];
yltext   = ['\Delta \Theta (',num2str(-zz(1)),' - ', ...
            num2str(-zz(2)),' m) measured [C]'];
     
% bin average
[~,bin_xi] = histcounts(dT_MOin,31,'BinLimits',tboundn);
Sy         = pinBin(bin_xi',dT_MOin,dT_obsn,'left');
Sx         = pinBin(bin_xi',dT_MOin,dT_MOin,'left');

Sy.qm(Sy.n < 60) = NaN;

% [bn(i,j,:),bicn(i,j,:,:)] = regress(Sy.qm,[ones(size(Sx.qm)) Sx.qm]);
tmp          = fitlm(Sx.qm,Sy.qm,'linear','RobustOpts','on');
bn(i,j,:)    = tmp.Coefficients.Estimate;
[yref,yCI]   = predict(tmp,xref');

fit_str   = ['slope = ',num2str(round(bn(i,j,2),2)),...
             ', offset = ',num2str(round(bn(i,j,1),4))];

figure('position',[0 0 455 400])
scatter(dT_MOin,dT_obsn,20,'filled','MarkerFaceAlpha',.2)
hold on
fill([xref flip(xref)]',[yCI(:,1);flip(yCI(:,2))],rgb('coral'),'LineStyle','none','FaceAlpha',.2)
bfit = plot(xref,yref,'Color',rgb('coral'),'lineWidth',1.8);
errorbar(Sx.qm,Sy.qm,Sy.qmL,Sy.qmU,'o','MarkerSize',10,'CapSize',10,...
         'Color',rgb('ultramarine'),'linewidth',1.2,'linestyle','none');
scatter(Sx.qm,Sy.qm,85,'filled','MarkerFaceColor',rgb('yellow'))
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[0   0],':k'); plot([0   0],[-.1 .2],':k')
box on; grid on; set(gca,'gridlinestyle',':'); xlim(tboundn); ylim(tboundn)
xlabel(xltext,'fontsize',14); ylabel(yltext,'fontsize',14)
legend(bfit,fit_str,'fontsize',16,'location','best','color','none')
title(night_title,'fontsize',14)
        
    end
end



%%

show_otc = 1;

if show_otc

figure('position',[0 0 1200 400])
subplot('position',[0.1 0.15 0.5 0.75])
pos_neg(idatm,dT_MOi,0);plot(idatm,dT_k,'color','#7E2F8E','linewidth',1.8)
legend({'MOST prediction - unstable','MOST prediction - stable',['Mooring obs. ',DEPname(iD,:)]},'fontsize',14,'location','best')
xlim(dateBound(iD,:)); ylim([-.01 .05])
ylabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m)'],'fontsize',14)

% TEST code
% figure;plot(zeta(Iana))
% figure;plot(itime,dT_obs,itime,dT_MOi)
% figure;plot(itime(1:end-2),dT_obs(3:end),itime,IntPhi)
% figure;scatter(Int(1:end-2),dT_obs(3:end),20,'filled');box on

goodQ = and(~isnan(dT_MOi),~isnan(dT_obs));
sl    = dT_MOi(goodQ) \ dT_obs(goodQ);

subplot('position',[0.7 0.15 0.25 0.75])
scatter(dT_MOi,dT_obs,20,'filled')
hold on
% plot(xref,sl*xref,'Color',rgb('bright sky blue'),'LineWidth',2)
plot([-.1 .2],[-.1 .2],'--k','linewidth',2)
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle','--')
xlim([-.03 .08])
ylim([-.03 .08])
xlabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [C]'],'fontsize',14)
ylabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [C]'],'fontsize',14)
title([DEPname(iD,:),' all'],'fontsize',14)

end

%% Temperature difference vs MOST prediction - focus on night time

% nightly average
night_bin = itime(FS.Inighte);
dT_obsDN  = pinBin(night_bin,itime,dT_obs,'none');
dT_obsN   = dT_obsDN.qm(1:2:end);
dT_MOiDN  = pinBin(night_bin,itime,dT_MOi,'none');
dT_MOiN   = dT_MOiDN.qm(1:2:end);

%% Buoyancy difference vs MOST prediction

show_bs = 0;

if show_bs
    
dBprof  = -diff(Bprof)';
bref    = linspace(-1e-3,.02);

IntPhiB = -FS.Bstar .* dzeta .* trapz(intgrd,2);
dB_MOi  = IntPhiB;
dB_k    = dBprof(:,k);
dB_obs  = dB_k;

dB_MOi(Iexc) = NaN;
dB_obs(Iexc) = NaN;

dB_MOin = dB_MOi(FS.Inighti);
dB_obsn = dB_obs(FS.Inighti);

goodQn  = and(~isnan(dB_MOin),~isnan(dB_obsn));
Bsln    = polyfit(dB_MOin(goodQn),dB_obsn(goodQn),1);

switch k
    case 1
        bin_xi = (-.0053:.0003:-.0014)';   % for ze(1)
        tboundn = [-.007 .002];
        tboundn = [-.007 .002];
end

scatter(dB_MOi,dB_obs,20,'filled')
hold on
% plot(xref,sl*xref,'Color',rgb('bright sky blue'),'LineWidth',2)
plot([-.1 .2],[-.1 .2],'--k','linewidth',2)
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle','--')
xlim([-.0005 .002])
ylim([-.0005 .002])
xlabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [C]'],'fontsize',14)
ylabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [C]'],'fontsize',14)
title([DEPname(iD,:),' all'],'fontsize',14)


figure('position',[0 0 415 880])
subplot(2,1,1)
scatter(dB_MOin,dB_obsn,20,'filled')
hold on
plot(bref,polyval(Bsln,bref),'Color',rgb('bright sky blue'),'LineWidth',2)
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle',':')
xlim(tboundn)
ylim(tboundn)
% text(-.0085,.001,['slope = ',num2str(round(sln(1),2))],'color',rgb('bright sky blue'),'fontsize',16)
xlabel(['\Delta B (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [m/s^2]'],'fontsize',14)
ylabel(['\Delta B (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [m/s^2]'],'fontsize',14)
title(['Night time ',DEPname(iD,:)],'fontsize',14)

% bin average
bin_x   = (bin_xi(1:end-1) + bin_xi(2:end))/2;
Sy    = pinBin(bin_xi,dB_MOin,dB_obsn,'left');
Sx    = pinBin(bin_xi,dB_MOin,dB_MOin,'left');

goodQnBin = and(~isnan(Sx.qm),~isnan(Sy.qm));
Bsln_bin   = polyfit(Sx.qm(goodQnBin),Sy.qm(goodQnBin),1);

subplot(2,1,2)
errorbar(Sx.qm,Sy.qm,Sy.qerr,'o','MarkerSize',6,...
    'MarkerFaceColor','r','MarkerEdgeColor','r','linestyle','none');
hold on
plot(bref,polyval(Bsln_bin,bref),'Color',rgb('coral'),'lineWidth',1.5)
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle',':')
xlim(tboundn)
ylim(tboundn)
fit_str = ['slope = ',num2str(round(Bsln_bin(1),2)),', offset = ',num2str(round(Bsln_bin(2),4))];
text(tboundn(1)+3e-4,5e-4,fit_str,'color',rgb('coral'),'fontsize',16)
xlabel(['\Delta B (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [m/s^2]'],'fontsize',14)
ylabel(['\Delta B (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [m/s^2]'],'fontsize',14)
title(['Night time ',DEPname(iD,:)],'fontsize',14)


%% Salinity difference vs MOST prediction

sref = linspace(-.1,.2);
zeS  = -(depth_s(1) + depth_s(3))/2;

% -------------------------------------------------------------------------
[Lo,zeta,FS] = MOSTpar_from_flux(idatm,zeS,'SPURSI');

% gradient Richardson number, shear from law of wall
SSQprof_pin = (FS.Ustar' ./ repmat(-zeS,1,ntm) / kappa).^2;
NSQprof_pin = center_diff(Bprof,-depth_t',1,'pin');
Rig_pin     = NSQprof_pin(2,:) ./ SSQprof_pin;

% integral approach
zetaU      = depth_t(1) ./ Lo;
zetaD      = depth_t(3) ./ Lo;
zeta_i     = linspaceNDim(zetaD,zetaU);
dzeta      = zeta_i(:,2) - zeta_i(:,1);
phi_ti     = get_theo_phi(zeta_i);
intgrd     = phi_ti ./ zeta_i;

% subset
zeta_cr    =  2;
zeta_cl    = -2;
Ishallow   = find(sld < depth_t(k+1));
Ibig       = find(zeta > zeta_cr | zeta < zeta_cl); % extremely stable/unstable
Ilaminar   = find(Rig_pin > .3); % laminar case
Iexc       = union(union(Ishallow,Ibig),Ilaminar);
% -------------------------------------------------------------------------

IntPhiS = -FS.Sstar .* dzeta .* trapz(intgrd,2);
dS_MOi  = IntPhiS;
dS_k    = SAprof(1,:) - SAprof(3,:);
dS_obs  = dS_k';

dS_MOi(Iexc) = NaN;
dS_obs(Iexc) = NaN;

dS_MOin = dS_MOi(FS.Inighti);
dS_obsn = dS_obs(FS.Inighti);

goodQn  = and(~isnan(dS_MOin),~isnan(dS_obsn));
Ssln    = polyfit(dS_MOin(goodQn),dS_obsn(goodQn),1);

% TEST
figure;plot(idatm,dS_k);hold on;plot(idatm,dS_MOi)


switch k
    case 1
        bin_xi = (-.0053:.0003:-.0014)';   % for ze(1)
        tboundn = [-.007 .002];
        tboundn = [-.007 .002];
end

scatter(dS_MOi,dS_obs,20,'filled')
hold on
% plot(xref,sl*xref,'Color',rgb('bright sky blue'),'LineWidth',2)
plot([-.1 .2],[-.1 .2],'--k','linewidth',2)
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle','--')
xlim([-.0005 .002])
ylim([-.0005 .002])
xlabel(['\Delta S (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [g/kg]'],'fontsize',14)
ylabel(['\Delta S (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [g/kg]'],'fontsize',14)
title([DEPname(iD,:),' all'],'fontsize',14)


figure('position',[0 0 415 880])
subplot(2,1,1)
scatter(dS_MOin,dS_obsn,20,'filled')
hold on
plot(sref,polyval(bn,sref),'Color',rgb('bright sky blue'),'LineWidth',2)
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle',':')
xlim(tboundn)
ylim(tboundn)
% text(-.0085,.001,['slope = ',num2str(round(sln(1),2))],'color',rgb('bright sky blue'),'fontsize',16)
xlabel(['\Delta S (',num2str(depth_t(1)),' - ',num2str(depth_t(3)),' m) from MOST [g/kg]'],'fontsize',14)
ylabel(['\Delta S (',num2str(depth_t(1)),' - ',num2str(depth_t(3)),' m) measured [g/kg]'],'fontsize',14)
title(['Night time ',DEPname(iD,:)],'fontsize',14)

% bin average

bin_x   = (bin_xi(1:end-1) + bin_xi(2:end))/2;
Sy    = pinBin(bin_xi,dS_MOin,dS_obsn,'left');
Sx    = pinBin(bin_xi,dS_MOin,dS_MOin,'left');

goodQnBin = and(~isnan(Sx.qm),~isnan(Sy.qm));
Ssln_bin    = polyfit(Sx.qm(goodQnBin),Sy.qm(goodQnBin),1);

subplot(2,1,2)
errorbar(Sx.qm,Sy.qm,Sy.qerr,'o','MarkerSize',6,...
    'MarkerFaceColor','r','MarkerEdgeColor','r','linestyle','none');
hold on
plot(tref,polyval(Ssln_bin,tref),'Color',rgb('coral'),'lineWidth',1.5)
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle',':')
xlim(tboundn)
ylim(tboundn)
fit_str = ['slope = ',num2str(round(Ssln_bin(1),2)),', offset = ',num2str(round(Ssln_bin(2),3))];
text(tboundn(1)+3e-4,5e-4,fit_str,'color',rgb('coral'),'fontsize',16)
xlabel(['\Delta S (',num2str(depth_t(1)),' - ',num2str(depth_t(3)),' m) from MOST [g/kg]'],'fontsize',14)
ylabel(['\Delta S (',num2str(depth_t(1)),' - ',num2str(depth_t(3)),' m) measured [g/kg]'],'fontsize',14)
title(['Night time ',DEPname(iD,:)],'fontsize',14)

end