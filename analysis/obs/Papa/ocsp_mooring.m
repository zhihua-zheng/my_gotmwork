%% ocsp_mooring

% Script to annalyze mooring data from OCSP

% Zhihua Zheng, UW-APL, July 18 2019

clear

%% Loading

ocsp_root = '~/GDrive/UW/Research/Data/OCSP/';
Mname     = fullfile(ocsp_root,'Mooring/ocsp_prof_hrBox.mat');

load(Mname)

nDEP    = 12;
DEPname = num2str((3:nDEP)','PA-%03d');
DEPdstr = {'13-Jun-2009'
           '14-Jun-2010'
           '11-Jun-2011'
           '02-Jun-2012'
           '17-Jun-2013'
           '18-Jun-2014'
           '15-Jun-2015'
           '15-Jun-2016'
           '14-Jun-2017'
           '22-Jul-2018'
           '13-Jun-2019'};
DEPtime   = datenum(DEPdstr);
DEPdatm   = datetime(DEPtime,'ConvertFrom','datenum');
dateBound = [datetime('20-Aug-2009') datetime('01-Sep-2009');
             datetime('01-Apr-2011') datetime('14-Apr-2011')
             datetime('13-Sep-2011') datetime('29-Sep-2011')
             datetime('28-Feb-2013') datetime('16-Mar-2013')
             datetime('17-Apr-2014') datetime('29-Apr-2014')
             datetime('01-Apr-2015') datetime('26-Apr-2015')
             datetime('04-Apr-2016') datetime('17-Apr-2016');
             datetime('01-Apr-2017') datetime('20-Apr-2017');
             datetime('18-Apr-2018') datetime('01-May-2018');
             datetime('13-Apr-2019') datetime('26-Apr-2019')];

% note there is no rain measurements in PA-0012, air-sea fluxes estimations
% have bigger uncertainties!

%% Choose period

for iD = 1:9
    
    indx   = PROF.datm >= DEPdatm(iD) & PROF.datm < DEPdatm(iD+1);
    iPROF  = PROF(indx,:);

%   iPROF3 = retime(iPROF,'regular','mean','TimeStep',hours(3));
    PTprof  = iPROF.PTprof';
    SAproff = iPROF.SAprof';
    PDprof  = iPROF.PDprof';
    NSQprof = iPROF.NSQprof';
    idatm   = iPROF.datm;
    itime   = datenum(idatm);
    
%   idatmE = (DEPdatm(iD):hours(3):DEPdatm(iD+1))';
%   idatm  = (DEPdatm(iD)+hours(1):hours(3):DEPdatm(iD+1)-hours(2))';
%   itime  = datenum(idatm);

%% Temperature profile

dPTprof = -diff(PTprof)';

mld = get_mld(flipud(PDprof),-flipud(depth_t),1,-1);
sld = mld/5;

Hov = 0;
if Hov
    [TM,Z] = meshgrid(itime,-depth_t);
    figure('position',[0 0 600 300])
    contourf(TM,Z,PTprof,'linestyle','none')
    hold on
    plot(itime,-mld,'--k')
    datetick('x','mm');
end

%% Temperature difference with old prediction

nz    = length(depth_t);
ntm   = length(idatm);
kappa = 0.4;

% note uncertainties in actual measurement depth
dz   = diff(depth_t);
ze   = -(depth_t(2:end) + depth_t(1:end-1))/2; % evaluation depth
tref = linspace(-.1,.2);

k = 1; % PA-009 is good example for k=2
% PA-006, 009, 010 behave a bit differently than others

dT_k       = dPTprof(:,k);
[Lo,Bf,FS] = MOSTpar_from_flux(idatm,ze(k),'Papa');
zeta       = abs(ze(k)) ./ Lo;
[phi_t,~]  = get_theo_phi(zeta); % midpoint approximation

% gradient Richardson number, shear from law of wall
SSQprof = (repmat(FS.Ustar',nz-1,1) ./ repmat(-ze,1,ntm) /kappa).^2;
Rig = NSQprof ./ SSQprof;

% integral approach
zetaU      = depth_t(k)   ./ Lo;
zetaD      = depth_t(k+1) ./ Lo;
zeta_i     = linspaceNDim(zetaD,zetaU);
dzeta      = zeta_i(:,2) - zeta_i(:,1);
phi_ti     = get_theo_phi(zeta_i);
intgrd     = phi_ti./zeta_i;
IntPhiT     = -FS.Tstar .* dzeta .* trapz(intgrd,2);

zeta_cr    =  2;
zeta_cl    = -2;
Ishallow   = find(sld < depth_t(k+1));
Ibig       = find(zeta > zeta_cr | zeta < zeta_cl); % extremely stable/unstable
Ilaminar   = find(Rig(k,:) > .25); % laminar case
Iexc       = union(union(Ishallow,Ibig),Ilaminar);

% dT_MOg = phi_t(Iana) .* FS.Tstar(Iana) ./ (-zt) * dz1;
dT_MOi = IntPhiT;
dT_obs = dT_k;
dT_MOi(Iexc) = NaN;
dT_obs(Iexc) = NaN;

%%

figure('position',[0 0 1200 400])
subplot('position',[0.1 0.15 0.5 0.75])
pos_neg(idatm,dT_MOi,0);plot(idatm,dT_k,'color','#7E2F8E','linewidth',1.8)
legend({'MOST prediction - unstable','MOST prediction - stable',['Mooring obs. ',DEPname(iD,:)]},'fontsize',14,'location','best')
xlim(dateBound(iD,:)); ylim([-.06 .12])
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
plot([-.1 .2],[-.1 .2],'--k','linewidth',2)
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle','--')
xlim([-.02 .04])
ylim([-.02 .04])
xlabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [C]'],'fontsize',14)
ylabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [C]'],'fontsize',14)
title([DEPname(iD,:),' all'],'fontsize',14)

% % bin average
% bin_xi = [(-.017:.002:.017) (.021:.004:.049) (0.057:.01:.087)]';
% % bin_xi = (-.021:.004:.087)';
% Sobs   = pinBin(bin_xi,dT_MOi,dT_obs);
% Sthe   = pinBin(bin_xi,dT_MOi,dT_MOi);
% sl_bin = Sthe.qm\Sobs.qm;
% 
% figure('position',[0 0 600 500])
% errorbar(Sthe.qm,Sobs.qm,Sobs.qerr,'s','MarkerSize',8,...
%     'MarkerFaceColor','r','MarkerEdgeColor','none','linestyle','none');
% hold on
% plot(bin_xi,sl_bin*bin_xi,'Color',rgb('strawberry'),'lineWidth',1.5)
% plot([-.1 .2],[-.1 .2],':k')
% plot([-.1 .2],[ 0   0],':k')
% plot([ 0   0],[-.1 .2],':k')
% box on
% grid on
% set(gca,'gridlinestyle',':')
% xlim([-.02 .09])
% ylim([-.02 .02])
% xlabel('\Delta \Theta from MOST [C]','fontsize',14)
% ylabel('\Delta \Theta measured [C]','fontsize',14)
% title('Temperature difference (1-5 m) Papa Mooring 2011','fontsize',14)

%% Temperature difference with old prediction - focus on night time

dT_MOin = dT_MOi(FS.Inighti);
dT_obsn = dT_obs(FS.Inighti);

goodQn  = and(~isnan(dT_MOin),~isnan(dT_obsn));
sln     = polyfit(dT_MOin(goodQn),dT_obsn(goodQn),1);

figure('position',[0 0 315 880])
subplot(3,1,1)
scatter(dT_MOin,dT_obsn,20,'filled')
hold on
plot(tref,polyval(sln,tref),'Color',rgb('bright sky blue'),'LineWidth',2)
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle',':')
xlim([-.03 .01])
ylim([-.03 .01])
text(-.028,.003,['slope = ',num2str(round(sln(1),2))],'color',rgb('bright sky blue'),'fontsize',16)
xlabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [C]'],'fontsize',14)
ylabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [C]'],'fontsize',14)
title(['Night time ',DEPname(iD,:)],'fontsize',14)

% bin average
switch k
    case 1
        bin_xi = (-.011:.001:.003)';   % for ze(1)
        xbound = [-.02 .007];
        ybound = [-.02 .007];
    case 2
        bin_xi = (-.005:.001:.002)';   % for ze(2)
        xbound = [-.007 .003];
        ybound = [-.007 .003];
    case 3
        bin_xi = (-.0035:.001:.0015)'; % for ze(3)
        xbound = [-.02 .005];
        ybound = [-.02 .005];
end

bin_x   = (bin_xi(1:end-1) + bin_xi(2:end))/2;
Sobs    = pinBin(bin_xi,dT_MOin,dT_obsn,'left');
Sthe    = pinBin(bin_xi,dT_MOin,dT_MOin,'left');

goodQn_bin = and(~isnan(Sthe.qm),~isnan(Sobs.qm));
sln_bin    = polyfit(Sthe.qm(goodQn_bin),Sobs.qm(goodQn_bin),1);

subplot(3,1,2)
errorbar(Sthe.qm,Sobs.qm,Sobs.qerr,'o','MarkerSize',6,...
    'MarkerFaceColor','r','MarkerEdgeColor','r','linestyle','none');
hold on
plot(tref,polyval(sln_bin,tref),'Color',rgb('coral'),'lineWidth',1.5)
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle',':')
xlim(xbound)
ylim(ybound)
text(-.018,.002,['slope = ',num2str(round(sln_bin(1),2))],'color',rgb('coral'),'fontsize',16)
xlabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [C]'],'fontsize',14)
ylabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [C]'],'fontsize',14)
title(['Night time ',DEPname(iD,:)],'fontsize',14)

% nightly average
night_bin = itime(FS.Inighte);
dT_obsDN  = pinBin(night_bin,itime,dT_obs,'none');
dT_obsN   = dT_obsDN.qm(1:2:end);
dT_MOiDN  = pinBin(night_bin,itime,dT_MOi,'none');
dT_MOiN   = dT_MOiDN.qm(1:2:end);

goodQN = and(~isnan(dT_MOiN),~isnan(dT_obsN));
slN    = polyfit(dT_MOiN(goodQN),dT_obsN(goodQN),1);

subplot(3,1,3)
scatter(dT_MOiN,dT_obsN,20,'filled')
hold on
plot(tref,polyval(slN,tref),'Color',rgb('bright sky blue'),'LineWidth',2)
plot([-.1 .2],[-.1 .2],'--k')
plot([-.1 .2],[ 0   0],':k')
plot([ 0   0],[-.1 .2],':k')
box on
grid on
set(gca,'gridlinestyle',':')
xlim([-.02 .005])
ylim([-.02 .005])
text(-.018,.002,['slope = ',num2str(round(slN(1),2))],'color',rgb('bright sky blue'),'fontsize',16)
xlabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) from MOST [C]'],'fontsize',14)
ylabel(['\Delta \Theta (',num2str(depth_t(k)),' - ',num2str(depth_t(k+1)),' m) measured [C]'],'fontsize',14)
title(['Nightly average ',DEPname(iD,:)],'fontsize',14)

end
%%