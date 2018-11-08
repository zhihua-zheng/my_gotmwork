%[spc,sys]=read_cdip_buoy('PAPA_waveriderspectra_Aug-Sep2012','all');
load PAPA_waveriderspectra_Aug-Sep2012.mat
WRstkxdf=-(2/9.8)*(2*pi*spc.fr).^3.*spc.en.*spc.b1.*spc.bw;
WRstkydf=-(2/9.8)*(2*pi*spc.fr).^3.*spc.en.*spc.a1.*spc.bw;
WRh2df=spc.bw.*spc.en;
WRt2012=sys.tme-datenum(2012,1,0);

% problem with non-distinct WRt2012:
[WRt2012,ilast,jlast] = unique(WRt2012,'last');
WRstkxdf=WRstkxdf(:,ilast);
WRstkydf=WRstkydf(:,ilast);
WRh2df=WRh2df(:,ilast);

t2012=WRt2012(1):median(diff(WRt2012)):WRt2012(end);

load Papa_forcing_C20.mat
pft2012=time-datenum(2012,1,0);
HF_LwSL=HF_net-HF_shortwave;
%HF_net(isnan(HF_net))=interp1(pft2012(~isnan(HF_net)),HF_net(~isnan(HF_net)),pft2012(isnan(HF_net)));
HF_LwSL(isnan(HF_LwSL))=interp1(pft2012(~isnan(HF_LwSL)),HF_LwSL(~isnan(HF_LwSL)),pft2012(isnan(HF_LwSL)));
wr_HF_LwSL=interp1(pft2012,HF_LwSL,t2012,'linear',mean(HF_LwSL));

ustar_bulk(isnan(ustar_bulk))=sqrt(interp1(pft2012(~isnan(ustar_bulk)),ustar_bulk(~isnan(ustar_bulk)).^2,pft2012(isnan(ustar_bulk))));
W(isnan(W))=interp1(pft2012(~isnan(W)),W(~isnan(W)),pft2012(isnan(W)));

[z,sorad_pf]=soradna1(pft2012-1,2012,145,50);
[z,sorad]=soradna1(t2012-1,2012,145,50);
small=min(sorad(sorad~=0))*1e-6;
%wr_HF_shortwave=(sorad'+small).*interp1(pft2012,HF_shortwave./(sorad_pf'+small),t2012,'linear',0);
wr_HF_shortwave=(sorad'+small).*interp1(pft2012(~isnan(HF_shortwave)),HF_shortwave(~isnan(HF_shortwave))./(sorad_pf(~isnan(HF_shortwave))'+small),t2012,'linear',mean(HF_shortwave(~isnan(HF_shortwave))./(sorad_pf(~isnan(HF_shortwave))'+small)));

%wr_HF_LwSL=interp1(pft2012,HF_net-HF_shortwave,t2012,'linear',0);

wr_ustar=sqrt(interp1(pft2012,ustar_bulk.^2,t2012,'linear',0));
wr_ustar_bulk=interp1(pft2012,ustar_bulk,t2012,'linear',0);
% direction of W screwed up after June 2011 %  wr_W=interp1(pft2012,W,t2012,'linear',0.00001+i*0.00001);

xcomp=interp1(WRt2012,sin(atan2(WRstkxdf(end,:),WRstkydf(end,:))),t2012,'linear',0);
ycomp=interp1(WRt2012,cos(atan2(WRstkxdf(end,:),WRstkydf(end,:))),t2012,'linear',0);

% direction of W screwed up after June 2011 %  M=[t2012; wr_HF_shortwave; wr_HF_LwSL; 1023*wr_ustar_bulk.^2.*real(wr_W)./abs(wr_W); 1023*wr_ustar_bulk.^2.*imag(wr_W)./abs(wr_W); zeros(size(t2012))]';

M=[t2012; wr_HF_shortwave; wr_HF_LwSL; 1023*wr_ustar_bulk.^2.*xcomp; 1023*wr_ustar_bulk.^2.*ycomp; zeros(size(t2012))]';
save('papa2012.dat','M','-ascii', '-double')

fctr=spc.fr(:,1);
fctrpad=[-999. fctr'];

stkxdf=(1e-20+ones(size(fctr))*wr_ustar_bulk).*interp1(WRt2012',(WRstkxdf./(ones(size(fctr))*interp1(t2012,1e-20+wr_ustar_bulk,WRt2012,'linear','extrap')))',t2012')';
stkydf=(1e-20+ones(size(fctr))*wr_ustar_bulk).*interp1(WRt2012',(WRstkydf./(ones(size(fctr))*interp1(t2012,1e-20+wr_ustar_bulk,WRt2012,'linear','extrap')))',t2012')';
h2df=(1e-20+ones(size(fctr))*wr_ustar_bulk).^4.*interp1(WRt2012',(WRh2df./(1e-20+ones(size(fctr))*interp1(t2012,1e-20+wr_ustar_bulk,WRt2012,'linear','extrap')).^4)',t2012')';
sdirad=atan2(stkxdf,stkydf);
stksdf=sqrt(stkxdf.^2+stkydf.^2);

t_h2df=[t2012; h2df];
fpad_t_h2df=[fctrpad; t_h2df'];
save('papa2012_h2df.dat','fpad_t_h2df','-ascii', '-double')

t_stkxdf=[t2012; stkxdf];
fpad_t_stkxdf=[fctrpad; t_stkxdf'];
save('papa2012_stkxdf.dat','fpad_t_stkxdf','-ascii', '-double')

t_stkydf=[t2012; stkydf];
fpad_t_stkydf=[fctrpad; t_stkydf'];
save('papa2012_stkydf.dat','fpad_t_stkydf','-ascii', '-double')

t_stksdf=[t2012; stksdf];
fpad_t_stksdf=[fctrpad; t_stksdf'];
save('papa2012_stksdf.dat','fpad_t_stksdf','-ascii', '-double')

t_sdirad=[t2012; sdirad];
fpad_t_sdirad=[fctrpad; t_sdirad'];
save('papa2012_sdirad.dat','fpad_t_sdirad','-ascii', '-double')

%>> plot(t2012,-atan2(-real(wr_W),-imag(wr_W)),'.')
%>> axis auto
%>> plot(t2012, sdirad(end,:),'r')
