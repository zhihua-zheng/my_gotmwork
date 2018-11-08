[spc,sys]=read_cdip_buoy('PAPA_waveriderspectra_09Feb2011-16Jun2011','all');
WRstkxdf=-(2/9.8)*(2*pi*spc.fr).^3.*spc.en.*spc.b1.*spc.bw;
%rfix:
WRstkydf=-(2/9.8)*(2*pi*spc.fr).^3.*spc.en.*spc.a1.*spc.bw;
%fix WRstkydf=(2/9.8)*(2*pi*spc.fr).^3.*spc.en.*spc.a1.*spc.bw;
%rfix:
WRh2df=spc.bw.*spc.en;
WRt2011=sys.tme-datenum(2011,1,0);

% problem with non-distinct WRt2011:
[WRt2011,ilast,jlast] = unique(WRt2011,'last');
WRstkxdf=WRstkxdf(:,ilast);
WRstkydf=WRstkydf(:,ilast);
WRh2df=WRh2df(:,ilast);

t2011=WRt2011(1):median(diff(WRt2011)):WRt2011(end);

load unfixed_Papa_forcing.mat
pft2011=time-datenum(2011,1,0);
HF_LwSL=HF_net-HF_shortwave;
%HF_net(isnan(HF_net))=interp1(pft2011(~isnan(HF_net)),HF_net(~isnan(HF_net)),pft2011(isnan(HF_net)));
HF_LwSL(isnan(HF_LwSL))=interp1(pft2011(~isnan(HF_LwSL)),HF_LwSL(~isnan(HF_LwSL)),pft2011(isnan(HF_LwSL)));
wr_HF_LwSL=interp1(pft2011,HF_LwSL,t2011,'linear',mean(HF_LwSL));

ustar_bulk(isnan(ustar_bulk))=sqrt(interp1(pft2011(~isnan(ustar_bulk)),ustar_bulk(~isnan(ustar_bulk)).^2,pft2011(isnan(ustar_bulk))));
W(isnan(W))=interp1(pft2011(~isnan(W)),W(~isnan(W)),pft2011(isnan(W)));

[z,sorad_pf]=soradna1(pft2011-1,2011,145,50);
[z,sorad]=soradna1(t2011-1,2011,145,50);
small=min(sorad(sorad~=0))*1e-6;
%wr_HF_shortwave=(sorad'+small).*interp1(pft2011,HF_shortwave./(sorad_pf'+small),t2011,'linear',0);
wr_HF_shortwave=(sorad'+small).*interp1(pft2011(~isnan(HF_shortwave)),HF_shortwave(~isnan(HF_shortwave))./(sorad_pf(~isnan(HF_shortwave))'+small),t2011,'linear',mean(HF_shortwave(~isnan(HF_shortwave))./(sorad_pf(~isnan(HF_shortwave))'+small)));

%wr_HF_LwSL=interp1(pft2011,HF_net-HF_shortwave,t2011,'linear',0);

wr_ustar=sqrt(interp1(pft2011,ustar_bulk.^2,t2011,'linear',0));
wr_ustar_bulk=interp1(pft2011,ustar_bulk,t2011,'linear',0);
wr_W=interp1(pft2011,W,t2011,'linear',0.00001+i*0.00001);

%fix: M=[t2011; wr_HF_shortwave; wr_HF_LwSL; 1023*wr_ustar_bulk.^2.*real(wr_W)./abs(wr_W); 1023*wr_ustar_bulk.^2.*imag(wr_W)./abs(wr_W); zeros(size(t2011))]';
M=[t2011; wr_HF_shortwave; wr_HF_LwSL; 1023*wr_ustar_bulk.^2.*real(wr_W)./abs(wr_W); 1023*wr_ustar_bulk.^2.*imag(-wr_W)./abs(wr_W); zeros(size(t2011))]';
save('rfix_papa2011.dat','M','-ascii', '-double')

fctr=spc.fr(:,1);
fctrpad=[-999. fctr'];

stkxdf=(1e-20+ones(size(fctr))*wr_ustar_bulk).*interp1(WRt2011',(WRstkxdf./(ones(size(fctr))*interp1(t2011,1e-20+wr_ustar_bulk,WRt2011,'linear','extrap')))',t2011')';
stkydf=(1e-20+ones(size(fctr))*wr_ustar_bulk).*interp1(WRt2011',(WRstkydf./(ones(size(fctr))*interp1(t2011,1e-20+wr_ustar_bulk,WRt2011,'linear','extrap')))',t2011')';
h2df=(1e-20+ones(size(fctr))*wr_ustar_bulk).^4.*interp1(WRt2011',(WRh2df./(1e-20+ones(size(fctr))*interp1(t2011,1e-20+wr_ustar_bulk,WRt2011,'linear','extrap')).^4)',t2011')';
sdirad=atan2(stkxdf,stkydf);
stksdf=sqrt(stkxdf.^2+stkydf.^2);

t_h2df=[t2011; h2df];
fpad_t_h2df=[fctrpad; t_h2df'];
save('rfix_papa2011_h2df.dat','fpad_t_h2df','-ascii', '-double')

t_stkxdf=[t2011; stkxdf];
fpad_t_stkxdf=[fctrpad; t_stkxdf'];
save('rfix_papa2011_stkxdf.dat','fpad_t_stkxdf','-ascii', '-double')

t_stkydf=[t2011; stkydf];
fpad_t_stkydf=[fctrpad; t_stkydf'];
save('rfix_papa2011_stkydf.dat','fpad_t_stkydf','-ascii', '-double')

t_stksdf=[t2011; stksdf];
fpad_t_stksdf=[fctrpad; t_stksdf'];
save('rfix_papa2011_stksdf.dat','fpad_t_stksdf','-ascii', '-double')

t_sdirad=[t2011; sdirad];
fpad_t_sdirad=[fctrpad; t_sdirad'];
save('rfix_papa2011_sdirad.dat','fpad_t_sdirad','-ascii', '-double')

%>> plot(t2011,atan2(real(wr_W),-imag(wr_W)),'.')
%>> axis auto
%>> plot(t2011, sdirad(end,:),'r')
