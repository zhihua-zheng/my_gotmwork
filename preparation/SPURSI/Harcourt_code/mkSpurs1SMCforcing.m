path='~/Documents/GitLab/GOTM_dev/gotmwork/preparation/SPURSI/mat';
files=dir(path);

%assumes directory contains only wave mat files plus '.' and '..' first.
files=files(3:end);

nsp=length(files);

%from the WORD file. last 3 entries of freq. bin centers are dummy 0's

fctr=[.02 
      .0325
      .0375
      .0425
      .0475
      .0525
      .0575
      .0625
      .0675
      .0725
      .0775
      .0825
      .0875
      .0925
      .1
      .11
      .12
      .13
      .14
      .15
      .16
      .17
      .18
      .19
      .2
      .21
      .22
      .23
      .24
      .25
      .26
      .27
      .28
      .29
      .3
      .31
      .32
      .33
      .34
      .35
      .365
      .385
      .405
      .425
      .445
      .465
      .485
      0
      0
      0];

nf=length(fctr);
nfg=nf-3;
fbin = zeros(length(fctr)+1,1);
g = 9.81;

% assume fctr=0.5*(fbin(1:end-1)+fbin(2:end));
% and that the 1st bin starts at
fbin(1)=0.01;
for ifr=1:nf, fbin(ifr+1)=2*fctr(ifr)-fbin(ifr); end
df=diff(fbin);
%1st band is for noise & uses instead 
df(1)=0.010; % I didn't apply this to freq.mat

%last 3 are dummies:
df(end-2:end)=0;
      
for isp=1:nsp
load([path '/' files(isp).name]);
% Directional Wave Spectrum	= C11(f) * D(f,A)
% D(f,A) = (1/PI)*(0.5+R1*COS(A-ALPHA1)+R2*COS(2*(A-ALPHA2)))

tmatspc(isp)=datenum(files(isp).name(7:end-4),'yyyymmddHHMMSS');

SPEC1D(:,isp)=c11;
% Energy from spectrum, assuming frequency is in m*m/Hz
eta2(isp)=sum(c11.*df);
Hsig(isp)=metob.WVHGT;
[maxval, ifmax]=max(c11);
fpk(isp)=fctr(ifmax);
STKDIR(:,isp)=alpha1;
SPECR1(:,isp)=R1;

% minus sign is because wave spectra is for direction from which waves are propagating.
STKX(:,isp) = -(16*pi^3/g)*c11.*R1.*sin(alpha1*pi/180) .* fctr.^3;
STKY(:,isp) = -(16*pi^3/g)*c11.*R1.*cos(alpha1*pi/180).*fctr.^3;

% mult. Stokes spectra by df for model forcing purposes.
DSTKX(:,isp) = -(16*pi^3/g)*df.*c11.*R1.*sin(alpha1*pi/180).*fctr.^3;
DSTKY(:,isp) = -(16*pi^3/g)*df.*c11.*R1.*cos(alpha1*pi/180).*fctr.^3;

end

load ~/Documents/GitLab/GOTM_dev/gotmwork/preparation/SPURSI/spurs1_met_1hr.mat
load ~/Documents/GitLab/GOTM_dev/gotmwork/preparation/SPURSI/spurs1_flux_1hr.mat

%Energy & Frequency scales:
%PM limit is eta2/Escale=3.64e-3 ; fpeak/Fscale=0.13

Escale1h=interp1(mday,ws_h,tmatspc).^4/g^2;
Fscale1h=interp1(mday,ws_h,tmatspc)/g;

%%%%%%%%%%%%%%%
%patch into older PAPA routine to force SMC
%%%%%%%%%%%%%%

t2012=tmatspc-datenum(2012,1,0);
stkxdf=(df*ones(size(tmatspc))).*STKX;
stkydf=(df*ones(size(tmatspc))).*STKY;
h2df=(df*ones(size(tmatspc))).*SPEC1D;
sdirad=atan2(stkxdf,stkydf); % swap x and y ?
stksdf=sqrt(stkxdf.^2+stkydf.^2);


mft2012=mday-datenum(2012,1,0);
HF_LwSL=QN-Qs;
% do fancy interp to spectral time:

ndb_HF_LwSL=interp1(mft2012,HF_LwSL,t2012,'pchip',mean(HF_LwSL));
ndb_HF_sw=interp1(mft2012,Qs,t2012,'pchip',mean(Qs));
ndb_tau=interp1(mft2012,taumag,t2012,'pchip',0);

xcomp=interp1(mft2012,sin(taudir*pi/180),t2012,'linear',0);
ycomp=interp1(mft2012,cos(taudir*pi/180),t2012,'linear',0);

[cumEvap,cumPrecip]=ep(prate/60,QH);
ndb_emp=interp1(mft2012,-cumEvap-cumPrecip,t2012,'pchip',-cumEvap(end)-cumPrecip(end)); % -E-P ?

M=[t2012; ndb_HF_sw; ndb_HF_LwSL; ndb_tau.*sin(atan2(xcomp,ycomp)); ndb_tau.*cos(atan2(xcomp,ycomp)); ndb_emp]';
save('spurs1.dat','M','-ascii', '-double')

save spurs1.mat t2012 ndb_HF_sw ndb_HF_LwSL ndb_tau xcomp ycomp

fctrpad=[-999. fctr(1:nfg)'];

t_h2df=[t2012; h2df(1:nfg,:)];
fpad_t_h2df=[fctrpad; t_h2df'];
save('spurs1_h2df.dat','fpad_t_h2df','-ascii', '-double')

t_stkxdf=[t2012; stkxdf(1:nfg,:)];
fpad_t_stkxdf=[fctrpad; t_stkxdf'];
save('spurs1_stkxdf.dat','fpad_t_stkxdf','-ascii', '-double')

t_stkydf=[t2012; stkydf(1:nfg,:)];
fpad_t_stkydf=[fctrpad; t_stkydf'];
save('spurs1_stkydf.dat','fpad_t_stkydf','-ascii', '-double')

t_stksdf=[t2012; stksdf(1:nfg,:)];
fpad_t_stksdf=[fctrpad; t_stksdf'];
save('spurs1_stksdf.dat','fpad_t_stksdf','-ascii', '-double')

t_sdirad=[t2012; sdirad(1:nfg,:)];
fpad_t_sdirad=[fctrpad; t_sdirad'];
save('spurs1_sdirad.dat','fpad_t_sdirad','-ascii', '-double')

%>> plot(t2012,atan2(real(ndb_W),-imag(ndb_W)),'.')
%>> axis auto
%>> plot(t2012, sdirad(end,:),'r')
