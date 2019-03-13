load ~/Documents/GitLab/GOTM_dev/gotmwork/preparation/SPURSI/SPURS_TS_grid

tinit=264.6;
itidaybin=find( abs((yday-datenum(2012,1,0))-tinit) < 1 );

Z1=Z';
T1=mean(T(itidaybin,:))';
S1=mean(S(itidaybin,:))';

M=[Z1 T1 S1 zeros(size(Z1)) zeros(size(Z1)) 1e-4*ones(size(Z1)) ones(size(Z1))];

save('SpursMoorIC264_6.dat','M','-ascii', '-double')

