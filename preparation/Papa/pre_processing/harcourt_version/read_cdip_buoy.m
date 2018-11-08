function [tspc,tsys]=read_cdip_buoy(fle,tmestr)
%
%   [spc,sys]=read_cdip_buoy(fle,tmestr)
%
% Function to read spectral wavebuoy data from a file downloaded from the CDIP 
% (Coastal Data Information Program) web page.
%
% The function returns two data structures (spc and sys) when supplied with the data
% file name (fle) and the date/time string (tmestr) of the desired spectrum. The
% date/time string MUST be of the form 'YYYYMMDDHHmm' or as noted below.
%
% The structure 'spc' contains the frequency (fr), bandwidth (bw),
% energy density (en) (m^2/Hz), mean direction (Dm), a1, b1, a2, 
% b2, and the check factor (ck) at each of the 64 frequency bands. 
%
% The structure variables are M x N arrays where the rows (M) represent 
% values for different frequencies within a spectra and the columns (N)
% represent different spectra.
% 
% The 'sys' structure contains the following information:
%
%   Hs - significant wave height (m)
%   lat - latitude
%   lon - longitude
%   Tp  - peak period (s)
%   Dp  - mean direction (deg)
%   Ta  - ? (s)
%   tme - spectrum time (in matlab datenum format)
%
% NOTES:
% The user may elect to read all of the spectra contained in the file by specifying 
% the date/time string as 'all'.
%
%  The directional moments a1, a2, b1, b2 as supplied by CDIP are in a left-hand
%  coordinate system where +X is south and +Y is west. These are converted to
%  our system (right handed) where +X is north and +Y is west.
%
%  Version: 1.5 - 4/2004,              Paul Jessen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Written by: Paul F. Jessen
%              Department of Oceanography
%              Naval Postgraduate School
%              Monterey, CA
%
%  First Version: 1.0 - 11/2000
%
%  Latest Version: 1.5 - 4/2004
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Modifications:
%
%   Version: 1.1 - Mix up of a2 and b1 colums fixed.
%                - bandwidth variable added (12/4/00 - pfj)
%
%   Version: 1.2 - Changed function so that file name is passed rather
%                  than file handle. File is now opened from within
%                  the function.
%                - now using matlab function "strread" rather than in
%                  house function "ext_field" to extract fields from
%                  mixed string/numeric data records. (2/14/2003 - pfj)
%
%   Version: 1.3 - Function modified to check version of matlab being used.
%                  If version is less than 6, then the external function
%                  "ext_field.m" is used because "strread.m" doesn't exist
%                  until version 6.0 (2/19/2003 - pfj)
%
%   Version: 1.4 - Change "spc" and "sys" structure definition to correspond
%                  to the one created by the "load_wv_mat" function. This 
%                  entailed changing function so structure elements are M X N
%                  arrays where M is different frequencies and N is different
%                  spectra and the structure itself is 1 x 1 rather than 
%                  structure elements being 1 X M where M is different freqs.
%                  and the structure itself being a 1 X N array where N is the
%                  different spectra (4/1/2004 - pfj)
%
%   Version: 1.5 - Change how directional moments are oriented. CDIP (and NDBC)
%                  report the moments in a left handed coordinate system where
%                  +X is south and +Y is west. To convert these to our system
%                  (right handed) where +X is north and +Y is west we have to
%                  change the signs of a1 and b2 (4/20/2004 - pfj)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% check matlab version
%
mvc=version;
mv=str2num(mvc(1:1));
vflg=0;
if mv < 6, vflg=1;, end
%
% create empty matices for the structure variables
%
  Hs_a=[];
  lat_a=[];
  lon_a=[];
  Tp_a=[];
  Dp_a=[];
  Ta_a=[];
  tme_a=[];
  fr_a=[];
  en_a=[];
  a1_a=[];
  a2_a=[];
  b1_a=[];
  b2_a=[];
  bw_a=[];
  Dm_a=[];
  ck_a=[];
%
% check to see if all spectra are to be extracted or only one
%
ext_opt=0;
if strcmp(tmestr,'all') == 1, 
  ext_opt=1;
  dtme=0;
else
  % have to decode the time string to get a number
  yr=str2num(tmestr(1:4));
  mon=str2num(tmestr(5:6));
  dy=str2num(tmestr(7:8));
  hr=str2num(tmestr(9:10));
  mn=str2num(tmestr(11:12));
  sec=0;
  dtme=datenum(yr,mon,dy,hr,mn,sec);
end
%
% open file
%
fid=fopen(fle);
if fid <= 0,
  errordlg('Error Opening Data File!');
  return
end
%
% read first line of file
%
lne=fgetl(fid);
fprintf(1,'\n%s',lne);
%
% main loop
%
while feof(fid) == 0,  
%
% decode time and buoy id from part of string, use my function 
% "ext_field.m" to extract desired fields for matlab versions less than 6
%
if vflg == 1;
  tstr=ext_field(lne,3);
else
  a=strread(lne,'%s','delimiter',' ');
  tstr=char(a{3}); 
end
  ib=find(isletter(tstr) == 0);
  buoy_time=tstr(ib);
%
% separate buoy id and time, time will be used later to see if this
% spectrum is to be used
%
  buoy_id=str2num(buoy_time(1:5));
  tme=buoy_time(6:end);
  yr=str2num(buoy_time(6:9));
  mon=str2num(buoy_time(10:11));
  dy=str2num(buoy_time(12:13));
  hr=str2num(buoy_time(14:15));
  mn=str2num(buoy_time(16:17));
  sec=0;
  % read rest of spectrum
  lne=fgetl(fid);
  lne=fgetl(fid);
  % decode lat, lon
  ia=findstr(lne,' ');
  latd=str2num(lne(ia(1):ia(2)));
  latm=str2num(lne(ia(2):ia(3)));
  lond=str2num(lne(ia(4):ia(5)));
  lonm=str2num(lne(ia(5):ia(6)));
  % skip a few lines
  lne=fgetl(fid);  
  lne=fgetl(fid);  
  lne=fgetl(fid);  
  % decode Hs, Tp, Dp, and Ta
  lne=fgetl(fid);
  %
  % use my function "ext_field" to extract the correct field from the
  % line (for matlab versions less than 6)
  %
  if vflg == 1,
    swh=str2num(ext_field(lne,2));
    peak_per=str2num(ext_field(lne,4));
    peak_dir=str2num(ext_field(lne,6));
    whatever=str2num(ext_field(lne,8));
  else
    a=strread(lne,'%s','delimiter',' ');
    swh=str2num(char(a{2}));
    peak_per=str2num(char(a{4}));
    peak_dir=str2num(char(a{6}));
    whatever=str2num(char(a{8}));
  end
  % skip empty line
  lne=fgetl(fid);
  lne=fgetl(fid);
  lne=fgetl(fid);
  %
  % now get the spectral data
  %
  f=fscanf(fid,'%f %f %f %d %f %f %f %f %f',[9,64]);  
  spectime=datenum(yr,mon,dy,hr,mn,sec);
%
% compare system time with desired time to see if this spectrum is to be 
% used.
% It will also be used if the 'ext_opt' variable was set to 1 at beginning 
% of function. Difference between times of less than 1 minute are ignored
% to avoid possible roundoff errors.
%
%keyboard
  if (abs(spectime - dtme) < 1/1440) | ext_opt == 1,
    tme_a=[tme_a,spectime];
    lat_a=[lat_a,latd+latm/60];
    lon_a=[lon_a,lond+lonm/60];
    Hs_a=[Hs_a,swh];
    Tp_a=[Tp_a,peak_per];
    Dp_a=[Dp_a,peak_dir];
    Ta_a=[Ta_a,whatever];
    %
    % now spectral information
    %
    fr=f(1,:);
    bw=f(2,:);
    en=f(3,:);
    Dm=f(4,:);
    a1=f(5,:);
    % convert to rh system
    a1=-1*a1;
    b1=f(6,:);
    a2=f(7,:);
    b2=f(8,:);
    % convert to rh system
    b2=-1*b2;
    ck=f(9,:);  
    fr_a=[fr_a,fr'];
    en_a=[en_a,en'];
    a1_a=[a1_a,a1'];  
    a2_a=[a2_a,a2'];
    b1_a=[b1_a,b1'];
    b2_a=[b2_a,b2'];
    ck_a=[ck_a,ck'];
    bw_a=[bw_a,bw'];
    Dm_a=[Dm_a,Dm'];
    clear fr en a1 a2 b1 b2 ck bw Dm
    %
    % if this is the only spectrum to extract go ahead and exit,
    % otherwise keep going.
    %
    if ext_opt == 0,
      %
      % create spectral and system structures
      %
      tspc=struct('fr',fr_a,'bw',bw_a,'en',en_a,'Dm',Dm_a,'a1',a1_a,'a2',a2_a,'b1',b1_a,...
       'b2',b2_a,'ck',ck_a);
      tsys=struct('Hs',Hs_a,'lat',lat_a,'lon',lon_a,'Tp',Tp_a,'Dp',Dp_a,'Ta',Ta_a,'tme',tme_a);
      return;
      %    else
%      tspc=[tspc spc];
%      tsys=[tsys sys];
%      clear spc sys
%      spc=struct('fr',zeros(64,1),'bw',zeros(64,1),'en',...
%        zeros(64,1),'Dm',zeros(64,1),'a1',zeros(64,1),'a2',zeros(64,1),...
%      'b1',zeros(64,1),'b2',zeros(64,1),'chk',zeros(64,1));
%      sys=struct('Hs',0,'lat',0,'lon',0,'Tp',0,'Dp',0,'Ta',0,'tme',0);
    end          
  %
  % now there is one empty line before the next spectrum
  %
  end  
  lne=fgetl(fid);
  lne=fgetl(fid);
  fprintf(1,'\n%s',lne);
end
fclose(fid);
%
% have read all the spectra, create structures and exit
%      
tspc=struct('fr',fr_a,'bw',bw_a,'en',en_a,'Dm',Dm_a,'a1',a1_a,'a2',a2_a,'b1',b1_a,...
 'b2',b2_a,'ck',ck_a);
tsys=struct('Hs',Hs_a,'lat',lat_a,'lon',lon_a,'Tp',Tp_a,'Dp',Dp_a,'Ta',Ta_a,'tme',tme_a);
return
