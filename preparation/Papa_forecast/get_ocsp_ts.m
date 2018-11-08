function ts = get_ocsp_ts()

% get_ocsp_ts
%==========================================================================
%
% USAGE:
%  ts = get_ocsp_ts()
%
% DESCRIPTION:
%  Function to retrieve subsurface temperature and salinity data from Ocean
%  Climate Station Papa mooring.
%
% INPUT:
%
%
% OUTPUT:
%
%  Struct containing subsurface temperature and salinity data
%
% AUTHOR:                      
%  September 17 2018. Zhihua Zheng                     [ zhihua@uw.edu ]   

%% old code

% t_stamp = datetime('now','TimeZone','UTC');
% t_stamp = datestr(t_stamp,'yyyy-mm-ddT'); % ISO8601 format time string
% t_stamp = [t_stamp,'00:00:00Z'];
% 
% url_head = 'http://osmc.noaa.gov/erddap/tabledap/';
% url_datasetID = 'papa_8bb5_dcd5_981f';
% url_fileType = 'mat';
% url_query = 'longitude,latitude,depth,time,TEMP,PSAL&time>=%s';
% 
% url_fmt = [url_head,url_datasetID,'.',url_fileType,'?',url_query];
% url = sprintf(url_fmt,t_stamp);
% ts = load(urlwrite('http://osmc.noaa.gov/erddap/tabledap/papa_8bb5_dcd5_981f.mat?longitude,latitude,depth,time,TEMP,SALT&time>=2018-09-27T00:00:00Z','test.mat'));

ts = struct();

%% ftp approach
ftpobj = ftp('data.ndbc.noaa.gov');
cd(ftpobj,'/data/oceansites/DATA_GRIDDED/PAPA');
fname = 'OS_PAPA_200706_M_TSVM_50N145W_dy.nc';
%fname = 'OS_PAPA_200706_M_TSVMBP_50N145W_hr.nc';
mget(ftpobj,fname);
close(ftpobj)

%% read variables

time = ncread(fname,'TIME');
t_ref = ncreadatt(fname,'TIME','units');
t_ref = t_ref(12:end);
t_ref = datenum(t_ref,'yyyy-mm-ddTHH:MM:SSZ');
time = time + t_ref; % days since t_ref
ts.time = time(end);
ts.date = string(datestr(ts.time,'yyyy-mm-dd HH:MM:SS'));

depth_t = double(ncread(fname,'DEPTH'));
ts.depth_t = depth_t(1:end-1); % eliminate the erroneous depth
t_prof = ncread(fname,'TEMP',[1 1 1 length(time)],[1 1 Inf 1]);
ts.t_prof = squeeze(squeeze(t_prof(1:end-1)));

depth_s = double(ncread(fname,'DEPPSAL'));
ts.depth_s = depth_s(1:end-1);
s_prof = ncread(fname,'PSAL',[1 1 1 length(time)],[1 1 Inf 1]);
ts.s_prof = squeeze(squeeze(s_prof(1:end-1)));

%% clean
delete(fname);

end

