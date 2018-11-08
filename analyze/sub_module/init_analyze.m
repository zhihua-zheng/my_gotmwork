%% init_analyze

% Subroutine to setup general configuration before analysis

% Zhihua Zheng, UW-APL, Sep. 17 2018

clear

%% gather data info using dialog box

prompt = {'GOTMRUN root directory path:','Specific output path:',...
    'Closure method used:','Time step (dt) used:',...
    'Saving period (nsave) used:'};
title = 'Specify Output Data Information';
dims = [1 100];

tmp = pwd;
run_dir = extractAfter(tmp,'run/');
clear tmp

definput = {'~/Documents/GitLab/GOTM_dev/run/',char(run_dir),...
    'SMCLT','60','60'};
data_info = inputdlg(prompt,title,dims,definput);

% path
cd ([data_info{1},data_info{2}])

% simulation info
turb_method = data_info{3};
dt = str2double(data_info{4});
nsave = str2double(data_info{5});


% check ./figs folder
if ~exist('figs/','dir') % if doesn't exist, create one
    mkdir figs
end

%% --------- general analyzing options ------------------------------------

mld_smooth = 1; % choose to smooth mixed layer depth or not


%% --------- load file ----------------------------------------------------
% find the netCDF file
dinfo = dir(fullfile('./*.nc'));
fname = fullfile('./',{dinfo.name});

% load
out = read_gotm_out(fname{:},2);

%% --------- read general variables ---------------------------------------

time = out.time;
% date = out.date;
% dateVec = datevec(char(date));

z = mean(out.z,2);
zi = mean(out.zi,2);
h = mean(out.h,2); % layer thickness

%% -------- general plotting specification info ---------------------------

if (time(end) - time(1)) < 4  % less than 4 days
    spec_info.timeformat = 'dd-hh';
elseif (time(end) - time(1)) <= 30 % less than a month
    spec_info.timeformat = 'mm-dd';
elseif (time(end) - time(1)) <= 366 % less than a year
    spec_info.timeformat = 'mm/yyyy';
else  % multiple years
    spec_info.timeformat = 'yyyy'; 
end
 