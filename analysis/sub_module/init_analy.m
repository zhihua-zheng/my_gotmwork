%% init_analy

% Subroutine to setup general configuration before analysis

% Zhihua Zheng, UW-APL, Sep. 17 2018

clear

%% --------- load file ----------------------------------------------------

[out,data_info] = load_gotm_out();

% check ./figs folder
if ~exist('figs/','dir') % if doesn't exist, create one
    mkdir figs
end

%% --------- general analyzing options ------------------------------------

mld_smooth = 0; % choose to smooth mixed layer depth or not

%% --------- read general variables ---------------------------------------

time = out.time;
dt = (time(2) - time(1))*24*3600; % time interval in output [s]
% date = out.date;
% dateVec = datevec(char(date));

% yearday from 'date2doy' (File Exchange)
yd = date2doy(time)-1; % day of year output from date2doy starts as 1

z = mean(out.z,2);
zi = mean(out.zi,2);
h = mean(out.h,2); % layer thickness

%% -------- general plotting specification info ---------------------------

if (time(end) - time(1)) < 4  % less than 4 days
    plot_info.timeformat = 'dd-hh';
elseif (time(end) - time(1)) <= 30 % less than a month
    plot_info.timeformat = 'dd';
elseif (time(end) - time(1)) <= 366 % less than a year
    plot_info.timeformat = 'mm';
else  % multiple years
    plot_info.timeformat = 'yyyy'; 
end
 