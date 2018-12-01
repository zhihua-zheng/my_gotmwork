function [out, data_info] = load_gotm_out()

% load_gotm_out
%==========================================================================
%
% USAGE:
%  [out, data_info] = load_gotm_out()
%
% DESCRIPTION:
%  Load GOTM simulation output from specific run directory, returned as a
%  cell containing structs for different output, if multiple data files are
%  present.
%
% INPUT:
%
% OUTPUT:
%
%  out - cell array for multiple data files, struct for single file.
%  data_info - cell array containing paths and turbulence method.
%
% AUTHOR:
%  November 30 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%

%% --------- gather data info using dialog box ----------------------------

prompt = {'GOTMRUN root directory path:','Specific output path:',...
    'Closure method used:'};
title = 'Specify Output Data Information';
dims = [1 100];

tmp = pwd;
run_dir = extractAfter(tmp,'run/');
tmp = strsplit(run_dir,{'_','/'});

% get the turb_method
for i = 1:length(tmp)
    [~,status] = str2num(char(tmp(i))); %#ok<ST2NM>
    if status
        method = tmp(i-1);
    end
end

definput = {'~/Documents/GitLab/GOTM_dev/run/',char(run_dir),...
    char(method)};
data_info = inputdlg(prompt,title,dims,definput);

%% --------- load file ----------------------------------------------------

% find the netCDF file
dinfo = dir(fullfile('./*.nc'));
fname = fullfile('./',{dinfo.name});

% load
if length(fname) < 2
    out = read_gotm_out(fname{:},2);
else
    out = cell(length(fname),1);
    
    for i = 1:length(fname) 
        out{i} = read_gotm_out(fname{i},2);
    end
end

end
