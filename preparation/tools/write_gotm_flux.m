function write_gotm_flux(Fname,A,time)
%
% write_gotm_flux
%==========================================================================
%
% USAGE:
%  write_gotm_flux(Fname,A,time)
%
% DESCRIPTION:
%  Writes variable into ASCII file as GOTM boundary condition (surface flux)
%
% INPUT:
%
%  Fname - ASCII file name as a string, including path
%  fmtSpec - Format of the field, specified using formatting operators
%  A - (t,m) variable to be written
%  time - 1D column vector, timestamp as a string (yyyy-mm-dd HH:MM:SS)
%
% OUTPUT:
%
%  status - Error flag
%
% AUTHOR:
%  June 12 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================


fileID = fopen(Fname,'w');

n = size(A,2);
switch n
    case 1    
        fmtSpec = '%s  % 8.6e\n';
    case 2
        fmtSpec = '%s  % 8.6e  % 8.6e\n';
end

% H = [cellstr(time) num2cell(A)];

fprintf(fileID,fmtSpec,[time'; A']);

fclose(fileID);

end

