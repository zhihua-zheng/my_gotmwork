function write_gotm_ini(Fname,A,time,z)
%
% write_gotm_ini
%==========================================================================
%
% USAGE:
%  write_gotm_ini(Fname,A,time,z)
%
% DESCRIPTION:
%  Writes variable into ASCII file as GOTM initial condition
%
% INPUT:
%
%  Fname - ASCII file name as a string, including path
%  fmtSpec - Format of the field, specified using formatting operators
%  A - Variable to be written (scalar, matrix or array)
%  time - Initial time as a string (yyyy-mm-dd HH:MM:SS)
%  z - 1D column vector for vertical coordinates of field [-,m]
%
% OUTPUT:
%
%  status - Error flag
%
% AUTHOR:
%  June 12 2019. Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================


fileID = fopen(Fname,'w');

m    = length(z);
Tfmt = ['%s  ',num2str(m),'  2\n'];

n = size(A,2);
switch n
    case 1    
        fmtSpec = '%8.2f   %9.6f\n';
    case 2
        fmtSpec = '%8.2f   %9.6f   %9.6f\n';
end

fprintf(fileID,Tfmt,time);
fprintf(fileID,fmtSpec,[z'; A']);

fclose(fileID);

end

