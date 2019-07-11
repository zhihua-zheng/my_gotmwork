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
%  A - 3D matrix (z,m,t), variable with dimension m to be written
%  time - 1D vector, timestamp as a string (yyyy-mm-dd HH:MM:SS)
%  z - 1D column vector, vertical coordinates (descending) of field [-,m]
%
% OUTPUT:
%
%  status - Error flag
%
% AUTHOR:
%  June 12 2019. Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================


fileID = fopen(Fname,'w');

iz = size(A,1);
im = size(A,2);
it = size(A,3);

Tfmt  = ['%s  ',num2str(iz),'  ',num2str(im+1),'\n'];
ifmt  = '   %9.6f';
zfmt  = '%8.2f';

VIfmt = repmat(ifmt,1,im);
Vfmt  = strcat(zfmt,VIfmt,'\n');

% loop through time-points
if it > 1

  for k = 1:it
    
      fprintf(fileID,Tfmt,time(k));
      fprintf(fileID,Vfmt,[z'; A(:,:,k)']);
  end

else

  fprintf(fileID,Tfmt,time);
  fprintf(fileID,Vfmt,[z'; A']);
end

fclose(fileID);

end

