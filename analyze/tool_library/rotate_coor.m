function new_vec = rotate_coor(out)

% rotate_coor
%==========================================================================
%
% USAGE:
%  new_vec = rotate_coor(out)

%
% DESCRIPTION:
%  Compute the coordinate of vectors under a new rotated reference system
%
% INPUT:
%
%  out - A struct containing all the model output from GOTM
%
% OUTPUT:
%
%  new_vec - A struct containing new coordinates for vectors under rotated
%            reference system
%
% AUTHOR:
%  October 8 2018. Zhihua Zheng                       [ zhihua@uw.edu ]
%


%% Read relevant variables

% The interpolation approach has been changed to use center_diff instead
u = out.u;
v = out.v;
u_stokes = out.u_stokes;
v_stokes = out.v_stokes;

%% Computation

new_vec.u = nan(size(u));
new_vec.v = nan(size(v));
new_vec.u_stokes = nan(size(u_stokes));
new_vec.v_stokes = nan(size(v_stokes));

for t = 1:length(out.time)
        
    w_angle = atan2(out.ty(t),out.tx(t)); % in rad unit (-pi, pi]
    
    % rotation matrix
    rot_mat = [ cos(w_angle), sin(w_angle);
               -sin(w_angle), cos(w_angle)];
           
    tmp = cat(2,u(:,t),v(:,t))';
    cur = rot_mat*tmp;
    new_vec.u(:,t) = cur(1,:);
    new_vec.v(:,t) = cur(2,:);

    tmp = cat(2,u_stokes(:,t),v_stokes(:,t))';
    cur = rot_mat*tmp;
    new_vec.u_stokes(:,t) = cur(1,:);
    new_vec.v_stokes(:,t) = cur(2,:);
end

end