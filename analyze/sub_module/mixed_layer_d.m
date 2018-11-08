%% mixed_layer_d

% Subroutine to ananlyze themixed layer depth of the water column simualted

% Zhihua Zheng, UW-APL, Sep. 24 2018


%% ------------------------------------------------------------------------

mld = out.mld_surf;

% inertial frequency
f_Coriolis = gsw_f(out.lat); % [radian/s]
% inertial period
t_Coriolis = 2*pi/f_Coriolis/3600; % [hour]

if mld_smooth
    % filter length (~ 3 inertial periods)
    filter_length = 3*t_Coriolis;
    win_size = ceil(filter_length*3600/(nsave*dt));
    b = (1/win_size)*ones(1,win_size);
    a = 1;
    mld = filter(b,a,mld);
end


figure('position', [0, 0, 900, 300])
line(time,-mld,'LineWidth',.4,'Color',[.2 .6 .9])

spec_info.grid_on = 1;
spec_info.x_lim = [];
spec_info.y_lim = [];
spec_info.xlabel = 'time';
spec_info.ylabel = 'mixed layer depth (m)';
spec_info.lgd = [];
spec_info.save_path = './figs/mld';

line_annotate(spec_info)
