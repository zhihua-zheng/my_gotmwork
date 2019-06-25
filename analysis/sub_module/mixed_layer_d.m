%% mixed_layer_d

% Subroutine to ananlyze themixed layer depth of the water column simualted

% Zhihua Zheng, UW-APL, Sep. 24 2018


%% Filter -----------------------------------------------------------------

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

%% MLD line plot ----------------------------------------------------------

figure('position', [0, 0, 900, 300])
line(time_model,-mld,'LineWidth',.4,'Color',rgb('pinkish'))

plot_info.grid_on = 1;
plot_info.x_lim = [time_model(1) time_model(end)];
plot_info.y_lim = [];
plot_info.xlabel = 'time';
plot_info.ylabel = 'mixed layer depth (m)';
plot_info.lgd = [];
plot_info.save_path = './figs/mld';

line_annotate(plot_info)
