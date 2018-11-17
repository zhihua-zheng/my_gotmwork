function hodogram_t(t, cur, save_switch)

% hodogram_t
%==========================================================================
%
% USAGE:
%  hodogram_t(t, cur, save_switch)
%
% DESCRIPTION:
%  Function to generate a temporal hodogram for input velocity data
%
% INPUT:
%
%  t - one dimensional array of date number
%  cur - one dimensional array of velocity in complex form
%  save_switch - 0 or 1 to specify if want the figure 
%
% OUTPUT:
%
%  Temporal evolution of velocity hodogram with detailed annotation
%
% AUTHOR:
%  September 6 2018. Zhihua Zheng                       [ zhihua@uw.edu ]


t_int = t(t==floor(t)); % label only integer days

color_mat = distinguishable_colors(length(t_int)); 
% define discrete colormap

figure('position', [0, 0, 480, 480])

u_m = 1.2*max(abs(real(cur)));
v_m = 1.2*max(abs(imag(cur)));

colormap(color_mat);
scatter(real(cur),imag(cur),30,t,'filled')
colorbar('Northoutside')
h = cbdate(t_int,'mmm-dd','horiz'); 
hold on
scatter(0,0,55,'k','p')
hold on
plot(cur,'Color',[.5 .6 .7],'LineStyle',':','LineWidth',.1)
hold off 

  pbaspect([1 1 1]) % same scale for x, y
  h.Label.Interpreter = 'latex';
  h.Label.FontName = 'computer modern';
  h.Label.FontSize = 14;
  set(h,'TickLabelInterpreter','latex','fontsize',9);
    
spec_info.grid_on = 0;
spec_info.xlabel = 'u (m/s)';
spec_info.ylabel = 'v (m/s)';
spec_info.x_lim = [-u_m u_m];
spec_info.y_lim = [-v_m v_m];
spec_info.lgd = [];
if save_switch
    spec_info.save_path = './figs/hodo';
else
    spec_info.save_path = [];
end

line_annotate(spec_info)

end

