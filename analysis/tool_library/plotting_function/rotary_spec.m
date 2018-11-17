function rotary_spec(f, p_cur, f_in, save_switch)


% rotary_spec
%==========================================================================
%
% USAGE:
%  rotary_spec(cur, save_switch)
%
% DESCRIPTION:
%  Function to generate a plot of velocity power spectrum
%
% INPUT:
%
%  f - the frequency domain from pwelch output, in cycle per day
%  p_cur - the power spectrum density esitmate from pwelch output
%  f_in - local inertial frequency, in cycle per day
%  save_switch - 0 or 1 to specify if want the figure 
%
% OUTPUT:
%
%  figure of velocity power spectrum with detailed annotation
%
% AUTHOR:
%  September 5 2018. Zhihua Zheng                       [ zhihua@uw.edu ]


n = length(f);

figure('position', [0, 0, 780, 350])

% counter-clockwise, negative
loglog(f(1:n/2),p_cur(1:n/2),'Color',[.1 .6 .7],'LineStyle','-','LineWidth',.3) 
hold on

% clockwise, positive
loglog(f(1:n/2),flip(p_cur(n/2+1:end)),'Color',[.9 .4 .6],'LineStyle','-','LineWidth',.3) 
hold on

% mark inertial frequency
loglog([f_in f_in],[.001*min(p_cur) 100*max(p_cur)],...
    'Color',[.3 .4 .5],'LineStyle','-.','LineWidth',.3);
hold off


% call 'line_annotate.m' to add figure details
spec_info.grid_on = 1;
spec_info.xlabel = 'frequency (cycle per day)';
spec_info.ylabel = 'PSD ($$(m/s)^2/cpd$$)';
spec_info.x_lim = [f(1) f(n/2)];
spec_info.y_lim = [];
spec_info.lgd = {'counter clockwise', 'clockwise'};
spec_info.lgd_pos = [];
if save_switch
    spec_info.save_path = './figs/cur_spec';
else
    spec_info.save_path = [];
end

line_annotate(spec_info)

end

