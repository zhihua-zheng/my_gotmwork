%% turb_compare

% Script to compare GOTM-develop's simulation results from different closure methods

% Zhihua Zheng, APL-UW, August 7 2018

%% load data

clear
main_dir = '~/Documents/GitLab/GOTM_dev/run/OCSPapa';
cd(main_dir)
dir_str = genpath(main_dir);

% regular expression for any substring doesn't contain path seperator
expression = ['[^',pathsep,']*'];
sep_idx = regexp(dir_str,pathsep); % index of path seperator

% truncate the string to ignore the main folder
sub_folders = regexp(dir_str(sep_idx(1):end),expression,'match');

turb_method = cell(size(sub_folders));
temp_all = cell(size(sub_folders));

% loop through subfolders
for i = 1:length(sub_folders)
    
    path = sub_folders{i};
    tmp = strsplit(path,'/');
    tmp = strsplit(tmp{end},'_');
    turb_method(i) = tmp(2); % get turbulence closure method name
    
    cd(path)
    dinfo = dir(fullfile('./*.nc'));
    fname = fullfile('./',{dinfo.name});
    
    % read variables
    temp_all{i} = squeeze(double(ncread(fname{:},'temp')));
    
    if i == 1
        sst_obs = squeeze(double(ncread(fname{:},'sst_obs')));
        z = mean(squeeze(double(ncread(fname{:},'z'))),2);
        time = squeeze(double(ncread(fname{:},'time')));
        
        % get the reference time
        t_ref = ncreadatt(fname,'time','units'); % get the attribute 'units' for 'time'
        t_ref = t_ref(15:end); % truncate to get the time string
        t_ref = datenum(t_ref, 'yyyy-mm-dd HH:MM:SS'); % datenumber for initialized time of simulation

        time = t_ref + time./3600/24;
        % date = string(datestr(time, 'yyyy/mm/dd HH:MM:SS'));
    end
end

clear t_ref fname dinfo
cd(main_dir)


%% KPP 
kpp_inx = find(contains(turb_method,'KPP'));

% set KPP-CVMix as bassis for comparison
base_inx = find(strcmp(turb_method,'KPP-CVMix'));

[T, Z] = meshgrid(time,z);
cnum = 15;

for i = kpp_inx
        
    figure('position', [0, 0, 900, 250])
    
    if i == base_inx     
        temp_plt = temp_all{base_inx}; % plot temp. for basis
        CL = [min(min(temp_plt)) max(max(temp_plt))];
        cmocean('matter')
        text_label = turb_method{base_inx};
    else
        temp_plt = temp_all{i} - temp_all{base_inx}; % plot temp. diff
        CL_limit = max(abs(min(min(temp_plt))),abs(max(max(temp_plt))));
        CL = [-CL_limit CL_limit];
        cmocean('balance')
        text_label = [turb_method{i},' --- ', turb_method{base_inx}];
    end
    
    conts = linspace(CL(1),CL(2),cnum);
    contourf(T,Z,temp_plt,conts,'LineWidth',0.01,'LineStyle','none')
    caxis(CL);
    box on
    datetick('x','mmm')
    text(time(200),z(15),text_label,'fontname','Optima','fontsize', 20)
    ylabel('depth ($$m$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
    setDateAxes(gca,'XLim',[time(1) time(end)],...
      'YLim',[-300 0],'fontsize',11,'fontname','computer modern','TickLabelInterpreter','latex')
    h = colorbar('EastOutside');
    h.Label.String = 'potential temperature ($$^{\circ}C$$)';
    h.Label.Interpreter = 'latex';
    h.Label.FontName = 'computer modern';
    h.Label.FontSize = 14;
    set(h,'TickLabelInterpreter','latex','fontsize',9);
    fig_name = ['./temp_diff_',turb_method{i}];
    set(gca(),'LooseInset', get(gca(),'TightInset')); % no blank edge
    saveas(gcf, fig_name, 'epsc');
    
%    export_fig (fig_name,'-eps','-transparent','-painters')
    
end

%% SMC
smc_inx = find(contains(turb_method,{'SMC','K-EPSILON'}));

% set basic SMC as basis for comparison
base_inx = find(strcmp(turb_method,'SMC'));

for i = smc_inx
    
    figure('position', [0, 0, 900, 250])

    if i == base_inx    
        temp_plt = temp_all{base_inx}; % plot temp. for basis
        CL = [min(min(temp_plt)) max(max(temp_plt))];
        cmocean('matter')
        text_label = turb_method{base_inx};
    else
        temp_plt = temp_all{i} - temp_all{base_inx}; % plot temp. diff
        CL_limit = max(abs(min(min(temp_plt))),abs(max(max(temp_plt))));
        CL = [-CL_limit CL_limit];
        cmocean('balance')
        text_label = [turb_method{i},' --- ', turb_method{base_inx}];
    end
    
    conts = linspace(CL(1),CL(2),cnum);
    contourf(T,Z,temp_plt,conts,'LineWidth',0.01,'LineStyle','none')
    caxis(CL);
    box on
    datetick('x','mmm')
    text(time(200),z(15),text_label,'fontname','Optima','fontsize', 20)
    ylabel('depth ($$m$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
    setDateAxes(gca,'XLim',[time(1) time(end)],...
      'YLim',[-300 0],'fontsize',11,'fontname','computer modern','TickLabelInterpreter','latex')
    h = colorbar('EastOutside');
    h.Label.String = 'potential temperature ($$^{\circ}C$$)';
    h.Label.Interpreter = 'latex';
    h.Label.FontName = 'computer modern';
    h.Label.FontSize = 14;
    set(h,'TickLabelInterpreter','latex','fontsize',9);
    fig_name = ['./temp_diff_',turb_method{i}];
    set(gca(), 'LooseInset', get(gca(), 'TightInset')); % no blank edge
    saveas(gcf, fig_name, 'epsc');

%    export_fig (fig_name,'-eps','-transparent','-painters')
    
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     fig_pos = fig.PaperPosition;
%     fig.PaperSize = [fig_pos(3) fig_pos(4)];
%     print(fig,fig_name,'-depsc')
end


%% ePBL
epbl_inx = find(contains(turb_method,'EPBL'));

% set EPBL as basis for comparison
base_inx = find(strcmp(turb_method,'EPBL'));

for i = epbl_inx
        
    figure('position', [0, 0, 900, 250])
    
    if i == base_inx   
        temp_plt = temp_all{base_inx}; % plot temp. for basis
        CL = [min(min(temp_plt)) max(max(temp_plt))];
        cmocean('matter')
        text_label = turb_method{base_inx};
    else
        temp_plt = temp_all{i} - temp_all{base_inx}; % plot temp. diff
        CL_limit = max(abs(min(min(temp_plt))),abs(max(max(temp_plt))));
        CL = [-CL_limit CL_limit];
        cmocean('balance')
        text_label = [turb_method{i},' --- ', turb_method{base_inx}];
    end
    
    conts = linspace(CL(1),CL(2),cnum);
    contourf(T,Z,temp_plt,conts,'LineWidth',0.01,'LineStyle','none')
    caxis(CL);
    box on
    datetick('x','mmm')
    text(time(200),z(15),text_label,'fontname','Optima','fontsize', 20)
    ylabel('depth ($$m$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
    setDateAxes(gca,'XLim',[time(1) time(end)],...
      'YLim',[-300 0],'fontsize',11,'fontname','computer modern','TickLabelInterpreter','latex')
    h = colorbar('EastOutside');
    h.Label.String = 'potential temperature ($$^{\circ}C$$)';
    h.Label.Interpreter = 'latex';
    h.Label.FontName = 'computer modern';
    h.Label.FontSize = 14;
    set(h,'TickLabelInterpreter','latex','fontsize',9);
    fig_name = ['./temp_diff_',turb_method{i}];
    set(gca(), 'LooseInset', get(gca(), 'TightInset')); % no blank edge
    saveas(gcf, fig_name, 'epsc');
    
    %export_fig (fig_name,'-eps','-transparent','-painters')
end

%% SST

cmp = distinguishable_colors(14); % distinguishable color schemes

figure('position', [0, 0, 1000, 500])

for i = 1:length(turb_method)
    
    temp = temp_all{i};
    sst = temp(end,:);
    line(time(2:end),sst(2:end),'LineWidth',.8,'Color',cmp(i,:))
end
line(time(2:end),sst_obs(2:end),'LineWidth',.8,'Color',cmp(i+1,:))

box on
datetick('x','mmm')
lgd = legend(['observation',turb_method],'Location','best');
set(lgd,'Interpreter','latex','fontsize', 14)
ylabel('sea surface temperature ($$^{\circ}C$$)', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
xlabel('time', 'fontname', 'computer modern', 'fontsize', 14,'Interpreter', 'latex')
setDateAxes(gca,'XLim',[time(2) time(end)],...
  'fontsize',11,'fontname','computer modern','TickLabelInterpreter', 'latex')

set(gca(), 'LooseInset', get(gca(), 'TightInset')); % no blank edge
saveas(gcf, './sst', 'epsc');


