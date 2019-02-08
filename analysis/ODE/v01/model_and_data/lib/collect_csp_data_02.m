function collect_csp_data_02(cara_path, cara_repo_name, projects, peaklists, export_data, plot_results)
% rawList = importdata('shell_cat_all_txts.txt');

global FIG;    
global SB;

addpath('/Volumes/Data/yar/Dropbox/Programming/Matlab/tests');

global DATA_IVTNMR;

if exist('./lib','dir')
    datasave_folder = strcat(DATA_IVTNMR, '/'); % if running in root (runall)
else
    datasave_folder = strcat('../', DATA_IVTNMR, '/'); % if running from lib directory
end

%% Params to run "locally"
%===========================
if nargin == 0

	% Just to make sure we're in the right directory (was complaining)
    [fld,~,~] = fileparts( mfilename('fullpath') );
    cd(fld);
    
    reimport_data = 1;
    export_data = 0;

    use_only_7_peaks = 1;
    use_same_yaxis_scale_for_all_peaks = 0;
    
    plot_results = 1;
    
    %---------------
    cara_path = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR';
    % cara_repo_name = '170425_UP1_A1R1_R0_trace.cara';
    cara_repo_name = '170714_R1_UP1_SMN1_SMN2_trace.cara';

    projects = {...
        'IN60a'
        'IN61a'
        'IN62a'
        'IN65b'
        'IN70a'
        'IN71a'
        'IN72a'
        'IN75b'    
        };

    peaklists = {...
        '170425_IN60a'
        '170713_IN61a'
        '170713_IN62a'
        '170718_IN65b'
        '170425_IN70a'
        '170713_IN71a'
        '170713_IN72a'
        '170718_IN75b'    
        };
    
    % in collect_csp_data_02 - converted to MATRIX instead of CELL
    colors = [...
        [0.4, 0.4, 0.4]; % grey
        [0, 0.9, 0]; % light green
        [1, 0, 0]; % red
        [0.4, 0.4, 1]; % light blue
        [0, 0, 0]; % black
        [0, 0.7, 0]; % dark green
        [0.7, 0, 0]; % dark red
        [0, 0, 1]; % blue
        ];
    
    datasave_name = strcat(mfilename,'.mat');
    
else
    
    n_projects = numel(projects);
    
    % Colormap would make a new figure if none is open.
    % So preventively opening a figure, but setting it to invisible.    
    figure('Visible','off'); 
    colors = colormap(hsv(n_projects));
    
    reimport_data = 1;
        
    datasave_name = 'csp_all.mat';
    
    %%% Defined in the input!
    % export_data = x;    
    % plot_results = x;
    %%% Automatically assumed FALSE in the below code
    % use_only_7_peaks = x;
    % use_same_yaxis_scale_for_all_peaks = x;    
end

plotLW = 2;
plot_subset = 2; % 0 - all, 1 - A1R1, 2 - UP1       
n_projects = numel(projects);

%% Import data
%===========================
if ~reimport_data
    
    load(strcat(datasave_folder, datasave_name));

else % reimportData

    fl = cell(n_projects,1);
    
    for i=1:n_projects
        fl{i} = get_cara_peakshifts( fullfile(cara_path, cara_repo_name), projects{i}, peaklists{i});
        fl{i}.name = projects{i};
        fl{i}.n_peaks = numel(fl{i}.peak_tags);
%         fl{i}.color = colors{i}; % was in collect_csp_data
        fl{i}.color = colors(i,:); % now in collect_csp_data_02
        
    %     n_peaks = numel(fl.peak_tags);
    %     peakNumbers = 1:n_peaks;
    %     csArray = [peakNumbers; fl.HN_csp']; % transpose HN_csp to match earlier versions
    %     disp(projects{i});
    %     disp(peaklists{i});
    end
end % reimportData


%% Export data
%===========================
if export_data
% dlmwrite('../chemshift_results_hsqc.txt',csArray, '	');
% dlmwrite(output_file,csArray, '	');

% Just saving as MAT file - not table or anything.
save_mat_name = sprintf('%s%s', datasave_folder, datasave_name);
fprintf(1,'Data saved into %s\n', save_mat_name);        
save(strcat(save_mat_name),'fl');

% disp(csArray);
end


%% Select some data
%===========================
switch plot_subset
    case 1
        fl([5:end]) = [];
    case 2
        fl([1:4]) = [];        
end
    
if exist('use_only_7_peaks','var') && use_only_7_peaks
    fl{1}.H(8:end,:) = [];
    fl{1}.N(8:end,:) = [];
    fl{1}.peak_tags(8:end) = [];
    fl{1}.H_csp(8:end,:) = [];
    fl{1}.N_csp(8:end,:) = [];
    fl{1}.HN_csp(8:end,:) = [];
    fl{1}.n_peaks = 7;
end

%% Visualize
%===========================

if plot_results

n_projects = numel(fl);

% cellfun(@(x) x.name, fl, 'unif', 0)
n_peaks = cell2mat( cellfun(@(x) numel(x.peak_tags), fl, 'unif', 0) );

rows = 4;
columns = max(n_peaks);
sbSide = 150;
fSize = 11;
SB=0; % global defined at the top of the function
sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

%%% TODO: Make a helper function - so that with x.H - it returns ylim
mean_H = cellfun(@(x) mean(x.H,2), fl, 'unif', 0);
mean_N = cellfun(@(x) mean(x.N,2), fl, 'unif', 0);
mean_HN = cellfun(@(x) mean(x.HN_csp,2), fl, 'unif', 0);
% this max is for all peaks in ALL DATASETS!
max_delta_H = max(cell2mat( cellfun(@(x) abs(min(x.H,[],2)-max(x.H,[],2)), fl, 'unif', 0) ));
max_delta_N = max(cell2mat( cellfun(@(x) abs(min(x.N,[],2)-max(x.N,[],2)), fl, 'unif', 0) ));
max_delta_HN = max(cell2mat( cellfun(@(x) abs(min(x.HN_csp,[],2)-max(x.HN_csp,[],2)), fl, 'unif', 0) ));

ymin_H = cellfun(@(x) x-max_delta_H/2, mean_H, 'unif', 0);
ymax_H = cellfun(@(x) x+max_delta_H/2, mean_H, 'unif', 0);
ymin_N = cellfun(@(x) x-max_delta_N/2, mean_N, 'unif', 0);
ymax_N = cellfun(@(x) x+max_delta_N/2, mean_N, 'unif', 0);
ymin_HN = cellfun(@(x) x-max_delta_HN/2, mean_HN, 'unif', 0);
ymax_HN = cellfun(@(x) x+max_delta_HN/2, mean_HN, 'unif', 0);

% FIG global defined at the top of the function
if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
FIG=FIG+1;
set(figure(FIG), 'Color', repmat(0.95,1,3), 'Position', [0 300 columns*sbSide rows*sbSide]);

i = 1;

legend_string = {};
legend_string = {legend_string{:}, fl{i}.name};

for pk=1:n_peaks(i)
    d = fl{i};
    SB=pk; sb(SB);
    plot(d.time, d.H(pk,:), 'LineWidth', plotLW, 'Color', fl{i}.color);
    
    if n_projects > 1
        hold on;    
        for pp=2:n_projects
            if pk <= fl{pp}.n_peaks % some proj may have fewer peaks.
                plot(fl{pp}.time, fl{pp}.H(pk,:), 'LineWidth', plotLW, 'Color', fl{pp}.color);
            end
            legend_string = {legend_string{:}, fl{pp}.name};
        end
    end
    
    if pk == 1
        legend(legend_string);
    end
    
    axis tight;
    if exist('use_same_yaxis_scale_for_all_peaks','var') && use_same_yaxis_scale_for_all_peaks
        ylim([ymin_H{i}(pk) ymax_H{i}(pk)]);
    end
    title(d.peak_tags{pk});
end

% Plot N
for pk=1:n_peaks(i)
    d = fl{i};
    SB=columns+pk; sb(SB);
    plot(d.time, d.N(pk,:), 'LineWidth', plotLW, 'Color', fl{i}.color);
    
    if n_projects > 1
        hold on;    
        for pp=2:n_projects
            if pk <= fl{pp}.n_peaks % some proj may have fewer peaks.
                plot(fl{pp}.time, fl{pp}.N(pk,:), 'LineWidth', plotLW, 'Color', fl{pp}.color);
            end
            legend_string = {legend_string{:}, fl{pp}.name};
        end
    end    
    
    axis tight;
    if exist('use_same_yaxis_scale_for_all_peaks','var') && use_same_yaxis_scale_for_all_peaks
        ylim([ymin_N{i}(pk) ymax_N{i}(pk)]);
    end
    title(d.peak_tags{pk});
end

% Plot HN
for pk=1:n_peaks(i)
    d = fl{i};
    SB=columns*2+pk; sb(SB);
    plot(d.time, d.HN_csp(pk,:), 'LineWidth', plotLW, 'Color', fl{i}.color);

    if n_projects > 1
        hold on;    
        for pp=2:n_projects
            if pk <= fl{pp}.n_peaks % some proj may have fewer peaks.
                plot(fl{pp}.time, fl{pp}.HN_csp(pk,:), 'LineWidth', plotLW, 'Color', fl{pp}.color);
            end
            legend_string = {legend_string{:}, fl{pp}.name};
        end
    end    
    
    axis tight;
    if exist('use_same_yaxis_scale_for_all_peaks','var') && use_same_yaxis_scale_for_all_peaks
        ylim([0 ymax_HN{i}(pk)]);
    end
    title(d.peak_tags{pk});
end

disp('');

end % plot_results
end % main func
