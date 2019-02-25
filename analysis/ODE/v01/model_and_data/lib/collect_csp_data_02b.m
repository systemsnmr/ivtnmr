function collect_csp_data_02b(cara_path, cara_repo_name, projects, peaklists, export_data, plot_results)
% rawList = importdata('shell_cat_all_txts.txt');
% collect_csp_data_02b - when adding import of RNA0, SMN2, EV2 replicates

global FIG;    
global SB;

addpath('/Volumes/Data/yar/Dropbox/Programming/Matlab/tests');

global DATA_IVTNMR NMR_DATA_PATH;

datasave_folder = DATA_IVTNMR;

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
%     cara_repo_name = '180820_170714_R1_UP1_SMN1_SMN2_trace.cara';
    cara_repo_name = '190121_UP1_replicates_f_paper.cara';

    projects = {...
        'IN70a'
        'IN70b'
        'IN70c'
        'IN71a'
        'IN71b'
        'IN71c'
        'IN72a'
        'IN72b'
        'IN72c'
        'IN72d'
        'IN75a'
        'IN75b'
        'IN75h'
        'IN75z'
        };

    peaklists = {...
        '170425_IN70a'
        '180820_IN70b'
        'IN70c'
        '170713_IN71a'
        'IN71b'
        'IN71c'
        '170713_IN72a'
        '180820_IN72b'
        'IN72c'
        'IN72d'
        '180820_IN75a'    
        '170718_IN75b'    
        'IN75h'
        'IN75z'
        };
    
%     select_subset = [4 5];
    
    % in collect_csp_data_02 - converted to MATRIX instead of CELL
    colors = [...
        [0.5, 0.5, 0.5]; % grey
        [0.25, 0.25, 0.25]; % grey2
        [0, 0, 0]; % black
        
        [0, 0.9, 0]; % light green
        [0, 0.75, 0]; % med green
        [0, 0.6, 0]; % dark green
        
        [1, 0, 0]; % red
        [0.85, 0, 0]; % med red
        [0.7, 0, 0]; % dark red
        [0.55, 0, 0]; % very dark red
        
        [0.45, 0.45, 1]; % light blue
        [0.3, 0.3, 1]; % blue
        [0.15, 0.15, 1]; % blue
        [0, 0, 1]; % blue
        ];


    if exist('select_subset','var') && ~isempty(select_subset)
        projects = projects(select_subset);
        peaklists = peaklists(select_subset);
        colors = colors(select_subset, :);
    end;

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

plotLW = 1;
plot_subset = 0; % 0 - all, 1 - A1R1, 2 - UP1       
n_projects = numel(projects);

datafile_path = fullfile(datasave_folder, datasave_name);

%% Import data
%===========================
if ~reimport_data
    
    load( datafile_path );

else % reimportData

    fl = cell(n_projects,1);
    
    for i=1:n_projects
        get_cara_peakshifts_optns.nmr_data_path = NMR_DATA_PATH;
        fl{i} = get_cara_peakshifts_04( fullfile(cara_path, cara_repo_name), projects{i}, peaklists{i}, get_cara_peakshifts_optns);
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
save_mat_name = datafile_path;
fprintf(1,'Data saved into %s\n', save_mat_name);        
save(save_mat_name, 'fl');

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
sbSide = 250;
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
%             legend_string = [legend_string, fl{pp}.name];
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
        end
    end    
    
    axis tight;
    if exist('use_same_yaxis_scale_for_all_peaks','var') && use_same_yaxis_scale_for_all_peaks
        ylim([ymin_N{i}(pk) ymax_N{i}(pk)]);
    end
    title(d.peak_tags{pk});
end

% % plot CSP of the last peaks - 20190204 - disabled
% cell2mat(cellfun(@(x) max(x.HN_csp([6 7],:),[],2), fl, 'un', 0))

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
