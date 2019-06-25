function fig_handle = plot_traj4c(peak_traj, optns)
% Plots trajectories (rows) on individual sub-plots.
% Uses common y-max scale for all plots.
% Automatically decides on the optimal number of rows/columns.
% peak_traj HAS TO BE IN A CELL (even if only one is used!)

% Versions
% plot_traj3:
% - can overlay multiple sets of trajectories
% plot_traj3b:
% - shifted title
% - can display legend
% - has option to normalize intensity before for display
% plot_traj4:
% - can provide a different xaxis vector via optns.xaxis and .xaxis_label
% - can provide .yaxis_label
% plot_traj4b 
% - slightly modified version of plot_traj4_report
% - title is treated as CELL - so can plot multi-line titles (including CS data).
% plot_traj4c
% - optns.xaxis is now optns.xaxis_data

% TODO:
% - Modify ymax - so that the one MAX FOR ALL SETS is taken!
% - USE an OPTIONS STRUCTURE!

%% Test parameters (when executing script with no arguments)
%===========================================================
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
%     tmp = load(fullfile('datasave','assign_IN60a.mat'));
    tmp = load('/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/180410_PR1_SLH/analysis/datasave/assign_IN60a.mat');
%    peak_traj = load(fullfile('datasave','assign_UP1.mat')); % test UP1
    
%     peak_traj = peak_traj(:,1:30); % use only first 30 time-points
    peak_traj = cell(2,1);
    peak_traj{1} = tmp.assign.traj_HNcsp;
    peak_traj{2} = tmp.assign.traj_HNcsp;
    
    optns.plotLW = 2;
    optns.same_yscale = 0;
    optns.same_min_xscale = 1;
    optns.titles = arrayfun(@num2str, [1:size(peak_traj{1},1)]', 'unif', 0);
    optns.plotSym = '-';
    optns.markerSize = 4;
    
    optns.legend = 'bla';
end

%% Plotting setup
%===========================
assert(isa(peak_traj, 'cell'), 'peak_traj must be a cell array. Aborting.');
n_sets = numel(peak_traj);
n_plots = size(peak_traj{1},1);

if n_sets > 1 && numel(optns.plotSym) == 1
    optns.plotSym = repmat({optns.plotSym}, 1, n_sets); % will make a cell
elseif ~iscell(optns.plotSym)
    optns.plotSym = {optns.plotSym}; % make a single cell
end

if n_sets > 1 && numel(optns.plotLW) == 1
    optns.plotLW = repmat({optns.plotLW}, 1, n_sets); % will make a cell
elseif ~iscell(optns.plotLW)
    optns.plotLW = {optns.plotLW}; % make a single cell
end

if isfield(optns, 'colors') && ~isempty(optns.colors)
    colors = optns.colors;
elseif n_sets <= 4 % Four colors for RNA0-SMN1/2-EV2
    colors = [...
        [0.7, 0, 0]; % dark red
        [0, 0, 1]; % blue
        [0, 0.8, 0]; % dark green
        [0, 0, 0]; % black
        ];
else
    colors = colormap(lines(n_sets));
end

% rough calc number of squares
scrsize = get(0,'screensize');
scr_w = scrsize(3);
scr_h = scrsize(4);
tot_area = scr_w * scr_h;
area_per_plot = tot_area/n_plots;
sbSide = round(sqrt(area_per_plot));
columns = ceil(scr_w/sbSide);
rows = ceil(scr_h/sbSide);

sbWidth = sbSide*0.8;
sbHeight = sbSide*0.8;
fSize = 11;

global FIG;
% global SB;

if isempty(FIG), FIG=0; end; % if this is the first figure - set FIG index to zero
FIG=FIG+1;
fig_handle = figure(FIG);
set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 0 columns*sbWidth rows*sbHeight]);

SB=0;
sb = @(x) subplot(rows,columns,x,'FontSize',fSize); % Helper function to draw subplots

%% Normalize data
%===========================
if isfield(optns,'normalize_intensity') && ~isempty(optns.normalize_intensity)
    % - subtract the minimal value
    % - divide by max value
    peak_traj = cellfun(@(x) bsxfun(@minus, x, min(x,[],2)), peak_traj, 'un', 0 );
    peak_traj = cellfun(@(x) bsxfun(@rdivide, x, max(x,[],2)), peak_traj, 'un', 0 );
    
end

%% Plot data
%===========================

y_max = max( cellfun(@(x) max(max(abs(x),[],2)), peak_traj) ); % get max y-value to use common y-lim


if isfield(optns, 'xaxis') && ~isempty(optns.xaxis_data)
    xaxis_provided = 1;
    x_max = min( cell2mat( cellfun(@(x) max(x), optns.xaxis_data, 'unif', 0) ) );
else
    xaxis_provided = 0;
    x_max = min( cell2mat( cellfun(@(x) size(x, 2), peak_traj, 'unif', 0) ) );
end

if isfield(optns, 'xmax') && ~isempty(optns.xmax)
    x_max = optns.xmax;
end;

for i=1:n_plots
    hold on;
    SB=SB+1; handle_sb = sb(SB);
    if ~xaxis_provided
        plot(peak_traj{1}(i,:), optns.plotSym{1}, 'Color', colors(1,:), 'LineWidth', optns.plotLW{1}, 'MarkerSize', optns.markerSize);
    else
        plot(optns.xaxis_data{1}, peak_traj{1}(i,:), optns.plotSym{1}, 'Color', colors(1,:), 'LineWidth', optns.plotLW{1}, 'MarkerSize', optns.markerSize);
    end;
    if n_sets > 1   
        hold on;
        for iSet=2:n_sets
            if ~xaxis_provided
                plot(peak_traj{iSet}(i,:), optns.plotSym{iSet}, 'Color', colors(iSet,:), 'LineWidth', optns.plotLW{iSet}, 'MarkerSize', optns.markerSize);
            else
                plot(optns.xaxis_data{iSet}, peak_traj{iSet}(i,:), optns.plotSym{iSet}, 'Color', colors(iSet,:), 'LineWidth', optns.plotLW{iSet}, 'MarkerSize', optns.markerSize);
            end;            
        end
    end
    
    axis tight;
    if optns.same_yscale
        ylim([0 y_max]);
    end
    
    if optns.same_min_xscale
        xlim([0 x_max]);
    end
    
    if i <= numel(optns.titles)
%         title(optns.titles(i));
%         title(optns.titles{i}, 'horizontalAlignment', 'right', 'VerticalAlignment', 'top'); % if treating as cell - can plot multi-line
        title(optns.titles{i}, 'VerticalAlignment', 'top'); % if treating as cell - can plot multi-line
    else
        title('unassigned');
    end

    % remove XTicks in all    
    if isfield(optns, 'xaxis_label') && ~isempty(optns.xaxis_label) && i==n_plots
        xlabel(optns.xaxis_label);
    else
        set(gca,...
            'XTick', []...
        );
    end;
    
    first_col_ids = 1:columns:n_plots;
    % remove YTicks in all except first column    
    if i ~= first_col_ids
        set(gca,...
            'YTick', []...
        );
    end

    if i == first_col_ids(round(end/2))
        if isfield(optns, 'yaxis_label') && ~isempty(optns.yaxis_label)
            ylabel(optns.yaxis_label, 'Interpreter', 'none');
        end;
    end

    if i == n_plots && isfield(optns,'legend') && ~isempty(optns.legend)
        % to avoid shrinking the graph by legend:
        axP = get(gca,'Position'); 
        legend(optns.legend, 'Location', 'eo');
        set(gca, 'Position', axP);
    end
    
    %     'XTick', [2 4 6 8],...
    %     'YTick', [0 1],...
%         'XDir','reverse',...
%         'Ticklength', [0.002 0.002],...
%         'FontSize', (fSize-2)...

    % Maximize the size of subplots
    ax = get(handle_sb,'Position');
    ax(4) = ax(4)+ax(4)*0.13; % ax(4) - height
    ax(3) = ax(3)+ax(3)*0.2; % ax(3) - width
    set(handle_sb,'Position',ax);    
    
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    set(get(gca,'title'), 'Position', [x_lim(2)/2 (y_lim(2)-y_lim(1))*0.8]); % position - x & y coord
end

end

%============================================================
% function save_figure - use rather via external function!
%============================================================