function fig_handle = plot_traj3b(peak_traj, optns)
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

% TODO:
% - Modify ymax - so that the one MAX FOR ALL SETS is taken!
% - USE an OPTIONS STRUCTURE!

%% Test parameters (when executing script with no arguments)
%===========================================================
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    tmp = load(fullfile('datasave','assign_IN60a.mat'));
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

if n_sets <= 4 % Four colors for RNA0-SMN1/2-EV2
    colors = [...
        [0.7, 0, 0]; % dark red
        [0, 0, 1]; % blue
        [0, 0.8, 0]; % dark green
        [0, 0, 0]; % black
        ];
elseif isfield(optns, 'colors') && ~isempty(optns.colors)
    colors = optns.colors;
else
    colors = colormap(jet(n_sets));
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
x_max = min( cell2mat( cellfun(@(x) size(x, 2), peak_traj, 'unif', 0) ) );

for i=1:n_plots
    SB=SB+1; handle_sb = sb(SB);
    plot(peak_traj{1}(i,:), optns.plotSym, 'Color', colors(1,:), 'LineWidth', optns.plotLW, 'MarkerSize', optns.markerSize);
    if n_sets > 1
        hold on;
        for iSet=2:n_sets
            plot(peak_traj{iSet}(i,:), optns.plotSym, 'Color', colors(iSet,:), 'LineWidth', optns.plotLW, 'MarkerSize', optns.markerSize);
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
        title(optns.titles(i));
    else
        title('unassigned');
    end

    % remove XTicks in all
    set(gca,...
        'XTick', []...
    );

    first_col_ids = 1:columns:n_plots;
    % remove YTicks in all except first column    
    if i ~= first_col_ids
        set(gca,...
            'YTick', []...
        );
    end

    if i == first_col_ids(round(end/2))
        ylabel('HN_CSP [ppm]', 'Interpreter', 'none');
    end

    if i == n_plots && isfield(optns,'legend') && ~isempty(optns.legend)
        legend(optns.legend, 'Location', 'eo');
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