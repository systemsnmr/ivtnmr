function a_visualize
% Visualize / analyse.
% This just plots constants - removed plotting of data.
clear; close all;

%% Select results file
%=====================
%%% By default takes first results file (assumes its only one)
model_sub_folder = 'model_and_data';
result_files = dir( fullfile(model_sub_folder, 'fit_results', '*.mat') );
datafile = fullfile(model_sub_folder, 'fit_results', result_files(1).name);
% datafile = fullfile(model_sub_folder, 'results', 'a_v01_multi.mat');

%% Settings / flags
%=====================
global FIG;

% select_sets = {'IN72b' 'IN75a'}; % can select which datasets to plot

f_show_data_points = 1;

sbSide = 180;
fSize = 11;    
fitLW = 1;

% maps dataset name/ids to actual RNA names
dset_to_rna = containers.Map({'IN70', 'IN71', 'IN72', 'IN75'},...
    {'RNA0','SMN1','SMN2','EV2'});

k_names = {'KeqS' 'kcatnorm' 'kdephosNTP' 'KdPR'};
k_units = {'Keq, MgHPO4 [mM]' 'kcat [nt s-1]' 'kdephosph [s-1]' 'KD [uM]'};
n_constants = numel(k_names);

addpath(genpath('./lib'))

%% Load
%=====================
load(datafile);
n_sets = numel(fit_data);

% select only subset
if exist('select_sets','var') && ~isempty(select_sets)
    [~, sel_idx] = ismember(select_sets, cellfun(@(x) x.name(1:5), fit_data(:,1), 'un', 0));
    sel_idx = sel_idx(sel_idx>0);
else
    sel_idx = 1:n_sets;
end
fit_data = fit_data(sel_idx,:);

% get names of datasets, and RNA names
dset_names = cellfun(@(x) x.name(1:5), fit_data(:,1), 'un', 0);
rna_names = cellfun(@(x) dset_to_rna(x(1:4)), dset_names, 'un', 0);

rnas = sort(unique(rna_names));
% if want to fix the order:
% rnas = {'RNA0' 'SMN1' 'SMN2' 'EV2'};
n_rnas = numel(rnas);

k_mean_table = nan(n_constants, n_rnas);
k_sd_table = nan(n_constants, n_rnas);

if f_show_data_points
    dset_numbers = cellfun(@(x) str2num(x(3:4)), dset_names);
    [unique_counts,~] = hist(dset_numbers, unique(dset_numbers));
    n_max_replicates = max(unique_counts);

    % cells for each constant. inside each cell an array with max number of
    % replicates.
    k_all_vals = repmat({nan(n_rnas,n_max_replicates)}, n_constants, 1); % need just one vector per Constant to display!
end


%%% Generate means and std for same RNA
%=======================================
for iRNA=1:n_rnas
    subset_idx = strcmp(rnas(iRNA), rna_names);
    
    for iConst=1:n_constants
        kval = cellfun(@(x) x.( k_names{iConst} ), fit_data(subset_idx));
                
        if f_show_data_points
%             fprintf(1,'iRNA=%i, iConst=%i\n', iRNA, iConst);
            % pad kval vector with NaNs if needed
            kval = [kval; nan(n_max_replicates-numel(kval),1)];
            k_all_vals{iConst}(iRNA,:) = kval; % fill one row for each RNA
        end
        
        k_mean_table(iConst,iRNA) = nanmean(kval);
        k_sd_table(iConst,iRNA) = nanstd(kval);        
    end;
end

% k_mean_table
% k_sd_table

%% Visualize
%=====================
rows = 2; % using 2 rows by default to get enough white space for labels
columns = n_constants; 

if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
FIG=FIG+1;
fig_to_save = figure(FIG); set(figure(FIG), 'Color', ones(1,3), 'Position', [0 300 columns*sbSide*1.2 rows*sbSide]);

SB=0;
sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

% can use extra rows to plot some other stuff.
% for iRow = 1:rows
    for iConst = 1:n_constants
        SB=SB+1;
        sb(SB);        
        
        errorbar(k_mean_table(iConst,:), k_sd_table(iConst,:), 's', 'LineWidth', fitLW);

        %%% Plot individual data points
        if f_show_data_points
            hold on;
            plot(k_all_vals{iConst}, 'ob');
        end; % f_show_data_points
        
        
        axis tight;
        yax = get(gca,'YLim');
        ylim([0 max(yax)*1.3]);
        ylabel( k_units{iConst} );
        
        % log scale for Kds
        if strcmp(k_names{iConst},'KdPR')
%             ylim([0 4000]);
            set(gca,'YTick',[1 10 100 1000]);
        end
        
        set(gca,'xtick',1:n_rnas);
        set(gca,'xticklabel',rnas);
        
        title(k_names{iConst});
        
    end; clear iConst;
% end; clear iRow;
    
file_name = fullfile(model_sub_folder, 'figure_images', sprintf('k_%s.pdf', sprintf('%s_', dset_names{:})));
save_figure(fig_to_save, file_name);


end % main func

function save_figure(fig_handle,file_name)
%% Save figure as pdf
%=============================================================
    set(fig_handle,'Units','Inches');
    pos = get(fig_handle,'Position');
    set(fig_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_handle,file_name,'-dpdf','-r1200')
%     print(fig_handle,fullfile(file_name),'-dpdf','-r1200')

% CONVENIENT TO USE with mfilename variable - which returns the name of the current file!
% 
% -r option sets the resolution (also for the VECTOR IMAGES!) -- set to 300-1200 if contour approximation is bad.
% -r0
% -r300
% -r1200
end