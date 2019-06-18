function viz_prot_perturb3b()

% TODO:
% - Automatically? find out which expnos used for the mapper?!? (Is this possible?)

% Versions:
% viz_prot_perturb:   initial version, based on 
%   Projects/2015_NMR_peak_tracking/site_mapper/190411_TDP_FUS_SMN2_LLPS/analysis/runall_TDP_190411.m

% viz_prot_perturb2:
%   - includes 31P data import and plotting against it
%   - partially based on CUG / MBNL analysis /Volumes/Data/yar/Dropbox/_eth2/project_CUG_MBNL/code/180926_mapper_with_31P/analysis/runall_var_MBNL_31P_v180930_f_report.m
%   - option to plot both CSP and Intensity
%   - option to plot chemical shifts of peaks
%   - fixed opt_to_duplicate

% viz_prot_perturb3:
%   - simultaneous plot of several mapper sets - with CSP and traj
%   - assumes xaxis has to be defined for each set (i.e. not just duplicating!)
%   - uses folder-namespace for '+lib'; (instead of genpath...)
%   - Testing with SLH datasets.

% viz_prot_perturb3b:
%   - Moved to github repo - made work with +lib and utils/rbnmr.

%% Settings
%=================================
tic;
addpath(genpath(fullfile(pwd,'utils')));

analyse_PTBRRM1_SLH = 1;

% Sections to run:
flag_run_assignments = 1;
flag_load_assignments = 1;

f_run_kd_fits = 1;

% Switching Y-data
% 1,2 - will overlay data from all datasets.
% 3 (CSP+intens) - forces plotting only of the first dataset!
optns.yaxis = 1; % 1-CSP; 2-intensity; 3-both

% Switching X-Axis data
xaxis = 1; % 1 - spectrum id; 2 - time; 3 - 31P(RNA)

%%% Select subsets (if have many defined, but want to plot only a few):
% select_dsets = [1 3 5]; % 1:2 4 - compares 4th nucleotide

% select_residues = [1:57 59:61];
select_residues = [1:3];

% "select_spectra" - in reverse - QnD solution to fix nc_proc
% drop_timepoints = [2 3 6 22:35]; % tmp hack - to match MAPPER vector to P50N525 time

optns.colors = [...
        [0, 0, 0]; % black
        [0.7, 0, 0]; % dark red
        [0, 0.6, 0]; % dark green
        [1, 0, 0]; % red
        [0, 0.9, 0]; % green        
        ];

if any(optns.yaxis == [2 3]) % If plotting intensity or both
    optns.normalize_intensity = 1; 
end;
optns.same_yscale = 1; % scales INTENSITY!?!
    
%============================================================
%%% Inputs: mapper data, peaklists, RNA integration, protein sequences
%============================================================
current_dir = fileparts(mfilename('fullpath'));

if ~analyse_PTBRRM1_SLH
    nmr_dataset_dirs = {
        '/Volumes/Data/yar/Dropbox/_eth2/project_Measles/NMR/spectra/190415_IN115a_pA20_co-P50N525_100uM_298K_600'
        '/Volumes/Data/yar/Dropbox/_eth2/project_Measles/NMR/spectra/190415_IN115a_pA20_co-P50N525_100uM_298K_600'
        };

    mapper_data_dirs = {
    %     fullfile(fileparts(mfilename('fullpath')), '..', 'mapper_data')
        fullfile(nmr_dataset_dirs{1}, 'IN115a_mapper')
        fullfile(nmr_dataset_dirs{1}, 'IN115a_mapper')
        };
    % s3_tracking dir - also contains intensities
    mapper_data_paths = cellfun(@(x) fullfile(x,'s3_tracking','trajectories_filtered.csv'), mapper_data_dirs, 'un', 0);

    % Can provide specific expnos, if not all experiments used in tracking!
    % This is primarily needed if need to read TIME or RNA concentration - for xaxis.
    % If note provided - script tries to read ALL expnos from the nmr_dataset_dir.
    mapper_expnos = [4003 4008 4013 4018 4023 4028 4033 4038 4043 4078 4083 4088];

    %%% Peaklists. If different for each set - provide multiple entries here.
    % If same for all sets - leave just one.
    %============
    peaklist_dirs = {fullfile( nmr_dataset_dirs{1}, 'peaklists' )};
    peaklist_names = {'peaks_start_mapper.peaks'};
    peaklist_paths = cellfun(@(x,y) fullfile(x,y), peaklist_dirs, peaklist_names, 'un', 0);

    %%% RNA data (optional, if not trying to plot time/RNA on xaxis)
    %============
    RNA_data_dirs = {...
        'datasave'
        'datasave'
        };
    RNA_data_filenames = {...
        'RNA_A20_P50N525.csv'
        'RNA_A20_P50N525.csv'
        };
    RNA_data_paths = cellfun(@(x,y) fullfile(x,y), RNA_data_dirs, RNA_data_filenames, 'un', 0);

end % ~analyse_PTBRRM1_SLH

if analyse_PTBRRM1_SLH
    nmr_dataset_dirs = {
	'/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra/180308_IN93a_SLH_co-NPR1_303K_600'
	'/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra/180315_IN95a_SLHar3_co-NPR1_303K_600'
    };
               
    mapper_data_dirs = {
    '/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/190321_PR1_SLH_C11/IN93a'
    '/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/190321_PR1_SLH_C11/IN95a'        
        };
    
    % s3_tracking dir - also contains intensities
    mapper_data_paths = cellfun(@(x) fullfile(x,'s3_tracking','trajectories_filtered.csv'), mapper_data_dirs, 'un', 0);

    % can define one if same for all
    peaklist_paths = {
        '/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/190321_PR1_SLH_C11/analysis/data_test/PR1_peaklist.peaks'
        };    
    
    %%% RNA data (optional, if not trying to plot time/RNA on xaxis)
    %============
    RNA_data_dirs = {...
        'datasave'
        'datasave'
        };
    RNA_data_filenames = {...
        'RNA_IN93a.csv'
        'RNA_IN95a.csv'
        };
    RNA_data_paths = cellfun(@(x,y) fullfile(x,y), RNA_data_dirs, RNA_data_filenames, 'un', 0);
    
end

protein_conc = 0.15; % uM

%%
%================================================================
%%%% EDITING BELOW THIS LINE IS NORMALLY NOT REQUIRED
%================================================================
n_sets = numel(mapper_data_paths);

flag_show_peak_CS = 1;

%% Input preprocessing
%=====================
% takes the second underscore-separated element as the ID of expt.
[~,dset_list] = cellfun(@(x) fileparts(x), nmr_dataset_dirs, 'un', 0);
dset_split_arr = cellfun(@(x) regexp(x, '_', 'split'), dset_list, 'un', 0);    
dset_id = cellfun(@(x) x{2}, dset_split_arr, 'un', 0);
dset_id_long = cellfun(@(x) sprintf('%s_%s_%s', x{2}, x{3}, x{4}), dset_split_arr, 'un', 0);
clear dset_split_arr;

dset_names = dset_id;

% Try getting expnos (if NMR paths are provided)
if exist('nmr_dataset_dirs','var') && ~isempty(nmr_dataset_dirs);
    nmr_dataset_dirs_exist = all(cell2mat(cellfun(@(x) isdir(x), nmr_dataset_dirs, 'un', 0)));
else
    nmr_dataset_dirs_exist = 0;
end;
if nmr_dataset_dirs_exist
    dsets_time0 = cellfun(@(x) lib.getTime0b(x), nmr_dataset_dirs, 'un', 0);
    dset_exp_31P = lib.getNMRExpnos(nmr_dataset_dirs, 5000, 5995);
    dset_exp_2DHN = lib.getNMRExpnos(nmr_dataset_dirs, 4000, 4995);
end;

need_to_access_expnos = xaxis > 1;
mapper_expnos_provided = exist('mapper_expnos','var') && ~isempty(mapper_expnos);
if need_to_access_expnos
    % Repeat mapper expno lists if needed
    if mapper_expnos_provided
        % if there are fewer mapper_expnos than n_sets:
        if iscell(mapper_expnos) && numel(mapper_expnos)==1 && n_sets>1
            % Suggests that mapper_expnos is defined once for all datasets.
            % If already a cell - creating xaxis assuming different number
            % of expnos per dataset.
            mapper_expnos = repmat(mapper_expnos, n_sets, 1);
        elseif ~iscell(mapper_expnos)
            mapper_expnos = repmat({mapper_expnos}, n_sets, 1);
        end;      
    else % If mapper expnos not provided manually - use all
        mapper_expnos = dset_exp_2DHN;
    end;
end; % need_to_access_expnos

% Repeat peaklists if needed.
if numel(peaklist_paths)==1 && n_sets>1
    peaklist_paths = repmat(peaklist_paths, n_sets, 1);
end;

% Select only some datasets
if exist('select_dsets','var') && ~isempty(select_dsets)
    dset_names = dset_names(select_dsets);
    mapper_data_paths = mapper_data_paths(select_dsets);
    peaklist_paths = peaklist_paths(select_dsets);
    RNA_data_paths = RNA_data_paths(select_dsets);
    n_sets = numel(mapper_data_paths);
end

% fprintf('\n= Settings: %.2f s\n', toc); tic;

%% Assign trajectories to peaks
%=================================
if flag_run_assignments
    assign = cell(n_sets,1);
    for iSet=1:n_sets                        
        assign{iSet} = mapper.match_assign_with_traj_v04_intensity(dset_names{iSet},...
            peaklist_paths{iSet},...
            mapper_data_paths{iSet},...
            ''... % prot_seq_paths{iSet}
            );
    end; clear iSet;
end

% fprintf('\n= Run assign: %.2f s\n', toc); tic;

if flag_load_assignments
    assign = cell(n_sets,1);
    for iSet=1:n_sets
        tmp = load( fullfile('datasave', sprintf('assign_%s.mat', dset_names{iSet})) );
        assign{iSet} = tmp.assign;
    end
    clear tmp iSet;
end

% fprintf('\n= Load assign: %.2f s\n', toc); tic;

%% Import RNA conc and time. Interpolate RNA conc to mapper 2DHN times.
%======================================================================
if xaxis > 1 
    f_rna_data_exists = (exist('RNA_data_paths','var') && ~isempty(RNA_data_paths));
    assert(f_rna_data_exists, 'Abort: file with RNA & time data seems missing.');    
    
    tmp = cellfun(@(x) mapper.import_31P_n_time_from_csv(x), RNA_data_paths, 'un', 0);
    RNA_time = cellfun(@(x) x(:,1), tmp, 'un', 0);
    RNA = cellfun(@(x) x(:,2), tmp, 'un', 0);
    clear tmp;    
    HN_time = cellfun(@(d,e,t) lib.getNMRTime_02(d,e,t), nmr_dataset_dirs, mapper_expnos, dsets_time0, 'un', 0);    
    
    % This was only in the IN93a so far
    if logical(cell2mat(strfind(dset_names,'IN93a')))    
        non_zero_first_points_HN = logical(cell2mat(cellfun(@(x) x(1), HN_time, 'un', 0)));
        if any(non_zero_first_points_HN)
            warning('Correcting first timepoint in some HN data from non-zero to ZERO');
            HN_time{non_zero_first_points_HN}(1) = 0;
        end;
    end
    
    RNA_in_HN_time = cellfun(@(x,y,z) interp1(x,y,z, 'linear', 'extrap'), RNA_time, RNA, HN_time, 'un', 0);
    
end;
% fprintf('\n= RNA import: %.2f s\n', toc); tic;


%% Data to use for x-axis
%================
switch xaxis
    case 1
        optns.xaxis_label = 'spectrum id';        
        optns.xaxis = cellfun(@(x) [1:size(x.traj_HNcsp, 2)]', assign, 'un', 0);
    case 2
        optns.xaxis_label = 'time [min]';        
%         optns.xaxis = repmat(HN_time, n_sets, 1); % Needs to be a CELL ARRAY (if have many datasets!)
        optns.xaxis = HN_time;
    case 3
        optns.xaxis_label = 'RNA conc [mM]';
%         optns.xaxis = repmat(RNA_in_HN_time, n_sets, 1);
        optns.xaxis = RNA_in_HN_time;
end;


%% Visualization
%=================================
optns.plotLW = 2;
optns.plotSym = '-';
optns.markerSize = 4;
optns.same_min_xscale = 1;
optns.legend = dset_names;

% Select what is plotted.
base_name = [sprintf('%s_',dset_names{1:end-1}),dset_names{end}];
switch optns.yaxis
    case 1
        optns.yaxis_label = 'CSP';
        figure_name = sprintf('%s_csp', base_name);
        traj_select = cellfun(@(x) x.traj_HNcsp, assign, 'unif', 0);        
    case 2        
        optns.yaxis_label = 'peak amplitude';
        figure_name = sprintf('%s_intens', base_name);
        traj_select = cellfun(@(x) x.intens, assign, 'unif', 0);
    case 3
        optns.yaxis_label = 'CSP / amplitude';
        figure_name = sprintf('%s_csp_intens', base_name);

        traj_select{1} = assign{1}.traj_HNcsp;
    %     traj_select{2} = assign{1}.traj_HNcsp;
        traj_select{2} = assign{1}.intens;    

        n_sets = 2;    
        optns.legend = cellfun(@(x) strcat(optns.legend{1}, x), {': CSP', ': int'}, 'un', 0)';

        % Scale intensities such that they'd show normalized:
        traj_select{2} = bsxfun(@rdivide, traj_select{2}, max(traj_select{2},[],1));
        traj_select{2} = traj_select{2} .* max(max(traj_select{1},[],2)); % max value of all CSPs.    

        optns.xaxis = repmat(optns.xaxis(1), n_sets, 1);
        
    otherwise
        error('Abort: unrecognized Yaxis type');
end

% arrayfun(@(x) find([8 9]==x,1), assign{1}.ids, 'unif', 0)
% [tf, idx] = ismember([8,12], assign{1}.ids)    

% seq_end = numel(assign{1}.prot_seq);
% peak_names = assign{1}.prot_seq(assign{1}.ids(assign{1}.ids < seq_end));
% optns.titles = peak_names;
optns.titles = assign{1}.names;

%%% Replicate plotting options -- TODO: remove? (don't really need to
%%% define this here?
opt_to_duplicate = {'plotLW', 'plotSym'};
for iOpt=1:numel(opt_to_duplicate)
    if ~iscell(optns.(opt_to_duplicate{iOpt})) % encapsulate into cell if not a cell yet.
        optns.(opt_to_duplicate{iOpt}) = {optns.(opt_to_duplicate{iOpt})};
    end;            
    optns.(opt_to_duplicate{iOpt}) = repmat(optns.(opt_to_duplicate{iOpt}),n_sets,1);        
    
end;

name_already_has_ppm_digits = any(strfind(assign{1}.names{1},'.'));
if flag_show_peak_CS && ~name_already_has_ppm_digits
    optns.titles = arrayfun(@(name,H,N) {name{:}, sprintf('%.2f-%.1f', H, N)}, optns.titles, assign{1}.traj_H(:,1), assign{1}.traj_N(:,1), 'un', 0);
end;

%%% TMP TMP
optns.plotLW{2} = 1;

% fprintf('\n= Prepare traj: %.2f s\n', toc); tic;

%% Filter residues and/or time-points
%==========================================
% Select only subset of residues -- works only for TRAJECTORIES at the
% moment!
if exist('select_residues','var') && ~isempty(select_residues)
    %%% Use IDS if the peaklist had them defined
    %%% otherwise just peak numbers "as is"
    if isfield(assign{1},'ids')
        [~, idx] = ismember(select_residues, assign{1}.ids);        
    else
        idx = select_residues;
    end;
    
    optns.titles = optns.titles( idx );
    traj_select = cellfun(@(x) x(idx,:), traj_select, 'unif', 0);    
end

if exist('drop_timepoints','var') && ~isempty(drop_timepoints)
    % TODO - vectorize this
    for iSet = 1:n_sets
        traj_select{iSet}(:,drop_timepoints) = [];
        optns.xaxis{iSet}(drop_timepoints) = [];
    end;
end

% fprintf('\n= Filter traj: %.2f s\n', toc); tic;

%% Run Kd fits
%==========================================
if f_run_kd_fits
    
    
end; % f_run_kd_fits


%% Plot and save
%==========================================
% fig_handle = plot_traj3(traj_select, optns);
% fig_handle = mapper.plot_traj3b(traj_select, optns);
fig_handle = mapper.plot_traj4b(traj_select, optns);
mapper.save_figure(fig_handle, figure_name);

end