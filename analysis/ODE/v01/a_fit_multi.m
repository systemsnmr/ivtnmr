function a_fit_multi
% Controller script - to fit multiple IVTNMR datasets.
% When changing, check:
% - model_sub_folder
% - fit_func_name

close all;

%% Settings
%=======================
model_sub_folder = 'model_and_data';
fit_func_name = 'a_fit';

data_sets = {...
'170310_IN70a_pR02_co-UP1_303K_600'
% '180702_IN70b_pR02_co-NUP1_303K_600'
% '170313_IN71a_p114_co-UP1_303K_600'
% '180914_IN71b_SMN1_co-NUP1_303K_600'
% '170426_IN72a_SMN2_co-UP1_303K_600'
'180314_IN72b_SMN214_co-NUP1_303K_600'
'170322_IN75a_EV2_co-UP1_303K_600'  
% '170515_IN75b_EV2_co-UP1_303K_600'
% '190111_IN70c_R02_co-NUP1_303K_600'
% '190110_IN71c_SMN1_co-NUP1_303K_600'
% '190111_IN72c_SMN2_co-NUP1_303K_600'
% '190117_IN72d_SMN2_co-NUP1_303K_600'
% '190112_IN75h_EV2_co-NUP1_303K_600'
% '190113_IN75z_EV2_co-NUP1_303K_600'
};

% hn_peaks_for_fit = {'H33' 'R75'};
hn_peaks_for_fit = {'all'};

data_suffix = 'intern_HN';

%%% Use these options if (1) protein is not saturated with RNA during reaction time
%%% and (2) you don't have a good estimate for maximum chemical
%%% shift perturbations of specific residues at saturation. Then can make an assumption
% optns.sets_for_maxcsp = {'IN70a' 'IN70b' 'IN70c' 'IN71a' 'IN71b' 'IN71c'...
%     'IN72a' 'IN72b' 'IN72c' 'IN72d' 'IN75a' 'IN75b' 'IN75h' 'IN75z'};
optns.peaks_for_maxcsp = []; % if empty - will use all peaks
%%% if you have a linear CSP in the course of reaction, assuming that maximum
%%% shift is < than 1.9 is unrealistic. maxcsp_scaling = 2 seems a reasonable
%%% lower limit.
optns.maxcsp_scaling = 2; 

% This overrides the above settings and uses fixed CSP values -
% need to be measured experimentally/independently.
optns.maxcsp_fixed = 1;
optns.maxcsp_list = containers.Map({...
        'H33','R75'},...
        [0.083, 0.046]);

optns.rna_fraction = 0.3; % is used to normalize RNA/Aborts fractions
optns.T7DNAcomplex_concentration = 33*1e-6; % mM
optns.Mg_total = 24; % mM
optns.NTP_initial_conc = 20; % mM

% Bootstrap stuff (not used anymore!)
optns.useBootstrap = 0; % if set to 0 - runs regular fit with full dataset
optns.nboot = 200; % number of fits


%% Script
%=======================
n_sets = numel(data_sets);
n_peaks = numel(hn_peaks_for_fit);

fit_data = cell(n_sets,n_peaks); % preallocate for sets and subsets

tic;
for i=1:n_sets
    dsetname = regexp(data_sets{i}, '_', 'split');
    dsetname = dsetname{2};
    
    for p=1:n_peaks
        fprintf(1,'\n\n== Fitting %s with peak %s...\n\n', dsetname, hn_peaks_for_fit{p});
        
        optns.HN_peak_for_fit = hn_peaks_for_fit{p};
        data_name = sprintf('%s_%s',dsetname,data_suffix);
        
        fit_data{i,p} = feval(fit_func_name, data_name, optns);
        
    end
end; clear i;

fprintf(1,'\n\n== Time for %i sets = %.2f min\n\n', n_sets, toc/60);

if ~optns.useBootstrap % do not need to save all together - each bootstrap is saved by itself
    results_dir = fullfile(model_sub_folder, 'fit_results');
    if ~exist( results_dir ,'dir'), mkdir(results_dir), end;

    save_mat_name = fullfile(results_dir, sprintf('%s_multi.mat', fit_func_name));
    fprintf(1,'Data saved into %s\n', save_mat_name);        
    save(strcat(save_mat_name),'fit_data');
end

end
