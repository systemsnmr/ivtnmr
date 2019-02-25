clear; close all;

%% Settings to check
%==========================
global FIG; FIG=0;

global NMR_DATA_PATH DATA_IVTNMR;
analysis_root = fullfile(pwd, '..', '..', '..');
NMR_DATA_PATH = fullfile(analysis_root, '..', 'spectra');
DATA_IVTNMR = fullfile(pwd, 'data_ivtnmr_full');

optns.cara_path = fullfile(analysis_root, 'csp');
optns.cara_repo_name = '190121_UP1_repl.cara';

init_NTP = 20; % mM

rna_lengths = containers.Map({...
    'IN60','IN61','IN62','IN63','IN65','IN70','IN71','IN72','IN73','IN75','IN78'},...
    [5.8, 28, 28, 33, 49, 5.8, 28, 28, 33, 49, 5.8]);

dset_names = {...
'170310_IN70a_pR02_co-UP1_303K_600'
'180702_IN70b_pR02_co-NUP1_303K_600'
'190111_IN70c_R02_co-NUP1_303K_600'

'170313_IN71a_p114_co-UP1_303K_600'
'180914_IN71b_SMN1_co-NUP1_303K_600'
'190110_IN71c_SMN1_co-NUP1_303K_600'

'170426_IN72a_SMN2_co-UP1_303K_600'
'180314_IN72b_SMN214_co-NUP1_303K_600'
'190111_IN72c_SMN2_co-NUP1_303K_600'
'190117_IN72d_SMN2_co-NUP1_303K_600'

'170322_IN75a_EV2_co-UP1_303K_600'  
'170515_IN75b_EV2_co-UP1_303K_600'
'190112_IN75h_EV2_co-NUP1_303K_600'
'190113_IN75z_EV2_co-NUP1_303K_600'
};

analyze_subset = [8 11];
% analyze_subset = [1 8 11];
% analyze_subset = [1];


%% Which sections to run:
%==========================
run_all = 1; % makes as if all below are TRUE
% ---
run_cds_initial = 0;
    cds_optns.plot_hsqc = 0;
    cds_optns.plot_p31 = 0;
    cds_optns.plot_final = 0;    
    cds_optns.data_for_fit = fullfile(pwd,'data_for_fit');

run_collect_csp = 0;
    export_data = 1;
    plot_csp_results = 0;
    
run_gen_ivtnmr_from_dsel = 1;
    plot_gen_ivtnmr_results = 1;
    save_to_global_ivtnmr = 1; % will save into NMR-spectra directory too.

run_add_csp = 1;

run_interpolate_p31 = 1;
    interp_optns.save_to_global_ivtnmr = 1;
    interp_optns.plot_results = 0;

run_cds_final = 0;

%% Normally 'stable' settings
%==========================
dset_id_names = cellfun(@(x) regexp(x,'_','split'), dset_names, 'un', 0);
dset_id_names = cellfun(@(x) x{2}, dset_id_names, 'un', 0);

% Assuming projects and peaklists are named by experiment ID (IN##) (can define otherwise
% if need to!)
projects = dset_id_names;
peaklists = dset_id_names;

if exist('analyze_subset', 'var') && ~isempty(analyze_subset)
    dset_names = dset_names(analyze_subset);
    projects = projects(analyze_subset);
    peaklists = peaklists(analyze_subset);
end;


%% Main code
%==========================
addpath( genpath(fullfile(pwd,'lib')) );

% Create dirs if they don't exist
if ~exist( cds_optns.data_for_fit ,'dir'), mkdir( cds_optns.data_for_fit ), end;
if ~exist( DATA_IVTNMR ,'dir'), mkdir( DATA_IVTNMR ), end;


%% 
%==================================
if run_cds_initial || run_all
    fprintf(1,'\n\n===\n== CDS initial...\n===\n\n');
    % Initial creation of 'data_for_fit' - needed as a base to create full
    % 'ivtnmr' dataset files
    cds_optns.flag_exlude_hsqc_for_initial_dsel_generation = 1;
    cds_optns.dset_names = dset_names;
    cds_optns.nmr_data_path = NMR_DATA_PATH;
    cds_multi(cds_optns);
end;


%% Collect CSP data from cara repo
%==================================
if run_collect_csp || run_all
    fprintf(1,'\n\n===\n== Collecting CSP from Cara...\n===\n\n');
    collect_csp_data_02b(optns.cara_path, optns.cara_repo_name, projects, peaklists, export_data, plot_csp_results);    
end


%% Generate ivtnmr
%===================
if run_gen_ivtnmr_from_dsel || run_all
    fprintf(1,'\n\n===\n== Generating full IVTNMR object/datastructure...\n===\n\n');

    dsel_path = cds_optns.data_for_fit;
    
    dsel_suffixes = {'intern'};
    
    gen_ivtnmr_from_dsel_02(dsel_path, dsel_suffixes, init_NTP, ...
        NMR_DATA_PATH, dset_names, plot_gen_ivtnmr_results, save_to_global_ivtnmr, rna_lengths)
end


%% Add CSP to ivtnmr
%===================
if run_add_csp || run_all
    fprintf(1,'\n\n===\n== Adding CSP to IVTNMR datastructure...\n===\n\n');
    datasave_folder = DATA_IVTNMR;
    csp_file = 'csp_all.mat';
    csp_file_location = fullfile(datasave_folder, csp_file);
    add_csp_to_ivtnmr_03(csp_file_location);
end

%% Interpolate
%===================
if run_interpolate_p31 || run_all
    fprintf(1,'\n\n===\n== Interpolating 31P data to HN (if want manual fits)...\n===\n\n');
    [~, ~] = interpolate_p31_obs_hn_time_03(NMR_DATA_PATH, dset_names, interp_optns);
end


%% Make final data structure for model fitting - including 31P and 2DHN
%==================================
if run_cds_final || run_all
    fprintf(1,'\n\n===\n== CDS final: exporting data for model fits...\n===\n\n');
    % Second round of create_data_structure - now with adding CSP data from
    % the full 'ivtnmr datasets'
    cds_optns.flag_exlude_hsqc_for_initial_dsel_generation = 0;
    cds_multi(cds_optns);
end;



