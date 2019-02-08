function cds_multi(optns)
% Controller script - to generate multiple IVTNMR data matrices for model
% fitting.

global DATA_IVTNMR;

tic;
if nargin == 0
    optns.flag_exlude_hsqc_for_initial_dsel_generation = 1;

    select_sets = {'IN72b' 'IN75a'};

    data_sets = {...
    '170313_IN71a_p114_co-UP1_303K_600'
    '170426_IN72a_SMN2_co-UP1_303K_600'
    '170515_IN75b_EV2_co-UP1_303K_600'
    '170310_IN70a_pR02_co-UP1_303K_600'
    '180702_IN70b_pR02_co-NUP1_303K_600'
    '180314_IN72b_SMN214_co-NUP1_303K_600'
    '170322_IN75a_EV2_co-UP1_303K_600'
    '180914_IN71b_SMN1_co-NUP1_303K_600'
    '190110_IN71c_SMN1_co-NUP1_303K_600'
    '190111_IN70c_R02_co-NUP1_303K_600'
    '190111_IN72c_SMN2_co-NUP1_303K_600'
    '190112_IN75h_EV2_co-NUP1_303K_600'
    '190113_IN75z_EV2_co-NUP1_303K_600'
    '190117_IN72d_SMN2_co-NUP1_303K_600'
    };

    if exist('select_sets','var') && ~isempty(select_sets)
        [~, sel_idx] = ismember(select_sets, cellfun(@(x) x(8:12), data_sets, 'un', 0));
        data_sets = data_sets(sel_idx);
    end
else
    data_sets = optns.dset_names;
end;

dset_id_names = cellfun(@(x) x(8:12), data_sets, 'un', 0);

hn_peak_names = {'H33', 'R75'};
hn_optns.p_conc = 0.150; % mM - protein concentration used

n_sets = numel(data_sets);
export_files = 1;

% ivtnmr_path_root = '../../ivtnmr_datasets';
ivtnmr_path_root = DATA_IVTNMR;
ivtnmr_paths = cellfun(@(x) fullfile(ivtnmr_path_root, sprintf('ivtnmr_%s.mat', x)), dset_id_names, 'un', 0)';

for i=1:n_sets
    fprintf(1,'\n\n== Processing %s...\n\n', data_sets{i});
    hn_optns.ivtnmr_path = ivtnmr_paths{i};
    hn_optns.hn_peak_names = hn_peak_names;
    
    if optns.flag_exlude_hsqc_for_initial_dsel_generation
        hn_optns = [];
    end
    
%     create_data_str4(data_sets(i), 0, 'intern', export_files, hn_optns);
    create_data_str5(data_sets(i), 0, 'intern', export_files, hn_optns, optns);
end

toc

end
