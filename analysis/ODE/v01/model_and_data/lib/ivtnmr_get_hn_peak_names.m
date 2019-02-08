function hn_names = ivtnmr_get_hn_peak_names(nmr_data_path, dset_names)
% Returns all HN peak names in a group of ivtnmr datasets

%% Params to run "locally"
%===========================
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    nmr_data_path = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra';
    dset_names = {...
        '170308_IN60a_pR02_co-A1R1_303K_600'
        '170310_IN70a_pR02_co-UP1_303K_600'
        '170503_IN61a_SMN1_co-A1R1_303K_600'
        '170313_IN71a_p114_co-UP1_303K_600'
        '170309_IN62a_p214_co-A1R1_303K_600'
        '170426_IN72a_SMN2_co-UP1_303K_600'
        '170607_IN65b_EV2_co-A1R1_303K_600'
        '170515_IN75b_EV2_co-UP1_303K_600'
        };
end

%% Actual script
%===========================
n_sets = numel(dset_names);
d = cell(n_sets,1);

% name_split = cellfun(@(x) regexp(x, '_', 'split'), dset_names, 'unif', 0);
% dset_snames = cellfun(@(x) x{2}, name_split, 'unif', 0); % extract short name

for i=1:n_sets
    load(fullfile(nmr_data_path,dset_names{i},'ivtnmr.mat'));
    d{i} = s;
end

hn_names = cellfun(@(x) x.hn.names, d, 'un', 0);
hn_names = [hn_names{:}];

disp('');
    
end

