function maxcsp_map = ivtnmr_get_maxcsp_for_peaks(nmr_data_path, dset_names, hn_name_list)
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
    
    hn_name_list = [];
end

%% Actual script
%===========================
if isempty(hn_name_list) % get all peak names
    hn_name_list = unique( ivtnmr_get_hn_peak_names(nmr_data_path, dset_names) );    
end

n_sets = numel(dset_names);
n_peaks = numel(hn_name_list);
d = cell(n_sets,1);

% name_split = cellfun(@(x) regexp(x, '_', 'split'), dset_names, 'unif', 0);
% dset_snames = cellfun(@(x) x{2}, name_split, 'unif', 0); % extract short name

% hn_names = cellfun(@(x) x.hn.names, d, 'un', 0);
% hn_names = [hn_names{:}];

% Load sets
for i=1:n_sets
    load(fullfile(nmr_data_path,dset_names{i},'ivtnmr.mat'));
    d{i} = s;
end

% Get max csp for each peak
maxcsp_list = nan(n_peaks,1);

for p=1:n_peaks
    pkname = hn_name_list{p};
    fprintf(1, 'peak = %s\n', pkname);
    pkidx_list = cellfun(@(x) strcmp(pkname, x.hn.names), d, 'un', 0); % indexes of the peak across all dsets
    pkcsp_list = cellfun(@(x,y) x.hn.csp(y,:), d, pkidx_list, 'un', 0); % csp of the peak in all dsets
    non_empty_sets = cell2mat( cellfun(@(x) ~isempty(x), pkcsp_list, 'un', 0) );
    maxcsp = max( cell2mat( cellfun(@(x) max(x), pkcsp_list(non_empty_sets), 'un', 0) ) );
    maxcsp_list(p) = maxcsp;
end

maxcsp_map = containers.Map(hn_name_list, maxcsp_list);
% maxcsp_map.keys';
% maxcsp_map.values';

disp('');
    
end

