function expnos = getNMRExpnos(nmr_dataset_path, expno_range_start, expno_range_end)
% Extract expnos from NMR dataset. Can process multiple datasets - if
% provided in cell array.

%% Params for testing (when function is run "locally")
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);    
    nmr_dir = '/Volumes/Data/yar/_eth2/data_NMR/spectra';
    dset_name = '190311_IN111a_pA20_free_BK';
    nmr_dataset_path = fullfile(nmr_dir,dset_name);
    expno_range_start = 5000;
    expno_range_end = 5950;
end

%% Actual script
%===================
first_digit = num2str(expno_range_start);
first_digit = first_digit(1);

if ~iscell(nmr_dataset_path)
    nmr_dataset_path = {nmr_dataset_path};
end;

n_sets = numel(nmr_dataset_path);

expnos = cell(n_sets,1);

for iSet=1:n_sets
    tmp = cell2mat(arrayfun(@(x) str2num(x.name), dir(fullfile(nmr_dataset_path{iSet}, strcat(first_digit,'*'))), 'un', 0));
    expnos{iSet} = tmp(tmp >= expno_range_start & tmp < expno_range_end);
    clear tmp;
end;

end