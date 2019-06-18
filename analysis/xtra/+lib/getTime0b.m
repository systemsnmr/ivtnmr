function time0 = getTime0b(nmr_data_path,dset_name)

% getTime0b:
% - can take either full path or two parts (path + dset_name).

%% Params for testing (when function is run "locally")
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);    
    nmr_data_path = '/Volumes/Data/yar/_eth2/data_NMR/spectra';
    dset_name = '190311_IN111a_pA20_free_BK';
    
    % Test for single-input
    if 0
        nmr_data_path = fullfile(nmr_data_path,dset_name);
        dset_name = [];    
    end;
end

%% Actual script
if ~exist('dset_name','var') || ~isempty(dset_name)
    notes_path = fullfile(nmr_data_path,'notes.txt');
else
    notes_path = fullfile(nmr_data_path,dset_name,'notes.txt');
end;
    
filetext = fileread(notes_path);
pattern = '<time0>(.*?)</time0>';
time0 = regexp(filetext,pattern,'tokens');
assert(numel(time0) == 1, 'Abort: more than one time0 matches found in notes.txt.');    
time0 = char(time0{:}); % convert from cell to a string % 2017-07-19: added char()
clear filetext;

end