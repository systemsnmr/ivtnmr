function time0 = getTime0(nmr_data_path,dset_name)

%% Params for testing (when function is run "locally")
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);    
    nmr_data_path = '/Volumes/Data/yar/_eth2/data_NMR/spectra';
    dset_name = '160816_IN47a_p213-B3E_noRE_NUP1_303K_600';
end

%% Actual script
notes_path = fullfile(nmr_data_path,dset_name,'notes.txt');

filetext = fileread(notes_path);
pattern = '<time0>(.*?)</time0>';
time0 = regexp(filetext,pattern,'tokens');
assert(numel(time0) == 1, 'Abort: more than one time0 matches found in notes.txt.');    
time0 = char(time0{:}); % convert from cell to a string % 2017-07-19: added char()
clear filetext;

end