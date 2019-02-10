function expnos = getExpnos(nmr_data_path,dset_name,data_type)

%% Params for testing (when function is run "locally")
if nargin == 0
    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);    
    nmr_data_path = '/Volumes/Data/yar/_eth2/data_NMR/spectra';
    dset_name = '170308_IN60a_pR02_co-A1R1_303K_600';
    data_type = 'hsqc';
    
end


%% Actual script
%===================
switch data_type
    case {'1H', 'zgwg', '1D1H'}
        expno_start = 2000;
    case {'TOCSY','tocsy'}
        expno_start = 3000;
    case {'HSQC', 'hsqc', '2DHN', 'HN'}
        expno_start = 4000;
    case {'31P','P31','phos'}
        expno_start = 5000;
    case {'im', 'imino', 'iminos'}
        expno_start = 6000;
    otherwise
        % User can provide the name of the file instead
        expno_start = data_type;
end

disp('');

first_digit = num2str(expno_start);
first_digit = first_digit(1);

%%% Separate
% dirlist = dir(fullfile(nmr_data_path, dset_name, strcat(first_digit,'*')));
% e_str = arrayfun(@(x) str2num(x.name), dirlist, 'unif', 0);
% expnos = cell2mat(e_str);
%%% Combined
expnos = cell2mat(arrayfun(@(x) str2num(x.name), dir(fullfile(nmr_data_path, dset_name, strcat(first_digit,'*'))), 'unif', 0));

% Filter out only experiments within x000 till x950
expnos = expnos(expnos >= expno_start & expnos <= expno_start+950);

%%% Stuff from getTime0
% notes_path = fullfile(nmr_data_path,dset_name,'notes.txt');
% filetext = fileread(notes_path);
% pattern = '<time0>(.*?)</time0>';
% time0 = regexp(filetext,pattern,'tokens');
% assert(numel(time0) == 1, 'Abort: more than one time0 matches found in notes.txt.');    
% time0 = char(time0{:}); % convert from cell to a string % 2017-07-19: added char()
% clear filetext;

end