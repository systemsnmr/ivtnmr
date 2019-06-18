function o = find_mapper_dir(nmr_dataset_dir)
% returns first dir with 'mapper' in the name (assumes there is only one)

%% Params to run "locally"
%===========================
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    
    current_dir = fileparts(mfilename('fullpath'));
    nmr_dataset_dir = fullfile(current_dir,'../..');
end


% disp(current_dir);
files = dir(fullfile(nmr_dataset_dir,'*mapper*'));

% https://ch.mathworks.com/matlabcentral/answers/166629-is-there-any-way-to-list-all-folders-only-in-the-level-directly-below-a-selected-directory
names = {files.name};
dir_flags = [files.isdir];
mapper_dirs = names(dir_flags);
o = mapper_dirs{1}; % take the first one, assuming its the only

end