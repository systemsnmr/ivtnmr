function intens = import_intensities(intens_filepath)
% based on import_trajectories

%% Params to run "locally"
%===========================
if nargin == 0
%     project_root = '~/Dropbox/Science/Projects/2015_NMR_peak_tracking/team/mapper_tt';
    project_root = '/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/180323_SMN_AGU';
    filepath = 'IN72b/s3_tracking/trajectories_filtered_intensities.csv';  
%     filepath = 'IN72b/s3_tracking/trajectories_filtered.csv';  
    intens_filepath = fullfile(project_root, filepath);    
else
    
end

%% Directory to keep intermediate results
%===========================
if ~exist('datasave', 'dir')
	mkdir('datasave');
end
datasave_folder = 'datasave';

%% Parsing
%===========================
% Had to modify the procedure here - just "importdata" did not work!
fid = fopen(intens_filepath);
tline = fgetl(fid);
raw_list = cell(0,1);
while ischar(tline)
    raw_list{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

n_peaks = numel(raw_list);

% all = regexp(raw_list,'(\d+\.\d+,\d+\.\d+)|(None)','match');
all = regexp(raw_list,'(\d+\.\d+)|(None)','match');

% could also do smth like:
% regexp(sel_peaks{1},';','split')

n_spectra = numel(all{1});

for i=1:n_peaks
    all{i} = regexprep(all{i},'None','nan'); % replace Nones with NaN
    all{i} = regexp(all{i},',','split'); % split into individ numbers
    all{i} = cellfun(@str2num, [all{i}{:}], 'unif', 0); % flatten array
    all{i} = [all{i}{:}];
end

intens = cell2mat(all);

%% Export
%===========================
if nargin == 0 % save stuff only if running locally
%     HNcsp_with_names = cat(2, peakNumbers', HNcsp);
%     dlmwrite(fullfile(datasave_folder, strcat(mfilename,'_intens','.txt')), HNcsp_with_names, '	'); % save as delimeted text
%     
%     save(fullfile(datasave_folder, strcat(mfilename,'.mat')),'traj');
end

end