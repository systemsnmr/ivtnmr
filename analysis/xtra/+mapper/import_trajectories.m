function traj = import_trajectories(traj_filepath)
% Previous imports: file:///Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra/170308_IN60a_pR02_co-A1R1_303K_600/chemshift_data_compilation/chemshift_titration_data_from_autom.m

%% Params to run "locally"
%===========================
if nargin == 0
%     project_root = '~/Dropbox/Science/Projects/2015_NMR_peak_tracking/team/mapper_tt';

%     project_root = '../..'; % "root" dir of the project - directory which contains code/ and data/ subdirs
%     filepath = 'data/171129_tracked_trajectories/IN60a/s3_tracking/trajectory.csv';    
%     traj_filepath = fullfile(project_root, filepath);    

    project_root = '/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/180323_SMN_AGU';
    filepath = 'IN72b/s3_tracking/trajectories_filtered.csv';  
    traj_filepath = fullfile(project_root, filepath);    
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

raw_list = importdata( traj_filepath );

% sel_pk_ids = [15 7 51 30 18 11 25 17 29 12 22 34];
% sel_peaks = raw_list(sel_pk_ids);
% names = {'G56','S22','G20','F59','G58','L21','H33','Q36','R75','V68','K15','T41'};

n_peaks = numel(raw_list);

all = regexp(raw_list,'(\d+\.\d+,\d+\.\d+)|(None)','match');
% could also do smth like:
% regexp(sel_peaks{1},';','split')

n_spectra = numel(all{1});


for i=1:n_peaks
    all{i} = regexprep(all{i},'None','nan,nan'); % replace Nones with NaN
    all{i} = regexp(all{i},',','split'); % split into individ numbers
    all{i} = cellfun(@str2num, [all{i}{:}], 'unif', 0); % flatten array
    all{i} = [all{i}{:}];
end

all = cell2mat(all);

peakNumbers = 1:n_peaks;

H = all(:,2:2:end); % get H shifts
N = all(:,1:2:end); % get N shifts
Htime0 = H(:,1); % shift at time0 - to calculate relative CSP later
Ntime0 = N(:,1); % shift at time0 - to calculate relative CSP later

Hdelta = bsxfun(@minus, H, H(:,1));
Ndelta = bsxfun(@minus, N, N(:,1));

HNcsp = sqrt( ( (Hdelta).^2+((Ndelta)./5).^2 )./2 ); % formula for combined HN CSP

% Matlab structure with shifts
traj.ids = 1:n_peaks;
traj.H = H;
traj.N = N;
traj.Hdelta = Hdelta;
traj.Ndelta = Ndelta;
traj.HNcsp = HNcsp;

disp('');

%% Cluster test:
%===========================
% kmeans(traj.HNcsp,3);


%% Export
%===========================
if nargin == 0 % save stuff only if running locally
    HNcsp_with_names = cat(2, peakNumbers', HNcsp);
    dlmwrite(fullfile(datasave_folder, strcat(mfilename,'_HN','.txt')), HNcsp_with_names, '	'); % save as delimeted text
    
    save(fullfile(datasave_folder, strcat(mfilename,'.mat')),'traj');
end

end