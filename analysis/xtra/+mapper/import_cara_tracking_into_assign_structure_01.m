function assign = import_cara_tracking_into_assign_structure_01(cara_repo_path, cara_project_name, cara_peaklist_name)


% Based on:
% match_assign_with_traj_v04_intensity.m -- creation of "assign" structure
% get_cara_peakshifts_04.m -- importing cara repo

DID NOT DO ANYTHING HERE YET!!!

% TODO:
% - combine the above "based on" files to make it generate the same
% "assign" structure to be used for plotting.

%% Params to run "locally"
%===========================
if nargin == 0
    dset_name = 'IN72b';
    project_root = '.';
    peaklist_filepath = fullfile(project_root,...
        'data_test',...
        '180319_UP1_v3_IN70_72_84_87_92.peaks');
    traj_filepath = fullfile(project_root,...
        'data_test',...
        'trajectories_filtered.csv');
    prot_seq_filepath = fullfile(project_root,...
        'data_test',...
        'UP1_1-196_1-letter.aa');
end

intens_filepath = strcat(traj_filepath(1:end-4),'_intensities.csv');

% Plot of CSP of first 10 peaks - checks whether there are indeed single minima, 
% (i.e. single assignment).
flag_visualize_peak_csp = 0;

%% Directory to keep intermediate results
%===========================
if ~exist('datasave', 'dir')
	mkdir('datasave');
end
datasave_folder = 'datasave';


%% Get peaks and trajectories
%===========================

% assign = import_NMR_peaklist_v02(peaklist_filepath, prot_seq_filepath);
% assign = import_NMR_peaklist_v02b(peaklist_filepath, prot_seq_filepath); % v02b - SORTS the peaklist by AA
% what mapper
assign = mapper.import_NMR_peaklist_v03b(peaklist_filepath, prot_seq_filepath); % v02b - SORTS the peaklist by AA
traj = mapper.import_trajectories(traj_filepath);
if exist(intens_filepath, 'file')
    intens = mapper.import_intensities(intens_filepath);
    flag_no_intensity_file_available = 0;
else
    fprintf(1,'No intensity file. Setting to empty.\n\n', mfilename);
    intens = '';
    flag_no_intensity_file_available = 1;
end


%% Match peaks to starting points of trajectories
%================================================
traj_time0 = cat(2, traj.H(:,1), traj.N(:,1)); % select HN at time0
traj_time0 = traj_time0(~sum(isnan(traj_time0),2),:); % drop all peaks with NaN

n_peaks = numel(assign.names);

if flag_visualize_peak_csp
    n_peaks_to_show = 10;
    colors = colormap(jet(n_peaks_to_show));
end

% For the moment - decided to use HN CSP formula - to assign based on
% lowest CSP value.
% (Maybe do someday (warning:perfectionism!!) - there probably is a way 
% to make this more efficient use vector dot-products).
trajH = traj_time0(:,1);
trajN = traj_time0(:,2);

assign.traj_id = nan(n_peaks,1);

for i=1:n_peaks
    % same formula as for combined HN CSP
    pk_csp = sqrt( ( (trajH-assign.H(i)).^2+((trajN-assign.N(i))./5).^2 )./2 ); 
    [~, traj_id] = min(pk_csp);
    assign.traj_id(i) = traj_id;
    
    if flag_visualize_peak_csp && i<=n_peaks_to_show
        semilogy(pk_csp,'Color',colors(i,:), 'LineWidth', 2);
        hold on;
        if i>1
            semilogy(pk_csp,'Color',colors(i,:), 'LineWidth', 2);
        end    
        xlabel('peak');
        ylabel('HN CSP');
    end   
    
end

% Save H/N/Hdelta/Ndelta/HNcsp info from trajectories to the matched peaks
assign.traj_H = traj.H(assign.traj_id, :);
assign.traj_N = traj.N(assign.traj_id, :);
assign.traj_Hdelta = traj.Hdelta(assign.traj_id, :);
assign.traj_Ndelta = traj.Ndelta(assign.traj_id, :);
assign.traj_HNcsp = traj.HNcsp(assign.traj_id, :);
% Calculate maxCSP:
assign.traj_maxHdelta = max(abs(assign.traj_Hdelta),[],2); % Need to make these ABSOLUTE!
assign.traj_maxNdelta = max(abs(assign.traj_Ndelta),[],2); % Need to make these ABSOLUTE!
assign.traj_maxHNcsp = max(assign.traj_HNcsp,[],2);

if ~flag_no_intensity_file_available
    assign.intens = intens(assign.traj_id, :); % also sorted by trajectory ID.
end

% Check unique assignments:
n_unique = length(unique(assign.traj_id));
if n_unique < n_peaks
    warning('\n=====\n %i / %i peaks have ambiguous assignments.\n=====\n',...
        n_peaks-n_unique, n_peaks);
end

% Save results
% if nargin == 0
    save(fullfile(datasave_folder, ...
        sprintf('assign_%s.mat',dset_name)),...
        'assign');
% end


end