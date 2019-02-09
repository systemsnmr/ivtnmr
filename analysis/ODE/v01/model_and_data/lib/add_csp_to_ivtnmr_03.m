function add_csp_to_ivtnmr_03(csp_file_location)
% Alternative direction to run this - could specify nmr_data_path, dset_name, but then also
% need to specify csp_file_location

% Adds HN CSP data to an ivtnmr dataset.
% Currently reads intermediate data generated/saved by 170718_IN6x_IN7x_CSP_n_Kd_fits/collect_csp_data.m
% but can potentially combine with that script - so that it reads from cara
% and adds directly to ivtnmr.mat w/o intermediate storage.

% TODO
% add CSP to ivtnmr - csp_file_location
% + load CSP
% + for each CSP - load ivtnmr.mat
% + add CSP data
% + save ivtnmr.mat

% Versions:
% _03 - takes NMR_DATA_PATH as a global parameter

global DATA_IVTNMR NMR_DATA_PATH;

%% Params to run "locally"
%===========================
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    
    csp_file = 'csp_all.mat';
    fprintf(1,'\n\n= This file (%s) needs re-definition of datasave_folder (DATA_IVTNMR) param\n',mfilename);
    csp_file_location = fullfile(datasave_folder, csp_file);

else
    datasave_folder = DATA_IVTNMR;
    nmr_data_path = NMR_DATA_PATH;
end


%% Actual script
%===========================
load( csp_file_location );

n_sets = numel(fl);
name_split = cellfun(@(x) regexp(x.dset_name, '_', 'split'), fl, 'unif', 0);
snames = cellfun(@(x) x{2}, name_split, 'unif', 0); % extract short name

for i=1:n_sets
    fprintf(1,'== Processing %s ...\n', fl{i}.dset_name);
    fl{i}
% %     nmr_data_path = fl{i}.nmr_data_path; % _03 fix
    dset_name = fl{i}.dset_name;
    load(fullfile(nmr_data_path,dset_name,'ivtnmr.mat'));
    s
    s.hn
    
    %%% Quality checks (potentially could do):
    assert( strcmp(fl{i}.dset_name, s.name), 'Abort: HN-CSP and 31P dset_names not identical.' );
%     assert( sum(s.hn.time-fl{i}.time) == 0, 'Abort: HN times not matching.');
    assert( strcmp(s.time0, fl{i}.time0), 'Abort: time0 not matching.');

    %%% Update HN
    s.hn.time_all = s.hn.time;
    s.hn.expnos_all = s.hn.expnos;    
    s.hn.time = fl{i}.time;
    s.hn.expnos = (4000:fl{i}.last_expno)';
    s.hn.names = fl{i}.peak_tags';
    s.hn.H = fl{i}.H;
    s.hn.N = fl{i}.N;
    s.hn.csp = fl{i}.HN_csp;
    s.hn.n_peaks = fl{i}.n_peaks;
    s.hn.color = fl{i}.color;
    s.hn
    
%     name_split = regexp(fl{i}.dset_name, '_', 'split');
%     dset_snames = cellfun(@(x) x{2}, name_split, 'unif', 0); % extract short name

    %%% Save
    save( fullfile(datasave_folder, sprintf('ivtnmr_%s.mat', snames{i})), 's');
    save( fullfile(nmr_data_path, dset_name, 'ivtnmr.mat') , 's');
   
end

end
