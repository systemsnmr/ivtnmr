function maxcsp_map = get_maxcsp_for_peaks(nmr_data_paths, hn_name_list)
% Returns all HN peak names in a group of "dsel" datasets
% Modified from ivtnmr_get_maxcsp_for_peaks.m - which works for "ivtnmr"
% data structures.

%% Params to run "locally"
%===========================
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    base_path = '/Volumes/Data/yar/Dropbox/_eth2/data_MathModeling/170614_IN60_IN70/v10_v9_Prot/data';
    dset_names = {'IN70a' 'IN71a' 'IN72a' 'IN75b'};
    nmr_data_paths = cellfun(@(x) fullfile(base_path, sprintf('data_nmr_selected_%s_intern_HN.mat', x)), dset_names, 'un', 0);    
    hn_name_list = [];
    
    plot_prot_vs_RNA = 1;
end

%% Actual script
%===========================
n_sets = numel(nmr_data_paths);
d = cell(n_sets,1);

% name_split = cellfun(@(x) regexp(x, '_', 'split'), dset_names, 'unif', 0);
% dset_snames = cellfun(@(x) x{2}, name_split, 'unif', 0); % extract short name

% hn_names = cellfun(@(x) x.hn.names, d, 'un', 0);
% hn_names = [hn_names{:}];

% Load sets
for i=1:n_sets
    load(nmr_data_paths{i});
    d{i} = dsel;
end

if isempty(hn_name_list) % checking that all sets have the same number and names of peaks
    pk_names = d{1}.hn.names;
    n_peaks = numel(pk_names);
    if sum(cell2mat(cellfun(@(x) numel(x.hn.names)~=n_peaks, d, 'un', 0)))
        error('Number of peaks must be same in all sets.');
    end
    if sum(~cellfun(@(x) isequal(pk_names,x.hn.names), d))
        error('Names of peaks must be same in all sets.');
    end    
    hn_name_list = pk_names;
end

n_peaks = numel(hn_name_list);

%%%%%%%%%%% Just checking if its really what expected
if exist('plot_prot_vs_RNA', 'var') && plot_prot_vs_RNA
    % (NTP_init - NTP)*rna_fraction / rna_length
    rna_fraction = 0.3;
    rna_conc = cellfun(@(x) (20-x.y(2,:)).*rna_fraction./x.hn.rna_length, d, 'un', 0);
    rows = 2;
    columns = n_peaks;
    sbSide = 200;
    fSize = 11;
    figure(1); set(figure(1), 'Color', ones(1,3), 'Position', [0 300 columns*sbSide rows*sbSide]);
    sb = @(x) subplot(rows,columns,x,'FontSize',fSize);
    colors = jet(n_sets);

    csp = cell(1,n_peaks);
    for p=1:n_peaks
        csp{p} = cellfun(@(x) x.hn.mean(:,p), d, 'un', 0);
        sb(p)
        for i=1:n_sets
            hold on;
            plot(rna_conc{i}, csp{p}{i}, 'o-', 'Color', colors(i,:));
            axis tight;
            ylim_curr = ylim;
            xlim([0 0.3]);
            ylim([0 ylim_curr(2)*2]);
        end
        if p==1
            legend(dset_names);
        end
        title(hn_name_list(p))
    end; clear p;
end
%%%%%%%%%%% Just checking if its really what expected


% Get max csp for each peak
maxcsp_list = nan(n_peaks,1);

for p=1:n_peaks
    pkname = hn_name_list{p};
    if nargin==0, fprintf(1, 'peak = %s\n', pkname), end;
    pkidx_list = cellfun(@(x) strcmp(pkname, x.hn.names), d, 'un', 0); % indexes of the peak across all dsets
    pkcsp_list = cellfun(@(x,y) x.hn.mean(:,y), d, pkidx_list, 'un', 0); % csp of the peak in all dsets
    non_empty_sets = cell2mat( cellfun(@(x) ~isempty(x), pkcsp_list, 'un', 0) );
    maxcsp = max( cell2mat( cellfun(@(x) max(x), pkcsp_list(non_empty_sets), 'un', 0) ) );
    maxcsp_list(p) = maxcsp;
end

maxcsp_map = containers.Map(hn_name_list, maxcsp_list);
% maxcsp_map.keys';
% maxcsp_map.values';

disp('');
    
end

