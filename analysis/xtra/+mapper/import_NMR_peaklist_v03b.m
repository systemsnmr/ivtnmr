function assign = import_NMR_peaklist_v03b(peaklist_filepath, prot_seq_filepath)
% rawList = importdata('shell_cat_all_txts.txt');
% started from collect_csp_data_02

% Versions:
% import_NMR_peaklist.m - initial
% import_NMR_peaklist_v02.m:
%   - added importing of protein sequence
% import_NMR_peaklist_v02b
%   - removed ID field & added sorting of initial peaklist by residue #
% import_NMR_peaklist_v03
%   - added a feature to check if the imported peak contains peak names (#... lines) or
%   just peak positions. If its with names - it imported as ONE cell.
%   - if there are no names - make a name of peak CS!
% import_NMR_peaklist_v03b - 20190411 - YarN
%   - added check for extra #-starting line at the beginning of the
%   peaklist -- see % f_peaklist_with_double_lines 


% TODO
% - 
% Maybe:
% - implement such that it determines the type of peaklist
% automatically (or expects as input).

%% Params to run "locally"
%===========================
if nargin == 0
    
	% Switch to the working directory
    [fld,~,~] = fileparts( mfilename('fullpath') );
    cd(fld);
%     project_root = '~/Dropbox/Science/Projects/2015_NMR_peak_tracking/team/mapper_tt';
     project_root = '../..'; % "root" dir of the project - directory which contains code/ and data/ subdirs
   
    peaklist_local_path = 'data/171129_A1R1_seq_n_assign/assignments/170718_A1R1_PeakVsResidue_IN60a.peaks';
    peaklist_filepath = fullfile(project_root, peaklist_local_path);
    prot_seq_local_path = 'data/171129_A1R1_seq_n_assign/sequence/A1R1_1-100_1-letter.aa';    
    prot_seq_filepath = fullfile(project_root,prot_seq_local_path);
    

    peaklist_filepath = fullfile('/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/180323_SMN_AGU/analysis/data_test/180319_UP1_v3_IN70_72_84_87_92.peaks');
    prot_seq_filepath = fullfile('/Volumes/Data/yar/Dropbox/Science/Projects/2015_NMR_peak_tracking/site_mapper/180323_SMN_AGU/analysis/data_test/UP1_1-196_1-letter.aa');
    
    datasave_name = strcat(mfilename,'.mat');
else
    
end

%% Directory to keep intermediate results
%===========================
if ~exist('datasave', 'dir')
	mkdir('datasave');
end
datasave_folder = 'datasave';

fprintf(1,'\n== CHECK import_NMR_peaklist_v03c - more robust in distinguish peaklist structure\n\n');

%% Import data
%===========================

fileID = fopen( peaklist_filepath );
c = textscan(fileID, '%s', 'delimiter', '\n'); % import all lines
headerlines = 3;
c = c{1}(headerlines+1:end);

if numel(c)==1 % Imported as one cell if XEASY peaklist has "#..." lines with peak names
    f_peaklist_with_double_lines = 1;
elseif strcmp(c{1,1}(1),'#') % if there is extra header line
    f_peaklist_with_double_lines = 1;
    c = c(2:end); 
else
    f_peaklist_with_double_lines = 0;
end;
    
if f_peaklist_with_double_lines
    n_peaks = size(c,1)/2;

    assign.ids = nan(n_peaks,1);
    assign.names = cell(n_peaks,1);
    assign.H = nan(n_peaks,1);
    assign.N = nan(n_peaks,1);

    for i=1:n_peaks
        pk_line = i*2-1;
    %%% Example of data in the peaks file:
    %%% first three elements in the first line, and the peak name in the second
    %%% line too.
    %      8   8.444 121.818 1 U   0.000E+00   0.000E+00 e 0    13    14
    % # K8
        % first line - peak info, second line - peak name
        pk_info = textscan(c{pk_line}, '%d   %f %f %*d %*s   %*f%*s%*d   %*f%*s%*d %*s %*d    %*d    %*d');
        pk_name = textscan(c{pk_line+1}, '# %s');

        [assign.ids(i), assign.H(i), assign.N(i), assign.names(i)] = ...
            deal(pk_info{1}, pk_info{2}, pk_info{3}, pk_name{1});    
    end; clear i;
else
    n_peaks = numel(c);

    assign.ids = nan(n_peaks,1);
    assign.names = cell(n_peaks,1);
    assign.H = nan(n_peaks,1);
    assign.N = nan(n_peaks,1);

    for i=1:n_peaks
        pk_line = i;
    %%% Example of data in the peaks file:
    %%% first three elements in the first line, and the peak name in the second
    %%% line too.
    %      1       8.1274     119.6781 1 -    0.000E00    0.000E00 -  0     0     0     0
    % # K8
        % first line - peak info, second line - peak name
        pk_info = textscan(c{pk_line}, '%d   %f %f %*d %*s   %*f%*s%*d   %*f%*s%*d %*s %*d    %*d    %*d');
        pk_name = {{sprintf('%.2f-%.1f', pk_info{2}, pk_info{3})}};

        [assign.ids(i), assign.H(i), assign.N(i), assign.names(i)] = ...
            deal(pk_info{1}, pk_info{2}, pk_info{3}, pk_name{1});    
    end; clear i;
    
end


%% Sort the peak by residue #
%===========================
assign = rmfield(assign,'ids');
% find all AA-style peak names
aa_ids = ~cellfun(@(x) isempty(x), regexp( assign.names ,'^[A-Z]\d+$'));

% TODO - can perhaps just get away by working with indices, w/o splitting
% into two arrays
% split peaklist in two parts:
aa.names = assign.names(aa_ids);
aa.H = assign.H(aa_ids);
aa.N = assign.N(aa_ids);
nonaa.names = assign.names(~aa_ids);
nonaa.H = assign.H(~aa_ids);
nonaa.N = assign.N(~aa_ids);

% sort the AA-style peaks
ids_to_sort = cell2mat( cellfun(@(x) str2num(x(2:end)), aa.names, 'un', 0) );
[~, sorted_indexes] = sort(ids_to_sort);
aa.names = aa.names(sorted_indexes);
aa.H = aa.H(sorted_indexes);
aa.N = aa.N(sorted_indexes);

% reassemble into one array
assign.names = [aa.names; nonaa.names];
assign.H = [aa.H; nonaa.H];
assign.N = [aa.N; nonaa.N];

%% Import prot seq
%===========================
if exist('prot_seq_filepath','var') && ~isempty(prot_seq_filepath)
    seq = mapper.import_prot_sequence(prot_seq_filepath);
    assign.prot_seq = seq.aa;
    assign.prot_seq_ids = seq.num;
end

%% Export data
%===========================
% Save Matlab structure with shifts - only if running locally
if nargin == 0
    save(fullfile(datasave_folder, strcat(mfilename,'.mat')),'assign');
end

end % main func
