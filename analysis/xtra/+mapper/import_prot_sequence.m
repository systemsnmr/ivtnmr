function seq = import_prot_sequence(sequence_filepath)
% rawList = importdata('shell_cat_all_txts.txt');
% started from import_NMR_peaklist

%% Params to run "locally"
%===========================
if nargin == 0    
	% Switch to the working directory
    [fld,~,~] = fileparts( mfilename('fullpath') );
    cd(fld);

%     project_root = '~/Dropbox/Science/Teaching/Students/2017_Slaven/NMR_mapper_shared';
%     data_dir = 'data';
%     filepath = 'Prot_sequences_assignments_PDB/Protein2_UP1/sequence/UP1_1-196_1-letter.aa';
%     sequence_filepath = fullfile(project_root,data_dir,filepath);        

    project_root = '~/Dropbox/Science/Projects/2015_NMR_peak_tracking/team/mapper_tt';
%     project_root = '../..'; % "root" dir of the project - directory which contains code/ and data/ subdirs
    filepath = 'data/171129_A1R1_seq_n_assign/sequence/A1R1_1-100_1-letter.aa';
    sequence_filepath = fullfile(project_root, filepath);        

else
    
end

%% Import data
%===========================

seq.aa = importdata(sequence_filepath);
seq.num = [1:numel(seq.aa)]';
num_string = arrayfun(@num2str, seq.num, 'unif', 0);
seq.aa = cellfun(@(x,y) sprintf('%s%s',x,y), seq.aa, num_string, 'unif', 0);

end % main func
