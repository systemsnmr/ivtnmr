function aborts_n_backfolding_02c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to design optimal RNA sequences for IVT / IVTNMR. Allows to do two
% things:
% 1) Algorithmically design abortive sequences (5'overhangs) - depending
% on filtering criteria - e.g. remove all 'GG's, etc.
% 2) Check that aborts do not back-fold on the main RNA sequence.

% For the second task - can also provide own list of 5'overhang sequences
% (i.e. do not need to run the algorithmically-designed-aborts part.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% More specifically one can:
% - Select for presence/absence of specific nucleotide patterns.
% - Remove sequences with too strong dimers in aborts or ones forming 
%   multiple hairpin species.
% ------
% - Dimer calculations are done based on DNA version of the sequence,
%   with an assumption that to a large-part energy of base-pairing is 
%   similar in DNA and RNA.
% - Calculations are rough - because neglect nearest-neighbor effects.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version history
% v01 - init, based on nUnD_v1_H68_170403 (nUnD_RNA_take3)
% v02 - added functions to algorithmically generate and filter abort seqs
%     - added purine-content filter (seqs, function which calcs and filters based on purines % (at least X%))
%     - added stable dimers filter (seqs, dG cutoff)
%     - added multi-hairpin filter (seqs)
% v02c 2018-09-17
%      - adds flag for enabling / disabling aborts check (run_property_check_for_aborts)
%      - adds more notes / documentation

% TODO:
%  - Add function to remove seqs with > 2 dimers! (hairpin & 5' abort)
%  - Add nearest-neighbor
%  - Add gap penalty
%  - Make a separate script - algorithmic generation of aborts

% @Yaroslav Nikolaev. 2017-..

% Settings
%====================
run_algorithmic_design_of_aborts = 1;

global DEBUG; DEBUG = 0;

% setting default values for the oligo analysis bioinformatic package.
global LENGTH_OF_DIMERS HPBaseValue HPLoopValue; 
LENGTH_OF_DIMERS = 3; % length=2 automatically impossible if cutoff larger than ~3 kcal/mol
HPBaseValue = 4; % Min paired bases at the neck of the hairpin. Default = 4.
HPLoopValue = 2; % Min number of bases that form the loop of a hairpin. Default = 2.

cutoff_dG = -9; % kcal/mol - if want to filter dimers by stability
cutoff_hairpin = 1; % max allowed hairpins in the sequence

run_property_check_for_aborts = 1; % disable if aborts are <= 2 nt

%% Algorithmic generation of aborts mixed with user input
%=====================================================================

%%% If want to check specific aborts sequences for back-folding onto the
%%% main RNA - enter them as a list here:
%=====================================================================
aborts = {...
'GGGGCCCC'
'GCACACCA'
'GCACCACA'
'GCCACACA'
'GCCGCACA'
    };

main_rnas = {...
%     'GCC CAG TGC TCT GAA TGT CAA AGT GAA GAA ATT CAA CCA AGC GCG GGT AAA CGG C'
%     'Gguaaccgggaaugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugg'
%     'GGGAUCAAUU CUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUG'
    'aaugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugcugg'
    'CUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUG'
    };


%%% If want to run algorithmic design of aborts, before checking them with
%%% the main sequence - put the filtering rules here:
%=====================================================================
if run_algorithmic_design_of_aborts
    seqs = generate_all_variants('AGC',8); % generate initial variants: nucleotides, length

    % Below filters are working based on REGEX matching:
    % For regex tests - check https://regex101.com/
    % seqs = filter_seqs('remove','[AG][AG]',all_variants); % removes all purine doublets    
    % seqs = filter_seqs('remove','[AG][AG]',all_variants); % removes all purine doublets
    seqs = filter_seqs('remove','GC',seqs); % removes all GCs
    seqs = filter_seqs('keep','^G',seqs); % keep only those starting with G
    seqs = filter_min_purine_percent(seqs, 30); % leave seqs with at least 30% A/G
    seqs = remove_seqs_with_stable_dimers(seqs, cutoff_dG); % calculates possible dimers, and removes sequences with stable dimers
    aborts = seqs;
    
    % Other examples of filtering:
%     seqs = filter_seqs('remove','[CU]{3,7}',seqs);  % remove stretches of 3 to 7 pyrimidines (C/U)    
%     seqs = filter_seqs('remove','T..T',seqs); % removes all U*U wobble - two Uridines spaced with two nucleotides
%    seqs = filter_seqs('remove','G$',seqs); % keep only those not ending with G
    
end


%% User-independent part of script
%==================================

% Process and display the input
aborts = rna2dna(aborts); % convert to dna
main_rnas = upper(main_rnas); % convert to upper case
main_rnas = regexprep(main_rnas,'[^\w'']',''); % remove spaces
main_rnas = rna2dna(main_rnas); % convert to dna

n_aborts = numel(aborts);
n_main_rnas = numel(main_rnas);

fprintf(1,'\n== Abortive sequences\n');
fprintf(1,'=====================\n');
% show aborts only if they are manually designed, or are not too many
if run_algorithmic_design_of_aborts && (n_aborts > 20)
    fprintf(1,'.. %i sequences total (not showing all)\n', n_aborts);
else
    disp(aborts); 
end

fprintf(1,'\n== Main RNA sequences\n');
fprintf(1,'=====================\n');
disp(main_rnas);


% Check properties of Aborts sequences
%==============================================

if run_property_check_for_aborts
    [abort_prop, abort_prop_extra] = calc_oligoprop(aborts); % determine properties
% aborts = remove_seqs_with_stable_dimers(aborts, cutoff_dG);
end

if ~run_algorithmic_design_of_aborts
    fprintf(1,'\n== Dimers (>= %i bp) in aborts\n',LENGTH_OF_DIMERS);
    fprintf(1,'================================\n');
    for i=1:n_aborts
        fprintf(1,'Seq = %s', aborts{i});
        if run_property_check_for_aborts && ~isempty(abort_prop(i).Dimers)
            abort_prop(i).Dimers
        else
            fprintf(1,' - None\n');
        end
    end
end


%% Check properties of FULL RNA sequence
%=======================================
% - add RNA sequences to aborts
% - calc oligo props
% - remove all with more than 1 hairpin
% - remove all with more than 2 dimers (self-anneal hairpin + 5'abort)
% - visually check the remaining ones

fprintf(1,'\n\nMaking all %i sequence combinations (Aborts)x(RNA)...\n\n', n_aborts*n_main_rnas);
% make all combinations of aborts+main_seq
[a, b] = ndgrid(1:numel(aborts),1:numel(main_rnas));
seq_full = strcat(aborts(a(:)), main_rnas(b(:)));
n_full = size(seq_full,1); % number of sequences

% Get properties of full sequences
[s2_prop, s2_prop_extra] = calc_oligoprop(seq_full);

if ~run_algorithmic_design_of_aborts
    fprintf('\n== Multi-hairpins (>= %i bp) in full RNAs\n',HPBaseValue);
    fprintf('==========================================\n');
    for i=1:n_full
        fprintf(1,'Seq # %i', i);    
    %     s2_prop(i).Hairpins    
        if ~isempty(s2_prop(i).Hairpins)
            s2_prop(i).Hairpins
        else
            fprintf(1,' - None\n');
        end
    end
end

% Remove sequences with > n hairpins
seq_full = remove_seqs_with_many_hairpins(seq_full, cutoff_hairpin);

% Remove sequences with > 2 dimers
% TODO - need to implement this function
% seq_full = remove_seqs_with_many_hairpins(seq_full, cutoff_hairpin);

% To check how many dimers:
% arrayfun(@(x) size(x.Dimers,1), s_prop, 'unif', 0)

fprintf('\n== Final seqs\n');
fprintf('===============\n');

if size(seq_full,1) > 0
    seq_full
%     s2_prop.Hairpins
%     s2_prop.Dimers
end

end % main func


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Generate sequence all possible variants
%================================================
function seqs = generate_all_variants(allowed_nucleotides,length)
n = length;
n_letters = numel(allowed_nucleotides);

x = 1:n_letters;

X = cell(1, n);
[X{:}] = ndgrid(x);
X = X(end : -1 : 1); 
y = cat(n+1, X{:});
y = reshape(y, [n_letters^n, n]);

seqs = cellstr( allowed_nucleotides(y) ); % generate and convert to cell array
n_seq = numel(seqs);

fprintf(1,'= All possible %i-nt seqs, from %i nucleotides (%s) : %i\n', n, n_letters, allowed_nucleotides, n_seq);
end


%% Filter sequences
%================================================
function seqs = filter_seqs(filter_flag,pattern,seqs)
n_seq0 = numel(seqs); % save number of starting sequences
idx = cellfun('isempty', regexp(seqs,pattern)); % get indexes of all sequences which match the pattern

if strcmp(filter_flag,'keep')
    idx = ~idx; % invert the indexes if want to keep, not remove
elseif ~strcmp(filter_flag,'remove')
    fprintf(1,'= %s - wrong filter_flag. Allowed two options: keep, remove.\n', filter_flag);
end

seqs = seqs(idx); % filter sequences based on the indexes found above
n_seq = numel(seqs); % calc final # of seq.

fprintf(1,'= Took %i/%i sequences by %s-ing %s\n', n_seq, n_seq0, filter_flag, pattern);
end



%% Set minimum purine fraction (should help the yield)
%================================================
function seqs = filter_min_purine_percent(seqs,cutoff_GA)
n_seq0 = numel(seqs); % starting # of seq

numSeq = cellfun(@(x) double(nt2int(x)), seqs, 'unif', 0); % let-to-numbers: % A => 1 , C => 2, G => 3,  T(U) => 4    
baseNum = cellfun(@(x) [sum(x == 1) sum(x == 2) sum(x == 3) sum(x == 4) sum(x == 15)], numSeq, 'unif', 0); % calc number matches
GA = cellfun(@(x,y) 100*((x(1) + x(3)) / length(y)), baseNum, numSeq); % calc GA content

min_GA_idx = GA >= cutoff_GA;
seqs = seqs(min_GA_idx);
n_seq = numel(seqs); % final # of seq

fprintf(1,'= Kept %i/%i seqs with mininum %u%% of purines \n', n_seq, n_seq0, cutoff_GA);
end



%% Remove sequences with stable dimers
%=====================================
function seqs = remove_seqs_with_stable_dimers(seqs, cutoff_dG)
n_seq0 = numel(seqs); % starting # of seq

[s_prop, s_prop_extra] = calc_oligoprop(seqs);

for i=1:n_seq0
    if numel(s_prop(i).Dimers) >= 1 % if there are ANY dimers
        if sum(s_prop_extra(i).dimers_dG < cutoff_dG)
            s_prop_extra(i).dimer_dG_above_cutoff = 1;
        end
    else
        % need to have at least empty field to avoid errors in below test.
        s_prop_extra(i).dimer_dG_above_cutoff = '';
    end
end

% Remove sequences which have dimers too stable
no_stable_dimers = cellfun(@isempty, arrayfun(@(x) x.dimer_dG_above_cutoff, s_prop_extra, 'unif' , 0));
seqs = seqs(no_stable_dimers);

n_seq = numel(seqs); % final # of seq

fprintf(1,'= Kept %i/%i seqs with dimers weaker than dG = %d \n', n_seq, n_seq0, cutoff_dG);
end



%% Remove sequences with too many hairpins
%=====================================
function seqs = remove_seqs_with_many_hairpins(seqs, cutoff_hairpin)
n_seq0 = numel(seqs); % starting # of seq

[s_prop, ~] = calc_oligoprop(seqs);

% Find and remove sequences with too many hairpins
more_than_n_hairpin = ~cellfun(@isempty, arrayfun(@(x) find(size(x.Hairpins,1) > cutoff_hairpin), s_prop, 'unif', 0));
seqs(more_than_n_hairpin) = [];

n_seq = numel(seqs); % final # of seq

fprintf(1,'= Kept %i/%i seqs with <= %d hairpins \n', n_seq, n_seq0, cutoff_hairpin);
end



%% Calculate oligo properties, including stability of dimers
%=============================================================
function [seqs_prop, seqs_prop_extra] = calc_oligoprop(seqs)
n_seq = numel(seqs);

global DEBUG;
global LENGTH_OF_DIMERS HPBaseValue HPLoopValue;

fprintf(1,'\n= Calculated oligo properties for %i seqs\n',n_seq);

% Data structure to keep basic results
seqs_prop(n_seq,1) = struct('GC',[],'GCdelta',[],'Hairpins',[],'Dimers',[],...
    'MolWeight',[],'MolWeightdelta',[],'Tm',[],...
    'Tmdelta',[],'Thermo',[],'Thermodelta',[]);

% Data structure to keep extra results
seqs_prop_extra(n_seq,1) = struct();

for i=1:n_seq
    % calc basic properties
    seqs_prop(i) = oligoprop(seqs{i}, 'Dimerlength',LENGTH_OF_DIMERS,'HPBase',HPBaseValue,'HPLoop',HPLoopValue);
    
    % calc stability of dimers - QUICK AND DIRTY (not Nearest neighbor!)
    seqs_prop_extra(i).n_dimers = size(seqs_prop(i).Dimers, 1);
    if seqs_prop_extra(i).n_dimers
        if DEBUG
            fprintf(1,'Calculating dG for %i dimers in seq #%i:\n', seqs_prop_extra(i).n_dimers, i);
        end
        seqs_prop_extra(i).dimers_dG = nan(seqs_prop_extra(i).n_dimers, 1);
        AU = regexp(cellstr(seqs_prop(i).Dimers),'[ATU]');
        GC = regexp(cellstr(seqs_prop(i).Dimers),'[GC]');

        %%% TODO: implement inclusion of GAP PENALTY - for when dimers are
        %%% from fragments:
        n_gaps = cell2mat(cellfun(@(x) numel(x), regexp(cellstr(seqs_prop(i).Dimers),'[AGCT][agct]|[agct][AGCT]'), 'un', 0));
        % 0 hit in the above: dimer fits both ends - use the below gap
        % 1 hit: one half of dimer fits - use the below penalty once
        % 2 hits?        
        % - find which penalty to use
        % - check if dimer is ending or starting the sequence (first or
        % last capital) - and if yes, then ...
        
        % potentially need nearest neighbor (https://www.ncbi.nlm.nih.gov/pubmed/10329189)
        % and nearest-neighbor for DNA:
        % https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
        % but for our quick-n-dirty - enough 1-to-1.

        % Based on: http://eu.idtdna.com/calc/analyzer
        % for RNA molecule!!
        % 
        % For GC pairs (hetero-dimer, self-dimer gives same)
        % 2xGC = 3.07/2 = 1.535 -- self-dimer
        % 2xG = 3.07/2 = 1.5
        % 3xG = 6.1/3 = 2.03
        % 4xG = 9.21/4 = 2.302
        % 100XG = 303.85/100 = 3.04
        % > I.e. first (initiating) pair is not even calculated!
        % 
        % For AU pairs
        % 2xA = 1.94/2 = 0.97
        % 100xA = 192.5/100 = 1.93
        % Using 3.04
        
        dG_AU = cellfun(@(x) size(x, 2), AU) .* -1.93; 
        dG_GC = cellfun(@(x) size(x, 2), GC) .* -3.04;
        seqs_prop_extra(i).dimers_dG = dG_GC + dG_AU - (-1.93-3.04)/2; % subtract average terminating energy
        if DEBUG
            disp(seqs_prop(i).Dimers);
            disp(seqs_prop_extra(i).dimers_dG);
        end
%         if sum(seqs_prop_extra(i).dimers_dG < cutoff_dG)
%             seqs_prop_extra(i).dimer_dG_above_cutoff = 1;
%         end
    end

end

end




