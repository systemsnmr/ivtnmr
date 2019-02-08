function plot_mean_kds(use_log_scale)

% Patched this file a bit on 20180821 - to visualize RNA-based Kds.

global FIG; 
global SB;

%% Params to run "locally"
%===========================
if nargin == 0    
    use_log_scale = 1;
    datafile = 'kd_arr.mat';        
else
    % ...
    datafile = 'kd_arr.mat';
end

global DATA_IVTNMR;

if exist('./lib','dir')
    datasave_folder = strcat(DATA_IVTNMR, '/'); % if running in root (runall)
else
    if exist('./sublib','dir')
        datasave_folder = strcat('../', DATA_IVTNMR, '/');
    else
        datasave_folder = strcat('../../', DATA_IVTNMR, '/'); % if running from sub-lib
    end
end

%% Load
%=========
load( fullfile(datasave_folder, datafile) ,'kd_arr');
d = kd_arr;

% Mapping function
aas = containers.Map({'21', '22', '33', '56', '75'},...
                {'L21','S22','H33','G56','R75'});

set_to_rna = containers.Map({'IN60', 'IN61', 'IN62', 'IN65', 'IN70', 'IN71', 'IN72', 'IN75'},...
                {'RNA0','SMN1','SMN2','EV2','RNA0','SMN1','SMN2','EV2'});
            

%% Filtering datasets
%======================
% drop_IN6x = 1; % A1R1 datasets
% if drop_IN6x
%     d = d(strcmp(d.prot,'UP1'),:);
% end
            
            
%% Analyze results (plots & outputs)
%=======================================
peaks = unique(d.peak); % unique peaks
sets = unique(d.sname); % unique names

% Weird QND way to change the variable in dataset :-)
if any(strcmp(d.rna,'SMN214'))
    n_entries_to_change = sum(strcmp(d.rna,'SMN214'));
    d(strcmp(d.rna,'SMN214'),'rna') = dataset({repmat({'SMN2'},n_entries_to_change,1), 'rna'});
end;

rnas = unique(d.rna);
% need to fix the order:
rnas = {'RNA0' 'SMN1' 'SMN2' 'EV2'};

n_peaks = numel(peaks);
n_sets = numel(sets);
n_rnas = numel(rnas);

            
%% Mean Kds
%=======================================
d = d(~strcmp(d.peak,'56'),:); % remove G56 peak (unspecific)
peaks = unique(d.peak); % unique peaks
n_peaks = numel(peaks);

kd_mean_list = nan(n_peaks,1);
kd_sd_list = nan(n_peaks,1);

% Now doing over RNAs - cuz each has few sets:
% for s = 1:n_sets
for s = 1:n_rnas
    % but comparing with dset name - easier
    subset = d(strcmp(d.rna, rnas{s} ),:);
    kd_mean_list(s) = mean(subset.kd_mM);
    kd_sd_list(s) = std(subset.kd_mM);
%     fprintf(1,'%s Kd = %.3f (%.3f) mM\n', aas(peaks{p}), kd_mean, kd_sd);
%     fprintf(1,'%s Kd = %.3f (%.3f) mM\n', sets{s}, kd_mean_list(s), kd_sd_list(s));
end

kd_mean_list
kd_sd_list

%% Plot hist
%===============
rows = 2;
columns = 2;
sbSide = 220; sbWidth = sbSide; sbHeight = sbSide;
fSize = 12;
plotLW = 2;


if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
FIG=FIG+1;
figure(FIG); set(figure(FIG), 'Color', repmat(0.97,1,3), 'Position', [0 300 columns*sbWidth rows*sbHeight]);

SB=0;
sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

SB=SB+1; sb(SB);

y = kd_mean_list;
sd = kd_sd_list;

%%% New ITC values - 20180813
%%% from:
% file:///Volumes/Data/yar/Dropbox/_eth2/data_MathModeling/170614_IN60_IN70/v12c_H33_R75_RNA_sc/lib/get_k_from_bootstrap_f_paper.m
ITC_vals = [1390.8 51.3 47.4 5.1];
ITC_errs = [330.78 2.53 19.74 1.9];
ref_y = ITC_vals ./ 1000;
ref_sd = ITC_errs ./ 1000; 

% xlab = cellfun(@(x) x(3:end-1), sets, 'un', 0);
in_num = cellfun(@(x) x(1:4), sets, 'un', 0);
xlab = cellfun(@(x) set_to_rna(x), in_num, 'un', 0);
    
    semilogy(y,'o');
%     plot(y,'o-');
    axis tight;
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    errorbar(y, sd, 'o--', 'LineWidth', plotLW);
    hold on;
    plot(ref_y,'or--', 'LineWidth', plotLW);
    errorbar(ref_y, ref_sd, 'or', 'LineWidth', plotLW);
    
    if use_log_scale
        set(gca,'YScale','log');
    end
    
    xlim(x_lim);
%     ylim(y_lim);
%     ylim([0.5e-4 y_lim(2)]);
%     ylim([0.5e-4 1e2]);
%     ylim([0.5e-4 1.5e1]);
    y_lim = [0 max([kd_mean_list' ITC_vals])*1.3];

    set(gca, 'XTick', 1:numel(xlab), 'XTickLabel', xlab);
    title('Mean Kd');
    
    ylabel('Kd [mM]');
    
end
