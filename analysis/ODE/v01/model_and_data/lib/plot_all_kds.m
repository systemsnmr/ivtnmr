function plot_all_kds(use_log_scale)
% Loads all fitted Kd data. Visualizes in different forms.

% v01 2017-08-22
%   - init
%   - added error plotting

% TODO:
%  - options: log scale

% Flags to run different parts of code

global FIG; 
global SB;

%% Params to run "locally"
%===========================
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);    
    use_log_scale = 1;
else
    % ...
end

global DATA_IVTNMR;

if exist('./lib','dir')
    datasave_folder = strcat(DATA_IVTNMR, '/'); % if running in root (runall)
else
    datasave_folder = strcat('../', DATA_IVTNMR, '/'); % if running from lib directory
end


%% Load
%=========
load( fullfile(datasave_folder, 'kd_arr.mat') ,'kd_arr');
d = kd_arr;

% Mapping function
aas = containers.Map({'21', '22', '33', '56', '75'},...
                {'L21','S22','H33','G56','R75'});
            
            
%% Analyze results (plots & outputs)
%=======================================
peaks = unique(d.peak); % unique peaks
sets = unique(d.sname); % unique names
rnas = unique(d.rna);

n_peaks = numel(peaks);
n_sets = numel(sets);
n_rnas = numel(rnas);


%% Plot hist
%===============
rows = 2;
columns = n_peaks;
sbSide = 150; sbWidth = sbSide*1.5; sbHeight = sbSide;
fSize = 11;
plotLW = 1;

if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
FIG=FIG+1;
figure(FIG); set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 300 columns*sbWidth rows*sbHeight]);

SB=0;
sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

for p = 1:n_peaks
    SB=SB+1; sb(SB);
    
    subd = d(strcmp(d.peak, peaks(p)), :);
    
    y = subd.kd_mM;

    snames = subd.sname;
    x = cellfun(@(x) x(3:end-1), snames, 'unif', 0);
    
    sd = subd.kd_sd;
        
%     rna_names = unique(subd.rna);
%     x = rna_names;
    
%     semilogy(y,'o');
    plot(y,'o-');
    axis tight;
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    errorbar(y, sd, 'o-', 'LineWidth', plotLW);
    if use_log_scale
        set(gca,'YScale','log');
    end
    
    xlim(x_lim);
    ylim(y_lim);
    
    set(gca, 'XTick', 1:numel(x), 'XTickLabel', x);
    title(aas(peaks{p}));
    
    if p == 1
        ylabel('Kd [mM]');
    end

    
%     h(h.protein=='aroG' & strcmp(h.assign,'HIS'),{'FSA'}) = ...
    
end

disp('');
    
end