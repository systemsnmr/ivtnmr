function plot_H_N_HN(protein)

%% Params to run "locally"
%===========================
if nargin == 0
    protein = 'A1R1';
else
end

switch protein
    case 'A1R1'
        x_lim = [1 100];        
    case 'UP1'
        x_lim = [1 200];
    otherwise
        error('Unknown protein name %s', protein);
end

%% Import data
%===========================
datasave_folder = 'datasave';
datafile = fullfile(datasave_folder, sprintf('assign_%s.mat',protein));
load(datafile);

%% Plotting setup
%===========================

global FIG; 
global SB;

 % roughly - some calculation like this (ratio width:height 4:3)
 % - rows = ceil(n_plots / 4)
 % - columns = ceil(n_plots / 3)
 % - check size of screen - divide width(pixels) / columns - will give
 % rough sbSide to fit all into one screen
 
rows = 3;
columns = 1; % columns = n_peaks /
sbSide = 750;
sbWidth = sbSide;
sbHeight = sbSide/4;
fSize = 11;

if isempty(FIG), FIG=0; end; % if this is the first figure - set FIG index to zero
FIG=FIG+1;
figure(FIG); set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 300 columns*sbWidth rows*sbHeight]);

SB=0;
sb = @(x) subplot(rows,columns,x,'FontSize',fSize); % Helper function to draw subplots

%% Example of how to run kmeans here on first 10 peaks
%====================================
kmeans(assign.traj_HNcsp(1:10,:)',3) % need to transpose, otherwise kmeans is unhappy
disp('');

%% Plot data
%===========================
SB=SB+1; sb(SB);
bar(assign.ids, assign.traj_maxHdelta);
xlim(x_lim);
ylabel('Hdelta [ppm]');

SB=SB+1; sb(SB);
bar(assign.ids, assign.traj_maxNdelta);
xlim(x_lim);
ylabel('Ndelta [ppm]');

SB=SB+1; sb(SB);
bar(assign.ids, assign.traj_maxHNcsp);
xlim(x_lim);
ylabel('HNcsp [ppm]');

end