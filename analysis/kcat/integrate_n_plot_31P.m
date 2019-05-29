function integrate_n_plot_31P
% integrate_n_plot_31P: reads 31P NMR data (multiple datasets possible); integrates
% peaks; calculates RNA synthesis rates (based on NTP consumption) and plots
% them. Saves intermediate data into *.mat files. And final results - into
% CSV and EPS figures. 

% rawList = importdata('shell_cat_all_txts.txt');
% collect_csp_data_02 - for IVTNMR analysis - uses time0 too
% collect_csp_data_02_5Pz - to analyse HNCO data
%%% 
% later based on: collect_csp_data_02_5Pz_FIG_v01.m
% later based on: r_collect_5pz_CSP_IN99.m
% 
global FIG; % 

% addpath('/Volumes/Data/yar/Dropbox/Programming/Matlab/tests');
% addpath(genpath('/local/home/allainkurs1203/Desktop/BK2019/code/lib/rbnmr'));
% addpath(genpath('./lib'));
addpath(genpath(fullfile(pwd,'lib')));

%% Settings
%===========================
flag_reload_raw = 1;
flag_save_to_scratch = 1; % inside the reload_raw
flag_load_from_scratch = 1;

flag_export_NTP_and_RNA_conc = 1;
datasave_dir = 'datasave';

% scratchdir = '/scratch/polyx';
scratchdir = 'datasave';
scratchfile_prefix = strcat(mfilename,'_2')

%% Data variables
%=======================================================
%%% Dir with NMR data
CURDIR = '/Users/yar/_eth2/data_NMR/spectra';
% CURDIR = '/local/home/allainkurs1203/Desktop/BK2019/nmr_spectra';

init_NTP_conc = 20; % mM
rna_length = 20*0.3+6*0.7; % average between 30% full RNA, and 70% ~6nt aborts.
x_end = 1100;

dset_list = {   
%     '190415_IN115a_pA20_co-P50N525_100uM_298K_600'
    '190314_IN107i_pA9_free_ATP_BK'
    '190306_IN107a_pA9_free_BK'
    '190311_IN111a_pA20_free_BK'
    };

dset_names = {...
%     'A20 (ATP) + P50N525'; 
    'A9 (ATP)';
    'A9 (4 NTPs)';
    'A20 (4 NTPs)';
    };

intrng = [-10.4 -11.2];

plotLW = 2;

colors = [
    [1, 0, 0]; % red
    [0, 204/255, 0]; % green
    [0, 0, 1]; % blue
    [0,0,0]; % black
    %[192/255,192/255,192/255]; % grey
    %[88/255,0,0]; % dark red
    %[0, 102/255, 0]; % dark green
    %[0, 0, 0.5]; % dark blue
    %[48/255,48/255,48/255]; % dark grey
    %[160/255,0,0]; % light red
    %[0, 255/255, 255/255]; %turquoise
    %[102/255, 255/255, 102/255]; % light green
    %[0.4, 0.4, 1]; % light blue
    %[200/255,200/255,200/255]; %light green
    ];


%% Preprocessing
%===========================
dset_time0 = cellfun(@(x) get_time0(x, CURDIR), dset_list); % uses helper func (see EOF)
n_sets = numel(dset_list);
dset_exp_31P = cell(n_sets,1);

for i=1:n_sets
    dset_exp_31P{i} = cell2mat(arrayfun(@(x) str2num(x.name), dir(fullfile(CURDIR, dset_list{i},'5*')), 'un', 0))';
    % Filter out only experiments > 5000 - i.e. drop any 5, 50, 500
    dset_exp_31P{i} = dset_exp_31P{i}(dset_exp_31P{i} >= 5000);
end; clear i;

% Create dirs if don't exist yet.
if ~exist( datasave_dir ,'dir'), mkdir(datasave_dir), end;
if ~exist( scratchdir ,'dir'), mkdir(scratchdir), end;


%% Integrate 31P
%===========================
% could have just integrated!

if flag_reload_raw
aNTP_conc = cell(n_sets,1);
rna_conc = cell(n_sets,1);
time_31p = cell(n_sets,1);

for i=1:n_sets
    dset_name = dset_list{i};
    expnos = dset_exp_31P{i};
    nmr_data = arrayfun(@(x) fullfile(CURDIR, dset_name, num2str(x), 'pdata/1/1r'), expnos, 'un', 0);
    integrals_31p = integrate1D(nmr_data,intrng); % uses RBNMR
    aNTP_conc{i} = (integrals_31p./integrals_31p(1)).*init_NTP_conc; % mM
    NTP_consumed = init_NTP_conc-aNTP_conc{i};
    rna_conc{i} = NTP_consumed./rna_length;
    
    % Get time vector
    time_31p{i} = getNMRTime_02(fullfile(CURDIR, dset_name), expnos, dset_time0{i});
%     % Interpolate
%     time_hn = getNMRTime(fullfile(nmr_data_path,dset_name), 4000:4053, time0);
%     rna_conc_hn_time = interp1(time_31p, rna_conc, time_hn, 'linear', 'extrap');    

end

% Assemble structure for save:
d.dset_list =      dset_list;
d.dset_time0 =     dset_time0;
d.dset_exp_31P =   dset_exp_31P;
d.aNTP_conc =      aNTP_conc;
d.rna_conc =       rna_conc;
d.time_31p =       time_31p;
d.names = dset_names;

if flag_save_to_scratch; save( fullfile(scratchdir, sprintf('%s.mat', scratchfile_prefix)), 'd'); end;

end; % flag_reload_raw

if flag_load_from_scratch; load(fullfile(scratchdir, sprintf('%s.mat', scratchfile_prefix))); end;

disp('');

%% Estimate kcats
%====================================
kcat = cell(n_sets,1);
for i=1:n_sets
    f = fit(d.time_31p{i}, d.aNTP_conc{i}, 'exp1');
    % evaluate cfit object to extract data.    
    kcat{i} = abs(f.b) / 60 / (33*1e-6); % /60 min>sec. *(33*1e-6) - normalize by template (nM in mM units).       

    if flag_export_NTP_and_RNA_conc
        save_csv_data(datasave_dir, sprintf('RNA_%s.csv',dset_names{i}), [d.time_31p{i}(:), d.rna_conc{i}(:)]);
        save_csv_data(datasave_dir, sprintf('NTP_%s.csv',dset_names{i}), [d.time_31p{i}(:), d.aNTP_conc{i}(:)]);
    end
end; clear i;

legend_with_kcat = cellfun(@(x,y) sprintf('%s (kcat=%.2f s-1)', x, y), dset_names, kcat, 'un', 0);       

%% Plot
%====================================
rows = 2;
columns = n_sets+2;
sbSide = 150;
sbWidth = sbSide;
sbHeight = sbSide;
fSize = 11;

FIG = 1;
fig_handle = figure(FIG); set(figure(FIG), 'Color', ones(1,3), 'Position', [0 300 columns*sbWidth rows*sbHeight]);

SB=0;
sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

% colors = lines(n_sets);
SB=SB+1;
sb([SB]);

x_end = x_end/60;

for iSet=1:n_sets
    x_dat = d.time_31p{iSet}./60;
    y_dat = d.rna_conc{iSet};
    plot(x_dat, y_dat, '-', 'Color', colors(iSet,:), 'LineWidth', plotLW);
    hold on;    
end
axis tight;
xlim([0 x_end]);

ylabel('RNA [mM] (via 31P)');
xlabel('time [hr]');

axP = get(gca,'Position'); 
% legend(dset_names, 'Location', 'eo');
legend(legend_with_kcat, 'Location', 'eo');
set(gca, 'Position', axP);

file_name = sprintf('%s_%s.eps', mfilename, datestr(now,'yymmdd'));
save_figure(fig_handle,file_name);

end % main func


%% Helper functions
%==============================
function time0 = get_time0(dset, CURDIR)
    notes_path = fullfile(CURDIR,dset,'notes.txt');
    filetext = fileread(notes_path);
    pattern = '<time0>(.*?)</time0>';
    time0 = regexp(filetext,pattern,'tokens');
    assert(numel(time0) == 1, 'Abort: more than one time0 matches found in notes.txt.');    
    time0 = time0{:}; % convert from cell to a string
    clear filetext;
end

%% Export data in csv
%==============================
function save_csv_data(dir, filename, data)

%     filename = fullfile(scratch_dir, 'testdata.csv');
%     A = [12.7 5.02 -98 63.9 0 -.2 56];
%     csvwrite(filename,A)

%     export_full_name2 = fullfile(export_csv_path, sprintf('%s_data.csv',export_name));
%     csvwrite(export_full_name2,[xd(:),yd(:)]);

    csvwrite(fullfile(dir,filename), data);
end

function save_figure(fig_handle,file_name)
%% Save figure as pdf
%=============================================================
    set(fig_handle,'Units','Inches');
    pos = get(fig_handle,'Position');
    set(fig_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    if ~exist( 'figure_images' ,'dir'), mkdir('figure_images'), end;
%    print(fig_handle,fullfile('figure_images',file_name),'-dpdf','-r1200')
	print('-depsc', fullfile('figure_images',file_name));

% CONVENIENT TO USE with mfilename variable - which returns the name of the current file!
% 
% -r option sets the resolution (also for the VECTOR IMAGES!) -- set to 300-1200 if contour approximation is bad.
% -r0
% -r300
% -r1200
end

