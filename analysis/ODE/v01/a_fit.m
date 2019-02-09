function fit_data = a_fit(data_name, optns)
close all;

modelSubFolder = 'model_and_data';

if nargin == 0
    data_name = 'IN70a_intern_HN';    
%     data_name = 'IN75a_intern_HN';    
%     data_name = 'IN72b_intern_HN';        

    optns.HN_peak_for_fit = 'all';
    
    %%% Use these options if (1) protein is not saturated with RNA during reaction time
    %%% and (2) you don't have a good estimate for maximum chemical
    %%% shift perturbations of specific residues at saturation. Then can make an assumption
    %%% that maxcsp is at least 2-fold the maximum observed shift.
    optns.sets_for_maxcsp = {'IN70a' 'IN70b' 'IN70c' 'IN71a' 'IN71b' 'IN71c'...
        'IN72a' 'IN72b' 'IN72c' 'IN72d' 'IN75a' 'IN75b' 'IN75h' 'IN75z'};
    optns.peaks_for_maxcsp = []; % if empty - will use all peaks
    %%% if you have a linear CSP in the course of reaction, assuming that maximum
    %%% shift is < than 1.9 is unrealistic. maxcsp_scaling = 2 seems a reasonable
    %%% lower limit.
    optns.maxcsp_scaling = 2;
    
    % This overrides the above settings and uses fixed CSP values -
    % need to be measured experimentally/independently.
    optns.maxcsp_fixed = 1;
    optns.maxcsp_list = containers.Map({...
            'H33','R75'},...
            [0.083, 0.046]);

    optns.rna_fraction = 0.3; % will be used below to normalize RNA/Aborts fractions
        
    optns.T7DNAcomplex_concentration = 33*1e-6; % mM
    optns.Mg_total = 24; % mM
    optns.NTP_initial_conc = 20; % mM
    
    % Bootstrap stuff (not used anymore!)
    optns.useBootstrap = 0; % if set to 0 - runs regular fit with full dataset
    optns.nboot = 200; % <= 7/minute
    
end 

%%% TODO - fix this - w/o using GLOBAL variables (make as passed parameters
%%% for the auxiliary functions)
global a_fit_T7DNAcomplex_concentration; a_fit_T7DNAcomplex_concentration = optns.T7DNAcomplex_concentration;
global a_fit_Mg_total; a_fit_Mg_total = optns.Mg_total;
global a_fit_NTP_initial_conc; a_fit_NTP_initial_conc = optns.NTP_initial_conc;


%%% Get maxcsp for peaks (move out to the multi_fit file in future?)
%=========================
addpath(genpath(fullfile('.', modelSubFolder, 'lib')));

if optns.maxcsp_fixed
    maxcsp_list = optns.maxcsp_list;
    optns.maxcsp_scaling = 1;
else
    base_path = fullfile(modelSubFolder, 'data');
    dset_names = optns.sets_for_maxcsp;
    dset_suffix = data_name(7:end);
    nmr_data_paths = cellfun(@(x) fullfile(base_path, sprintf('data_nmr_selected_%s_%s.mat', x, dset_suffix)), dset_names, 'un', 0);    
    maxcsp_list = get_maxcsp_for_peaks(nmr_data_paths, optns.peaks_for_maxcsp);
end

%=========================
f_load_data_and_model = 1;
f_use_time_below_hours = 24; % dropping several-day datasets
f_plot_solution = 1;
f_plot_residuals = 1; % works only when optns.useBootstrap = 0
f_show_NRMSD_etc_for_each_plot = 0;
f_disp_constants_n_stats = 1; % works only when optns.useBootstrap = 0

if nargin == 0, tic, end; % internal timing

%% Settings n code sections to run
%=========================
global FIG; FIG=0;
SB=0;

sbSide = 140;
sbWidth = 160;
sbHeight = sbSide;
fSize = 12;
flag_common_yscale = 0;

% Colors
%=======================
cWhite  = [1, 1, 1];
cBlack  = [0, 0, 0];     cDarkGrey = [0.4,0.4,0.4];     cLightGrey   = [0.6,0.6,0.6];
cBlue   = [0, 0, 1];     cDarkBlue = [0.1, 0.1, 0.56];  cLightBlue   = [0.3, 0.4, 1];
cRed    = [1, 0, 0];     cDarkRed    = [0.56, 0.17, 0]; cLightRed    = [1, 0.4, 0.2];  
cGreen  = [0, 1, 0];     cDarkGreen  = [0.17, 0.56, 0]; cLightGreen  = [0.2,0.6,0.2];
cOrange = [1,0.45,0.01]; cDarkOrange = [0.8,0.3,0.01];   cLightOrange = [1,0.8,0.6];


%% Load data and model
%=======================
datafile = fullfile(modelSubFolder, 'data_for_fit', sprintf('data_nmr_selected_%s.mat', data_name) );

if f_load_data_and_model

mdf_file_name = sprintf('s_%s/xIVT.mdf', data_name(1:4));
M = model4fit( fullfile(modelSubFolder, mdf_file_name) );
pin = struct(); % paramters initial

%% Setting initial parameters
%==============================
% General considerations:
% Diff-limits: from 1e5 - 1e10 M-1 s-1 (From geom.specificity of the encounter to electrostatically-enhanced).
% Diff-lim on-rates: diff-lim * 60(min-1)*10^-3(mM-1)*1e-1(1 o/m slower - common expt vs theory for diffusion)

%%% Define non-optimized parameters (when just keep some rxns fixed)
%=====================================================================
% The parameters will be removed from optimization vector later in the script.
% non_opt_param_names = {'kon' 'koff'}; % drop by name
% non_opt_param_names = {};
if exist('non_opt_param_names', 'var') && ~isempty(non_opt_param_names)
    [~, non_opt_param] = ismember(non_opt_param_names, M.par_names);
else
    non_opt_param = [];
end;
% non_opt_param = [2 3 4]; % drop non-optimized params by number

%%% Below many more parameters are defined than the ones actually used in
%%% specific models -- this was to accomodate different model complexities,
%%% e.g. including RNAP T7 degradation, RNAP T7 binding to DNA template,
%%% etc..

pin.kcat        = [1.94e-04  	, 0,   100    ]; % for ~30 nt RNA [min-1]
pin.kabortive   = [3.85e-04  	, 0,   100    ]; %
pin.kfold       = [60  			, 0,   1e5   ]; %
pin.kunfold     = [1  			, 0,   1e5    ]; %
pin.kcatPPi     = [1  			, 0.9, 1.1  ]; %
pin.kdegNTP     = [4.2280e-06  	, 0,   100    ]; % 
pin.kprecip     = [19.96		, 0,   100    ]; % 
pin.kdissolve   = [3.08 	   	, 0,   1000    ]; %
pin.kmgon       = [60*10^4 		, 0,   60*1e8 ]; % Eletrost.sm.mol: 1e10(M-1 s-1)* 60(min-1)*10^-3(mM-1)*1e-1(slower)
KdMgNTP = 0.3; % 0.1-0.5 mM
pin.kmgoff      = [60*10^3      , 0, 60*1e7]; % koff = Kd*kon. Upper limit - 2o/m the pinit

pin.kTDon       = [60*10^8*10^-3, 0,   1e9   ]; % Medium molecule: 1e8(M-1 s-1)* 60(min-1)* 1e-3(mM-1)* 1e-1(slower)
KdTD = 16* 1e-6; % 16nM* mM scale
pin.kTDoff      = [pin.kTDon(1)*KdTD, 0,   1e8   ]; % koff = Kd*kon. Upper limit - 2o/m the pinit

pin.kdegT7      = [0.001		, 0,   1e5   ]; % n specific reason f g value. Expasy T7 half-life: ~10h in vivo. k = ln(2)/(t1/2) = 0.0012 min-1

pin.kon         = [4.0271e+05          , 0,   1e9   ]; % Medium molecule: 1e8(M-1 s-1)* 60(min-1)* 1e-3(mM-1)* 1e-1(slower)
KdUP1 = 50* 1e-3; % 50uM*mM scale
pin.koff        = [pin.kon(1)*KdUP1	, 0,   1e8   ]; % koff = Kd*kon. Upper limit - 2o/m the pinit. pin.kTDon(1)*KdUP1
pin.koff        = [2.6055e+06	, 0,   1e8   ]; % koff = Kd*kon. Upper limit - 2o/m the pinit. pin.kTDon(1)*KdUP1

pin.kAbon         = pin.kon;
pin.kAboff        = pin.koff;
pin.kAboff   = pin.kAboff*100; % Decrease Kd for abortives

pin.kNTPon         = pin.kon;
pin.kNTPoff        = pin.koff;
pin.kNTPoff   = pin.kNTPoff*100; % Decrease Kd for NTPs

% Adjust (narrow) boundaries for T7-DNA Kd
pin.kTDon([2 3])  = pin.kTDon(1) * [0.1 10];
pin.kTDoff([2 3]) = pin.kTDoff(1) * [0.1 10];

% Assign params
for i = 1:M.npar
  pn = M.par_names{i};
  [M.par_init(i), M.par_min(i), M.par_max(i)] = deal(...
      pin.(pn)(1), pin.(pn)(2), pin.(pn)(3)); 
end; clear i

%%% LOAD DATA
%=======================
d = load(datafile); % for data generation see create_data_str.m
data = d.dsel; % Get data matrix from the structure: {'PO4','NTP','RNAtotal','SMNfolded','Abortives','NDP','MgPO4','HSQC';}
clear d.dsel;

%%% Normalize HN Prot observable
%=======================
prot_idx = numel(data.names); % assuming protein is LAST

if ~strcmp(optns.HN_peak_for_fit,'all') % do below only if its not for ALL peaks

    [~, hn_peak_idx] = ismember(optns.HN_peak_for_fit, data.hn.names);

    % Convert ppm to protein saturation.
    % DIVIDING by maxcsp_scaling, because at maxcsp_scaling saturation should be 1.
    ppm_to_saturation_scale = 1 ./ maxcsp_list(optns.HN_peak_for_fit) ./ optns.maxcsp_scaling;
    prot_y = data.hn.mean(:,hn_peak_idx) .* ppm_to_saturation_scale;
    prot_s = data.hn.std(:,hn_peak_idx) .* ppm_to_saturation_scale;

    % Convert protein saturation to concentration units
    prot_y = prot_y .* data.hn.p_conc;
    prot_s = prot_s .* data.hn.p_conc;

    data.y(prot_idx,:) = prot_y';
    data.s(prot_idx,:) = prot_s';

    data.names{prot_idx} = sprintf('%s(%s)', optns.HN_peak_for_fit, data.sname);
else    
    
    % Convert ppm to protein saturation.
    % DIVIDING by maxcsp_scaling, because at maxcsp_scaling saturation should be 1.
    n_HN_peaks = numel(data.hn.names);
    
    % preallocate NANs:
    data.y = [data.y(1:end-1,:); nan(n_HN_peaks,size(data.y,2))];
    data.s = [data.s(1:end-1,:); nan(n_HN_peaks,size(data.s,2))];
    data.names(prot_idx) = []; % remove completely - will add at each step
        
    for iPeak=1:n_HN_peaks
        %%% TODO could get maxcsp in form of vector - and then vectorize
        %%% whole thing.
        ppm_to_saturation_scale = 1 ./ maxcsp_list(data.hn.names{iPeak}) ./ optns.maxcsp_scaling;
        prot_y = data.hn.mean(:,iPeak) .* ppm_to_saturation_scale;
        prot_s = data.hn.std(:,iPeak) .* ppm_to_saturation_scale;

        % Convert protein saturation to concentration units
        prot_y = prot_y .* data.hn.p_conc;
        prot_s = prot_s .* data.hn.p_conc;

        data.y(prot_idx+(iPeak-1),:) = prot_y';
        data.s(prot_idx+(iPeak-1),:) = prot_s';
        data.names = [data.names sprintf('%s(%s)', data.hn.names{iPeak}, data.sname)];
        
    end; clear iPeak;
end


%%% Trim time-vector
%=======================
if exist('f_use_time_below_hours','var') && ~isempty(f_use_time_below_hours) && (max(data.t)>(f_use_time_below_hours*60))
    n_points_below_cutoff = numel(data.t(data.t < (f_use_time_below_hours*60)));
    data.t = data.t(1:n_points_below_cutoff);
    data.y = data.y(:, 1:n_points_below_cutoff);
    data.s = data.s(:, 1:n_points_below_cutoff);        
    % now this is saved below as "data0" (d.dsel not used anymore)
%     d.dsel = data; % re-save data, cuz d.dsel is used later
end



%%% Calc RNA from NTP
%=======================
%%% In most of the setups - for a >10nt-long RNA, with 5 min per 1D31P
%%% spectrum - RNA concentration determined directly from RNA peak will
%%% have low accuracy. More accurate is to use NTPs as a proxy to calculate
%%% how much RNA was synthesized.
% convert_RNA_from_NTP(data, source_name, target_name, target_new_name, init_NTP_conc)
data0 = data;
data = convert_RNA_from_NTP(data, 'NTP', 'RNAtotal', 'RNAfromNTP', optns.NTP_initial_conc);

% Correct MgPO4 observable after the above RNA conc change
%=======================
MgPO4_idx = strcmp(data.names, 'MgPO4');
RNA_idx = strcmp(data.names, 'RNAfromNTP');
RNA_difference = ( data0.y(RNA_idx,:) - data.y(RNA_idx,:) );
data.y(MgPO4_idx,:) = data.y(MgPO4_idx,:) + RNA_difference;



%%% Set the measured RNA/Aborts fractions
%=======================
% RNAtotal_idx = find(strcmp(data.names, 'RNAtotal'));
RNA_idx = strcmp(data.names, 'RNAfromNTP');
Abortives_idx = strcmp(data.names, 'Abortives');

data.y(Abortives_idx,:) = data.y(RNA_idx,:) .* (1-optns.rna_fraction);
data.s(Abortives_idx,:) = data.s(RNA_idx,:) .* (1-optns.rna_fraction);
data.y(RNA_idx,:) = data.y(RNA_idx,:) .* optns.rna_fraction;
data.s(RNA_idx,:) = data.s(RNA_idx,:) .* optns.rna_fraction;
data.names{Abortives_idx} = sprintf('Aborts (%.0f%%)', (1-optns.rna_fraction)*100);


%%% Scale/calculate RNA from 31P atoms to actual molecules
%=======================
data.y(RNA_idx,:) = data.y(RNA_idx,:) ./ data.hn.rna_length;
data.s(RNA_idx,:) = data.s(RNA_idx,:) ./ data.hn.rna_length;

%%% TODO - automatically filter out only the stuff which is in the
%%% observables -- 20180227 - BUT THEN NEED TO KEEP IN MIND OBSERVABLE
%%% NAMEs?!
%=======================
% Current export
% 'PO4'    'NTP'    'RNAtotal'    'Abortives'    'NDP'    'MgPO4' 'PROT not def.'
% TODO - can automatically decide here which 
%%% Currently have to set these manually, if using only some of
%%% observables.
obs_data_to_keep = [1 2 3 4 5 6 7];
colors = [cLightGrey; cOrange; cBlack; cBlack; cOrange; cLightGrey; cRed]; % set for OBSERVABLES
%%% Manually define colors for individual observables above - if want to.
%%% Otherwise will be adapted below automatically.
n_obs = numel(obs_data_to_keep); 
n_colors = size(colors,1);
if n_colors > n_obs % remove extra colors
    colors = colors(1:n_obs,:);
elseif n_colors < n_obs % add extra colors
    colors = [colors; repmat(cLightGrey, n_obs-n_colors, 1)];
end;

% Adds more observables - to use multiple protein peaks in a global fit.
if strcmp(optns.HN_peak_for_fit,'all')
    for iPeak=1:n_HN_peaks-1
        obs_data_to_keep = [obs_data_to_keep obs_data_to_keep(end)+1];
        colors = [colors; cRed];
    end; clear iPeak;
end

% Adapt to the observables
if exist('obs_data_to_keep','var') && ~isempty(obs_data_to_keep)    
    data.names = data.names(obs_data_to_keep);
    data.y = data.y(obs_data_to_keep,:);
    data.s = data.s(obs_data_to_keep,:);
end

% create time span for ODE simulations
tspan = 0:2:round(data.t(end)); % endpoint in time vector

obs_names = data.names;


% Check scale of observables - important that residuals for different
% observables are not TOO different! Otherwise better transition to use
% PottersWheel and LOG-scale optimization!
max_min_scale = max(mean(data.y,2)) / min(mean(data.y,2));
if max_min_scale > 1000
    fprintf(1,'\n= max/min data mean (%1.f) is > 1000:\n consider using log-scale optimization (e.g. PottersWheel)\n\n', max_min_scale);
end;

%=======================================
%=======================================
%%% Adjusting the BNGL model initial params
%=======================================
%=======================================

%%% NDP: set initial concentration to their actual initial values
%=======================================
[~, NDP_model_idx] = ismember('NDP()', M.var_names_full);
M.var_init(NDP_model_idx) = data.y( strcmp(data.names, 'NDP') ,1); % 

% Get ids and names of non-observed variables to display from the <VAR> section of *.mdf file
% Omitting composite observables - since also want to plot them individually
% Find matches, convert to logical index
obs_vars = M.obs( cellfun('isempty',strfind(M.obs,'+')) );
% Get matches, remove x, convert to double
obs_var_ids = cellfun(@(x) str2num(x(2:end)), obs_vars);
nonObsVars = 1:M.nvar;
nonObsVars(obs_var_ids) = [];
non_obs_names = M.var_names_full(nonObsVars);
n_non_obs_vars = numel(non_obs_names);

% Add names for vars: Observed - from data.y. Non-observed - from *.mdf index.
names = {obs_names{:}, non_obs_names{:}};

% Indexes of parameters to be optimized during model fitting:
% in v2 - added two more params [9 10] - for kmgon and kmgoff
% in v3 - added 11 12 - kT7on kT7off kdegT7

act = 1:M.npar;
act(non_opt_param) = []; % If want to drop some params from optimization.

end; % f_load_data_and_model



%% Find parameters -- this is where model fit happens.
%=============================
%=============================
options = optimset('TolX',1e-6);

% Find the parameters of the data
ss = fitModel_v2c(M,data,act,[],[],[],options);    


%% Below stuff not really critical for bootstrap
%=================================================
if f_plot_solution
    
n_fits = numel(ss);

subplots = M.nobs + n_non_obs_vars;

% Number of rows to display
rows = 2;
    
FIG=FIG+1;
figure(FIG);
if subplots < 10, fig_shift = 0; else fig_shift = -200; end;
set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [fig_shift 0 subplots*sbWidth rows*sbSide]) %#ok<RPMT1> % [x y width height]

if flag_common_yscale
    max_y = max(max(data.y,[],2));
end

% Iterate over observables
%=======================
for i=1:subplots
    subplot(rows,subplots,i);
    hold on;
    
    %%% Plot solutions (fits from final ODE solution)
    %=======================
    for bx=1:n_fits
        if i <= M.nobs        
            plot(ss(bx).int0.t,ss(bx).int0.y(i,:),'-','linewidth',1,'Color','black');
        else % plot non-observed variables
            plot(ss(bx).int0.tr,ss(bx).int0.xr(nonObsVars(i-M.nobs),:),'-')
            % Also set the axis limits for non-obs, since did n do this when
            % plotting experimental data
            ylim([0 max(ss(bx).int0.xr(nonObsVars(i-M.nobs),:))*1.2]);
            xlim([0 max(tspan)]);
        end    
    end; clear bx

    %%% Plot experimental data
    %=======================
    if i <= M.nobs        
        plot(data.t,data.y(i,:),'o','Color',colors(i,:));
        errorbar(data.t,data.y(i,:),data.s(i,:),'Color',colors(i,:));
        if flag_common_yscale
            ylim([0 max_y*1.2]);
        else            
            ylim([0 max(data.y(i,:)*1.2)]); % TMP TMP
        end        
    end
    xlim([0 max(tspan)]);
    title(names(i),'FontSize',fSize);    
%     xlabel('time (min)','FontSize',11);
    if any(i == 1:M.nobs)
        ylabel('mM','FontSize',fSize);
    end
        
    hold off;
end; clear i

end; % f_plot_solution


%=================================================
s = ss; % rename to match the old naming used below - for residuals and stats
%=================================================


%% Plot residuals
%=================================================
if f_plot_residuals % && ~optns.useBootstrap
totalNRMSD = 0;

% Convert residuals from vector into matrix of correct form.
% Takes the number of observables, and splits the vector accordingly.    
res = reshape((s.resManual)', M.nobs, length(data.y)); 

for i=1:M.nobs
    subplot(2,subplots,subplots+i)
    data_res = res(i,:);    
    plot(data.t, data_res,'o','Color',colors(i,:))
    ylim([-max(abs(data_res)) max(abs(data_res))]);
    xlim([0 max(tspan)]);
    rmsd = sqrt( sum(data_res.^2) / length(data_res) ); % sum of squared errors, divided by number of errors, square rooted

    % for normalization previously used int3.y(i,:) solution.
    % now taking this solution from insides of s/ss: s.int0.y(i,:)
    nrmsd = rmsd / max(s.int0.y(i,:)); % RMSD normalized by max variable range

    if ~f_show_NRMSD_etc_for_each_plot
        %%% Only NRMSD - from aSMN_v5m_IN36_IN37_v2b.m:
        title_string = ['NRMSD = ' num2str(nrmsd*100,2) '%'];
    else    
        %%% NRMSD with other stuff - from aSMN_v5m_IN36_IN37_v2b7_SD_CI_from_MSE_for_Parpan.m
        title_string = {...
            ['NRMSD = ' num2str(nrmsd*100,2) '%'],...
            ['Chi2red = ' num2str(s.i_stat.chi2red(i), '%.3g' )],...
            ['GOF = ' num2str(s.i_stat.gofP(i), '%.3g' )],...
            };
    end
    title( title_string ,'FontSize',fSize);
    
    totalNRMSD = totalNRMSD + nrmsd;
end; clear i

end % plot residuals


%% Display constants and statistics
%=================================================
if f_disp_constants_n_stats %% && ~optns.useBootstrap

fprintf(1,'\n= On CI/SD errors see Supplementary Protocol or '); 
fprintf(1,'\n XXUncertainty tag in "a_fit.m" function.\n\n');

%%% XXUncertainty
% fprintf(1,'\n\n== The below parameter uncertainties (CI/SD), based on'); 
% fprintf(1,'\n model fits appear to provide only rough approximations in our '); 
% fprintf(1,'\n setup - mostly because (1) we dont have proper error estimates on'); 
% fprintf(1,'\n data points, and (2) optimization is done not in the log-space!'); 
% fprintf(1,'\n Bootstrap analysis was also providing very narrow error ranges. '); 
% fprintf(1,'\n So the the best way to assess parameter(constants) uncertainty'); 
% fprintf(1,'\n in these IVTNMR datasets at the moment appears to be by using replicate '); 
% fprintf(1,'\n measurements (ideally with different DNA template/etc batches).');

% Get names and values of the final - ONLY OPTIMIZED PARAMS (since errors
% are only calculated for those!).
p = struct(); % create structure to hold parameters
opt_par_names = M.par_names(act);
opt_pars = s.popt(act);
n_opt_pars = numel(act);
% below s.sem / s.sd / s.conf - already have the size of ACT

results = cell(1, n_non_obs_vars); % same size as n of graphs w/o residuals
n_results_panels = n_non_obs_vars;

%%% Collect stats - ONLY FOR OPTIMIZED PARAMS!
for i=1:n_opt_pars
    % Indexes of result display go from the end: result(3), result(2), ..
    % Now automatically placing below results closest to residuals graphs.
    % Do this only if have enough space (n_results_panels>2)
    % If have more space (n_non_obs_vars > 2) and many params (n_opt>6)- split in two groups.
    % // if n_non_obs_vars == 1 - this info is not going to "results" at all.
    if n_results_panels > 1
        if n_results_panels == 2
            location_of_res = n_results_panels; % all in one place - 2nd from end (leftmost)
        elseif n_results_panels > 2
            if n_opt_pars < 6
                location_of_res = n_results_panels; % all in one place - leftmost
            elseif i <= (n_opt_pars/2)
                location_of_res = n_results_panels; % first half - leftmost
            elseif i > (n_opt_pars/2)
                location_of_res = n_results_panels-1; % second half - next from leftmost
            end                
        end 
        results(location_of_res) = strcat(results(location_of_res), opt_par_names{i}, {' = '}, sprintf('%0.2g', opt_pars(i)), {char(10)} );
    end

    % Get parameters
    p.( opt_par_names{i} ) = opt_pars(i);
    % Get errors
    ci.( opt_par_names{i} ) = s.conf(i);
end; clear i

% Another structure - with ALL params (not only the optimized ones)
pall = struct(); 
for i=1:M.npar
    pall.( M.par_names{i} ) = s.popt(i);
end; clear i

% % Do not try to calc KMgNTP if there are no related constants in the model
% if ~any(ismember({'kmgoff','kmgon'},opt_par_names))
%     KMgNTP = nan;
%     CI_KMgNTP = nan;
% else
%     KMgNTP = p.kmgoff/p.kmgon; % [mM]
%     CI_KMgNTP = sqrt( (ci.kmgoff/p.kmgoff)^2 + (ci.kmgon/p.kmgon)^2 ) * (KMgNTP);
% end

% Calc UP1-SMN Kd only if have constants in the model
if ~any(ismember({'kon','koff'},opt_par_names))
    KdPR = nan;
    CI_KdPR = nan;
else
    [KdPR, CI_KdPR] = calc_KdPR(p.kon, ci.kon, p.koff, ci.koff, 1); % do not normalize by RNA length
end

if ~ismember('kabortive',opt_par_names)
    [kcatnorm,CI_kcat] = calc_kcat(p.kcat, ci.kcat);
    kcatnorm_name = 'kcat norm';
else
    [kcatnorm,CI_kcat] = calc_kcat(p.kcat, ci.kcat, p.kabortive, ci.kabortive);
    kcatnorm_name = 'kcka norm'; 
end

[KeqS, CI_KeqS] = calc_KeqS(p.kprecip, ci.kprecip, p.kdissolve, ci.kdissolve);

%%% Indexes of result display go from the end: result(3), result(2), ..
%%% In the "results" structure - ALWAYS displaying summary at last position (results(1).
%%% If there is enough space - all parameters are also displayed -- 
%%% via results(n_results_panels-x) generated above.
results(1) = strcat(data_name, {char(10)},...
    {'Tot.NRMSD | Chi2r | GOF | R^2 '}, {char(10)},...
    {[num2str(totalNRMSD*100,3) '% | ' num2str(s.chi2red, '%.3g') ' | ' num2str(s.chi2.p, '%.3g' ) ' | ' num2str(s.R2,3)]}, {char(10)},...
    { sprintf('%s = ',kcatnorm_name)}, num2str(kcatnorm,2), '(',num2str(CI_kcat,2),')', {' s-1'}, {char(10)},...
    {'KeqS = '}, num2str(KeqS,2), '(',num2str(CI_KeqS,2),')', {' mM'}, {char(10)},...
    {sprintf('KdPR = %.2f(%.1f) uM', KdPR, CI_KdPR)}, {char(10)},...
    {'tlaps = '}, num2str(toc,'%.0f'), {' s'},...
    {});
    %     {'kc+ka norm = '}, num2str(kcka_norm,2), '(',num2str(CI_kcka,2),')', {' s-1'}, {char(10)},...
    %     {'KdPR = '}, num2str(KdPR,2),'(',num2str(CI_KdPR,2),')',  {' uM'}, {char(10)},...
    %     {'KdNTP = '}, num2str(KdNTP,2),'(',num2str(CI_KdNTP,2),')',  {' uM'}, {char(10)},...
    %     {'\DeltaG = '}, num2str(dG,2), {' kcal/mol'}, {char(10)},...
    %     {'KMgNTP = '}, num2str(KMgNTP,2),'(',num2str(CI_KMgNTP,2),')',  {' mM'}, {char(10)},...
    %     {'KdAbrt = '}, num2str(KdAbrt,2),'(',num2str(CI_KdAbrt,2),')',  {' uM'}, {char(10)},...

errLegend = '-';
for i=0:(length(results)-1)
    subplot(rows,subplots,subplots*2-i); plot(0,0); axis off; % create a dummy plot
    text(-1,1,results(i+1),'FontSize',fSize,'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'top', 'Interpreter', 'None');
    if i==(length(results)-2) % on the second plot from left side
        title( errLegend ,'FontSize',fSize); % Add notation for what is displayed as error
    end
end; clear i;

end

%% Export figure
%===========================
if ~optns.useBootstrap
    fig_folder = fullfile(modelSubFolder,'figure_images');
    fig_name = sprintf('%s_%s_%s',mfilename,data_name(1:6),datestr(now,'yymmdd_HHMM'));
    fig_path = fullfile(fig_folder,fig_name);
    save_figure(FIG,fig_path);        
end;

%% Summary of results - returned when "a_fit" is executed externally
%===========================
if ~optns.useBootstrap
    fit_data.name = data_name;
    fit_data.hn_peak = optns.HN_peak_for_fit;
    
    if ~strcmp(optns.HN_peak_for_fit,'all')
        fit_data.max_csp = maxcsp_list(optns.HN_peak_for_fit);
    else
        fit_data.max_csp = 'multi-peak';
    end
    
    fit_data.maxcsp_scaling = optns.maxcsp_scaling;    
    fit_data.kcatnorm = kcatnorm;
    fit_data.CI_kcatnorm = CI_kcat;
    fit_data.KeqS = KeqS;
    fit_data.CI_KeqS = CI_KeqS;
    if exist('KdPR','var')
        fit_data.KdPR = KdPR;
        fit_data.CI_KdPR = CI_KdPR;
    end;
    fit_data.kdephosNTP = pall.kdegNTP;
    fit_data.data = data;
    fit_data.ss = ss;

    toc;

    pall
else
    fit_data = 'not exported in bootstrap';
end;

end % whole function

%===========================%===========================%===========================
%===========================%===========================%===========================
%===========================%===========================%===========================

% Auxiliary functions to calculate the final constants from the model
% values.

function [kcatnorm, CI_kcat] = calc_kcat(kcat, kcat_ci, kabortive, kabortive_ci)
    fprintf(1,'kcatnorm calculated within independent function. \n');
    global a_fit_T7DNAcomplex_concentration;
    kcat_norm_factor = a_fit_T7DNAcomplex_concentration * 60; % mM T7-DNA cplx * 60 seconds in minute
    if nargin==2
        kcatnorm = kcat ./ kcat_norm_factor;
        CI_kcat = kcat_ci ./ kcat_norm_factor;        
    elseif nargin==4
        kcatnorm = (kcat+kabortive) ./ kcat_norm_factor;
        CI_kcat = (kcat_ci+kabortive_ci) ./ kcat_norm_factor;        
    else
        error('wrong number of arguments in calc_kcat');
    end
end

function [k, ci] = calc_KeqS(kprec, kprec_ci, kdis, kdis_ci)
    KeqS = kprec ./ kdis; % [mM !!] KeqSolubility - was Ksp
    global a_fit_Mg_total; Mg_total = a_fit_Mg_total;
    global a_fit_NTP_initial_conc; NTP_conc = a_fit_NTP_initial_conc;

    % Fraction of Mg in NTP complex.
    % from Kd=0.3mM (literature), Mg=24mM, NTP=20mM - gives 0.7871.
    PL_fraction_of_P = @(Kd,L,P) ((Kd + L + P) - sqrt((Kd+L+P)^2-(4*P*L))) / (2*P);
    NTPMg_bound_fraction = PL_fraction_of_P(0.3, NTP_conc, Mg_total);
    
    % We approximate that the fraction of NTP/RNA bound Mg is ~constant 
    % - i.e. when NTP is converted to RNA - Mg remains bound.
    Mg_free = (1 - NTPMg_bound_fraction)*Mg_total; 
    k = KeqS ./ Mg_free; % See AnchorKsp for derivation: Keq = Keq' / Mg_free
    ci = sqrt( (kprec_ci./kprec).^2 + (kdis_ci./kdis).^2 ) .* (kprec./kdis) ./ Mg_free;
end

function [KdPR, CI_KdPR] = calc_KdPR(kon, kon_ci, koff, koff_ci, rna_length)
    KdPR = koff./kon ./ rna_length .* 1000; % [uM]    
    CI_KdPR = sqrt( (koff_ci./koff).^2 + (kon_ci./kon).^2 ) .* (KdPR);
end


function save_figure(fig_handle,file_path)
%% Save figure as pdf
%=============================================================
    set(fig_handle,'Units','Inches');
    pos = get(fig_handle,'Position');
    set(fig_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(fig_handle,file_path,'-dpdf','-r1200');
%     print(fig_handle,fullfile('figure_images',file_name),'-dpdf','-r1200')
%     print(fig_handle,fullfile(file_name),'-dpdf','-r1200')

% CONVENIENT TO USE with mfilename variable - which returns the name of the current file!
% 
% -r option sets the resolution (also for the VECTOR IMAGES!) -- set to 300-1200 if contour approximation is bad.
% -r0
% -r300
% -r1200
end
