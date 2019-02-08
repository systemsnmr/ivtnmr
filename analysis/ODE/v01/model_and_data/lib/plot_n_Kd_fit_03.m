function plot_n_Kd_fit_03(nmr_data_path, dset_names, y_observable, select_peaks, scale_rna_aborts_fract, fix_dsat, scale_rna0_too, maxcsp_scaling)
% Imports, plots, fits Kd in IVTNMR datasets.

global FIG; 
global SB; 

% v01 - based on plot_RNA_vs_CSP.m
%     - see Evn 'Kd fits (for ...)' evernote:///view/350340/s4/221a920b-8b0d-498e-b225-65ec33475ffd/221a920b-8b0d-498e-b225-65ec33475ffd/
%     + Disp dataset short name outside
%     + Calc common max_Y (only for peaks existing in all sets)
%     + Call Kd fit function
%     + Display Kd fit results on each graph
%     + Display Kd fit results as a table
%     + Save all Kds locally as kd_arr dataset
% v02 - adapted for external running (from runall.m workflow)
% v03 - allows usage of fixed dsat

% TODO:
%     - Implement OPTIONS structure

% Flags to run different parts of code
run_kd_fits = 1;

save_kd_to_local_ivtnmr = 1;
save_kd_to_global_ivtnmr = 1;
reimport_ivtnmr_from_local = 1;
% overwrite_kd_struct_in_global_ivtnmr = 1; % not used at the moment

plot_data_n_etc = 1;
    overlay_kd_fits = 1; % nested in the plot_data_n_etc
    disp_kd_on_plot = 1; % nested in the plot_data_n_etc

make_kd_summary_dataset = 1;
disp_kd_as_table = 1;

addpath('./sublib');

%% Params to run "locally"
%===========================
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    nmr_data_path = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra';
    dset_names = {...
        '170308_IN60a_pR02_co-A1R1_303K_600'
        '170310_IN70a_pR02_co-UP1_303K_600'
        '170503_IN61a_SMN1_co-A1R1_303K_600'
        '170313_IN71a_p114_co-UP1_303K_600'
        '170309_IN62a_p214_co-A1R1_303K_600'
        '170426_IN72a_SMN2_co-UP1_303K_600'
        '170607_IN65b_EV2_co-A1R1_303K_600'
        '170515_IN75b_EV2_co-UP1_303K_600'
        };
    y_observable = 'RNA_from_NTP';
    
    select_peaks = {'21' '22' '33' '56' '75'};
    
    scale_rna_aborts_fract = 0.3; % set empty or zero if don't want to use
    
    fix_dsat = [];
    fix_dsat = 'global_maxcsp_for_each_peak';
    
    scale_rna0_too = 1;
else
    disp('');
end


global DATA_IVTNMR;

if exist('./lib','dir')
    datasave_folder = strcat(DATA_IVTNMR, '/'); % if running in root (runall)
else
    datasave_folder = strcat('../', DATA_IVTNMR, '/'); % if running from lib directory
end


global DEBUG; DEBUG = 0;

%% Actual script
%===========================
n_sets = numel(dset_names);
d = cell(n_sets,1);

name_split = cellfun(@(x) regexp(x, '_', 'split'), dset_names, 'unif', 0);
dset_snames = cellfun(@(x) x{2}, name_split, 'unif', 0); % extract short name

for i=1:n_sets
    load(fullfile(nmr_data_path,dset_names{i},'ivtnmr.mat'));
    d{i} = s;
end


%% Parse dsat setting
%=====================
%%% Parse dsat value: empty, auto_%d, getmax_for_same_peak, CONTAINER with peak names == maxcsp
dsatIsFixed = exist('fix_dsat','var') && ~isempty(fix_dsat);

if dsatIsFixed
    
    if ischar(fix_dsat) && strcmp(fix_dsat(1:4), 'auto')
        fprintf(1,'Leaving fix_dsat = %s\n', fix_dsat);

    elseif ischar(fix_dsat) && strcmp(fix_dsat, 'global_maxcsp_for_each_peak')
        fprintf(1,'fix_dsat=%s - Getting common max_csp for each peak\n', fix_dsat);
        hn_names = unique( ivtnmr_get_hn_peak_names(nmr_data_path, dset_names) );
%         hn_names = {'21','22','33','56','75'};
        % overwrites fix_dsat into container map
        fix_dsat = ivtnmr_get_maxcsp_for_peaks(nmr_data_path, dset_names, hn_names);        

    elseif isa(fix_dsat, 'containers.Map')
        fprintf(1,'fix_dsat already a container.Map - using directly to get peak max_csp\n');

    else
        error('Something wrong with fix_dsat');
    end
    
end % dsatIsFixed




%% Get all Kds
%==================
if run_kd_fits
for iSet=1:n_sets
    s = d{iSet};
    npk = s.hn.n_peaks;
    obs_idx = strcmp(s.p31.obs_names,y_observable);
    end_point = min(numel(s.hn.time), numel(s.p31.time));
    
%     s.hn.kd(npk,1) = struct(); % NOT PREALLOCATING - CUZ NEED TO UPDATE
%     FIELDS OTHERWISE!
    
    for iPeak=1:npk
        pkname = s.hn.names{iPeak};
        
        if DEBUG
            fprintf(1,'Fitting Kd for iSet=%i (%s), iPeak=%i (%s) ...\n', iSet, s.sname, iPeak, s.hn.names{iPeak});
        end
        % just example of what data is used in below plotting
%         plot(s.p31.obs_hn_time(obs_idx,:), s.hn.csp(p,:));        
        ligand_conc = s.p31.obs_hn_time(obs_idx,1:end_point);

        % Scale
        rna_fraction_is_defined = exist('scale_rna_aborts_fract','var') && ~isempty(scale_rna_aborts_fract) && (scale_rna_aborts_fract ~= 0);        
        % THIS can and needs to be optimized! TODO
        if rna_fraction_is_defined
            if ~strcmp(s.rna_name, 'RNA0') || scale_rna0_too
                disp('=======================');
                disp('=== RNA conc scaled ===');
                disp('=======================');
                ligand_conc = ligand_conc .* scale_rna_aborts_fract;
            end
        end
        
        target_csp = s.hn.csp(iPeak,1:end_point);
        target_conc = 0.15; % mM !!!
        plot_results = 0;        
        
        fprintf(1,'iPeak=%i (name = %s) \n', iPeak, pkname);
%         if overwrite_kd_struct_in_global_ivtnmr
           %%% IF RE-ACTIVATING THIS AGAIN: need to change such that not
           %%% WHOLE s.hn.kd sub-structure is blanked, but only specific
           %%% element of it: s.hn.kd(iPeak) -- otherwise each new peak
           %%% kills all previous info.
%             s.hn = rmfield(s.hn, 'kd');
            fprintf(1,'\n================\n iSet / iPeak = %i / %i \n================\n', iSet, iPeak);
                       
            if dsatIsFixed
                
                if isa(fix_dsat,'containers.Map')
                    s.hn.kd(iPeak) = fit_Kd_03(ligand_conc, target_csp, target_conc, plot_results, fix_dsat(pkname).*maxcsp_scaling );
                else
                    % Either 'auto' or some other value parsed downstream:
                    s.hn.kd(iPeak) = fit_Kd_03(ligand_conc, target_csp, target_conc, plot_results, fix_dsat);
                end
                
            else % Can call a different function here if dsat not fixed.
%                 s.hn.kd(iPeak) = fit_Kd_03(ligand_conc, target_csp, target_conc, plot_results, []);
                s.hn.kd(iPeak) = fit_Kd_03(ligand_conc, target_csp, target_conc, plot_results);
                
            end
%         else
%             error('Kd structure already exists. Modify overwrite_kd_struct_.. flag to overwrite.');
%         end
    end
    fprintf(1,'== Finished Kd fits for %s\n\n', s.sname);
    
    %%% Save
    if save_kd_to_local_ivtnmr
        save( strcat(datasave_folder, 'ivtnmr_', d{iSet}.sname,'.mat'), 's');
    end
    if save_kd_to_global_ivtnmr
        save( fullfile(nmr_data_path, d{iSet}.name, 'ivtnmr.mat') ,'s');
    end    
end
end % run_kd_fits

%% Reimport data with kd
%==========================
if reimport_ivtnmr_from_local
clear d;
d = cell(n_sets,1);
for iSet=1:n_sets
    load( strcat(datasave_folder, 'ivtnmr_', dset_snames{iSet},'.mat') );   
    d{iSet} = s;
end    
end


%% Select some peaks
%===================================
if exist('select_peaks','var') && ~isempty(select_peaks)
    for i=1:n_sets
        idx = ismember(d{i}.hn.names, select_peaks);
        d{i}.hn.names = d{i}.hn.names(idx);
        d{i}.hn.H = d{i}.hn.H(idx,:);
        d{i}.hn.N = d{i}.hn.N(idx,:);
        d{i}.hn.csp = d{i}.hn.csp(idx,:);
        d{i}.hn.kd = d{i}.hn.kd(idx);
        d{i}.hn.n_peaks = numel(select_peaks);
    end
end


%% Check num of peaks in diff sets
%===================================
n_peaks_array = cellfun(@(x) x.hn.n_peaks, d);
max_peaks = max( n_peaks_array );
min_peaks = min( n_peaks_array );

% Calc common max y-scale for all sets with same (min) number of peaks
max_csp = nan(min_peaks,1);
for iPeak=1:min_peaks
    max_csp(iPeak) = max( cell2mat( cellfun(@(x) max(x.hn.csp(iPeak,:)), d, 'unif', 0) ) );
end



%% Plot
%===================
if plot_data_n_etc
rows = n_sets;
columns = max_peaks;

sbSide = 150;
fSize = 11;
SB=0; % global defined at the top of the function
sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

% global defined at the top of the function
if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
FIG=FIG+1;
figure(FIG); set(figure(FIG), 'Color', repmat(0.97,1,3), 'Position', [0 0 columns*sbSide rows*sbSide]);

for i=1:n_sets
    % Set position at beginning of the row
    SB=(i-1)*columns;
    s = d{i};
    
    for p=1:d{i}.hn.n_peaks
        obs_idx = strcmp(d{i}.p31.obs_names,y_observable);
        end_point = min(numel(s.hn.time), numel(s.p31.time));

        SB=SB+1; sb(SB);
%         plot(s.p31.obs(obs_idx,1:end_point), s.hn.csp(p,1:end_point));
        % Still need endpoint, cuz in some cases have NaN at the end (even
        % after interpolation).
        plot(s.p31.obs_hn_time(obs_idx,1:end_point), s.hn.csp(p,1:end_point), 'o');
%         plot(s.p31.obs_hn_time(obs_idx,:), s.hn.csp(p,:));

        if overlay_kd_fits
            hold on;
            if DEBUG
                fprintf(1,'Plotting fitted Kd - set=%i (%s), peak=%i (%s) ...\n', i, s.sname, p, s.hn.names{p});
            end;
            x = s.p31.obs_hn_time(obs_idx,1:end_point);
            y = s.hn.kd(p).csp_fit;
            
            fprintf(1,'\n================\n set / peak = %i / %i \n================\n', i, p);
            plot(x, y, '--r','LineWidth', 2);
        end

        if p > 1
            title(s.hn.names(p));
        else
            title( sprintf('%s (%s)  %s', s.rna_name, s.sname, s.hn.names{p}) );
        end
        axis tight
        if p <= min_peaks % Set common scale for peaks existing in all sets
            ylim([0 max_csp(p)]);
        end
        
        % Display Kd info
        %=====================
        if dsatIsFixed
            kd_val = s.hn.kd(p).popt(1);
            kd_sd = s.hn.kd(p).sd(1);
            dsat_string = sprintf('ds=%.2f(fix)', fix_dsat(s.hn.names{p}) .* maxcsp_scaling);
        else
            kd_val = s.hn.kd(p).popt(2);
            kd_sd = s.hn.kd(p).sd(2);
            dsat_val = s.hn.kd(p).popt(1);
            dsat_sd = s.hn.kd(p).sd(1);
            sprintf('sat=%.2f (%.2f)', dsat_val, dsat_sd);
        end
        
        if disp_kd_on_plot
            fit_summary = {...
                sprintf('kd=%.2f (%.2f) mM', kd_val, kd_sd)...
                dsat_string...
            };
            
            x_lim = get(gca,'xlim');
            y_lim = get(gca,'ylim');
            text(x_lim(1),y_lim(2),...
              fit_summary,...
              'HorizontalAlignment','left',...
              'VerticalAlignment', 'top',...
              'FontSize',fSize);

        end;
                        
    end
    
end
end % plot_data_n_etc


%% Generate dataset with kds
%============================
% Reinitializing all variables - will likely spun-off this to external func
if make_kd_summary_dataset
%     rna_name, rna_length, prot_name, cotr, kd, kd_err, dsat, dsat_err
    n_sets = numel(d);
    n_peaks_array = cellfun(@(x) x.hn.n_peaks, d);
    total_peaks = sum(n_peaks_array);
    max_peaks = max( n_peaks_array );
    min_peaks = min( n_peaks_array );

    % TODO - was trying vectorization to avoid looping and indexes!
%     z = arrayfun(@(x) repmat('a',x,1), n_peaks_array, 'unif', 0)
        % But this did not work:
%     z = arrayfun(@(x,y) repmat(d.sname, y, 1), d, n_peaks_array, 'unif', 0)
    
    % Preallocate
    kd_arr = dataset(...
        {[1:total_peaks]','id'},...
        {cell(total_peaks,1),'sname'},...
        {cell(total_peaks,1),'rna'},...
        {nan(total_peaks,1),'rna_length'},...
        {cell(total_peaks,1),'prot'},...
        {cell(total_peaks,1),'cotr'},...
        {cell(total_peaks,1),'peak'},...
        {nan(total_peaks,1),'kd_mM'},...
        {nan(total_peaks,1),'dsat'},...
        {nan(total_peaks,1),'kd_sd'},...
        {nan(total_peaks,1),'kd_sem'},...
        {nan(total_peaks,1),'dsat_sd'},...
        {nan(total_peaks,1),'dsat_sem'}...
    );    
    
    for i=1:n_sets       
        start_id = sum( n_peaks_array(1:i-1) )+1;
        end_id = start_id+n_peaks_array(i)-1;
%         fprintf(1,'%i - %i\n', start_id, end_id);

        kd_arr.sname(start_id:end_id) = {d{i}.sname};
        kd_arr.rna(start_id:end_id) = {d{i}.rna_name};
        kd_arr.rna_length(start_id:end_id) = d{i}.rna_length;
        kd_arr.prot(start_id:end_id) = {d{i}.prot};
        kd_arr.cotr(start_id:end_id) = {d{i}.cotr};
        kd_arr.peak(start_id:end_id) = d{i}.hn.names;
        
%         arrayfun(@(x) x.popt(2), d{i}.hn.kd(1:end)) % test visualize KDs
%         arrayfun(@(x) x.popt(1), d{i}.hn.kd(1:end)) % test visualize

        % Set Kd values
        if dsatIsFixed
            kd_popt = 1;
        else
            kd_popt = 2;
        end
        
        kd_arr.kd_mM(start_id:end_id) = real( arrayfun(@(x) x.popt(kd_popt), d{i}.hn.kd(1:end)) ); % force REAL (drop IMAGINARY)        
        kd_arr.kd_sd(start_id:end_id) = real( arrayfun(@(x) x.sd(kd_popt), d{i}.hn.kd(1:end)) ); % force REAL (drop IMAGINARY)
        kd_arr.kd_sem(start_id:end_id) = real( arrayfun(@(x) x.sem(kd_popt), d{i}.hn.kd(1:end)) ); % force REAL (drop IMAGINARY)
                
        if ~dsatIsFixed % if dsat not fixed - also save dsat            
            kd_arr.dsat(start_id:end_id) = real( arrayfun(@(x) x.popt(1), d{i}.hn.kd(1:end)) ); % force REAL (drop IMAGINARY)
            kd_arr.dsat_sd(start_id:end_id) = real( arrayfun(@(x) x.sd(1), d{i}.hn.kd(1:end)) ); % force REAL (drop IMAGINARY)
            kd_arr.dsat_sem(start_id:end_id) = real( arrayfun(@(x) x.sem(1), d{i}.hn.kd(1:end)) ); % force REAL (drop IMAGINARY)
        end        

        for p=1:n_peaks_array(i)            
            if disp_kd_as_table
                if i==1 && p==1
                    fprintf(1,'set \t rna \t len. \t prot \t cotr \t peak \t Kd \t Kd_sd\n');
                end                
                    if dsatIsFixed
                        fprintf(1,'%s \t %s \t %i \t %s \t %s \t %s \t %d \t %d\n',...
                        d{i}.sname, d{i}.rna_name, d{i}.rna_length, d{i}.prot, d{i}.cotr, d{i}.hn.names{p}, d{i}.hn.kd(p).popt(1), d{i}.hn.kd(p).sd(1));
                    else
                        fprintf(1,'%s \t %s \t %i \t %s \t %s \t %s \t %d \t %d\n',...
                        d{i}.sname, d{i}.rna_name, d{i}.rna_length, d{i}.prot, d{i}.cotr, d{i}.hn.names{p}, d{i}.hn.kd(p).popt(2), d{i}.hn.kd(p).sd(2));
                    end
            end
        end
        
        
%         {cellfun(@(x) x.sname, d, 'unif', 0),'sname'}...
    end
    
    save(strcat(datasave_folder, 'kd_arr.mat'),'kd_arr');
    
end % make_kd_summary_dataset
    
end

