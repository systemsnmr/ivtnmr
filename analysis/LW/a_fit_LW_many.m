function a_fit_LW_many(nmr_data_path, dsets, imino_ppm_list, optns)
% a_fit_LW_list: largely same as a_fit_LW, but does many iminos at once.
% (import 1D1H spectra (expnos 5xxx); fit imino linediwths; plot
% time-resolved imino linewidths. Saves final result into EPS.
%
% Run with e.g.:
%   a_fit_LW('/opt/topspin/.../...', {'dataset1', 'dataset2'}, [14.2, 13.1], optns)
%
% Additional features:
% - has an option to visualize the LineWidth fits - to check the quality
% of the fitting.
%
% For details about setup and analysis protocol see github.com/systemsnmr/ivtnmr

%% Settings
%=================
if nargin == 0
    clear; close all;
    
    analysis_root = fullfile(pwd, '..');
    nmr_data_path = fullfile(analysis_root, '..', 'spectra');
%     nmr_data_path = '/Volumes/Data/yar/_eth2/data_NMR/spectra/';
%     nmr_data_path = '/Volumes/Data/yar/scratch/SystemsNMR_data_2019';
                
    dsets = {...
%         '180314_IN72b_SMN214_co-NUP1_303K_600'
%         '170322_IN75a_EV2_co-UP1_303K_600'
        '170322_IN75aa_EV2_co';
        '170402_IN75jj_EV2_post';
    };

    optns.fit_width = [0.3, 0.1, 0.07, 0.12]; % in ppm. (earlier was using half-width here and not /2 below).
    imino_ppm_list = [14.033, 13.518, 13.22, 11.937]; % 13.179
    % to select one particular imino
%     imino_ppm_list = imino_ppm_list(3);
%     optns.fit_width = optns.fit_width(3);

    optns.reread_or_load_spectra = [1 1 1]; % Three parameters: [reread, save, load]
    optns.refit_or_load_lw = [1 1 1]; % Three parameters: [refit, save, load]    
            
    optns.adjust_fitted_range_to_center_on_peak_max = 1; % not finished yet
    
    %%% Display: LW fits
    optns.displayFits = 1; % Shows fitting of spectra
    optns.ppm_window_to_display = 0.5; % ppm around the peak - NOT IMPLEMENTED YET
    optns.n_spectra_to_disp = 5; % Spectra shown when showing LW fits.    

    %%% Display: time-resolved LW
    optns.xaxis_start = 8; % h    
        
    %%% Less frequently changed    
    optns.temperature=303;
    optns.save_figure = 1;
    optns.display_lw = 1; % LW versus time
    optns.plotSym = 'o-';
        
    % takes the second underscore-separated element as the ID of expt.
    dset_split_arr = cellfun(@(x) regexp(x, '_', 'split'), dsets, 'un', 0);    
    dset_id = cellfun(@(x) x{2}, dset_split_arr, 'un', 0);
    dset_id_long = cellfun(@(x) sprintf('%s_%s_%s', x{2}, x{3}, x{4}), dset_split_arr, 'un', 0);
    clear dset_split_arr;
    optns.legend = dset_id_long;
    
%     optns.save_to_global_ivtnmr
end;

scratchfile_prefix = strcat(mfilename);
scratch_dir = 'datasave'; % for intermediate data
datasave_dir = 'datasave'; % for final data

%%% Defaults for optional parameters (are set all together at the moment)
% if ~isfield(optns, 'temperature') || isempty(optns.temperature); optns.temperature=303; end;

    
%% Extra options (usually unused if running by external call - thats why not)
%===========================
flag_unittest_spectra = 0;
flag_unittest_LW = 0;

flag_display_fit_details = 0;


%% "Fixed" Settings
%=====================
fSize = 9;
scale = 1;
sbSide = 220;
sbWidth = sbSide * scale;
sbHeight = sbSide * scale;
global FIG; FIG=0;


%% Initial prep / preallocation
%=====================
spec_flags = num2cell(optns.reread_or_load_spectra);
lwfit_flags = num2cell(optns.refit_or_load_lw);
[flag_reread_spectra, flag_save_spectra_to_scratch, flag_load_spectra_from_scratch] = ...
    deal(spec_flags{:});
[flag_refit_lw, flag_save_to_scratch, flag_load_from_scratch] = ...
    deal(lwfit_flags{:});

addpath( genpath(fullfile(pwd,'lib')) );
dset_full_paths = cellfun(@(x) fullfile(nmr_data_path,x), dsets, 'un', 0);
n_sets = numel(dsets);
n_iminos = numel(imino_ppm_list);

if numel(optns.fit_width)==1; optns.fit_width=repmat(optns.fit_width, n_iminos, 1); end;

colors = lines(n_sets);
plotLW_list = ones(1, n_sets);

%%% Preallocate to store variables for several sets:
% (use cells, cuz different dsets may h diff amount of data
% lw = cell(n_sets,1);
lw = cell(n_sets,n_iminos);
fitted_ppm_center = cell(n_sets,n_iminos);

% to store spectra params:
expt_names = cell(n_sets,1);
time = cell(n_sets,1);
S = cell(n_sets,1);

% Create dirs if don't exist yet.
if ~exist( datasave_dir ,'dir'), mkdir(datasave_dir), end;
if ~exist( scratch_dir ,'dir'), mkdir(scratch_dir), end;

%% Actual analysis
%=====================

%=====================
if flag_reread_spectra
for ds=1:n_sets    
        files = cell2mat(arrayfun(@(x) str2num(x.name), dir(fullfile(nmr_data_path, dsets{ds},'6*')), 'un', 0));
        expnos = files(files>=6000 & files<=6999)';

        time0 = get_time0(nmr_data_path, dsets{ds});
        time{ds} = getNMRTime_02(dset_full_paths{ds},expnos,time0); % time = getNMRTime(datasetName,experiments,time0)

        % TODO:
        % expt_names - derive from TIME
        expt_names{ds} = arrayfun(@num2str, expnos, 'unif', 0);

        rbSpectra = cellfun(@(x) fullfile(nmr_data_path, dsets{ds},x,'pdata/1/1r'), ...
            arrayfun(@num2str, expnos, 'un', 0), 'un', 0);
        
        S{ds} = rbnmr(rbSpectra);
        if flag_save_spectra_to_scratch; save( fullfile(scratch_dir, sprintf('%s_sp.mat', scratchfile_prefix)), 'S', 'time', 'expt_names'); end;                                
end;
end; % flag_reload spectra

if flag_load_spectra_from_scratch; load(fullfile(scratch_dir, sprintf('%s_sp.mat', scratchfile_prefix))); end;

if flag_unittest_spectra
    UT_ref = load(fullfile(scratch_dir,'a_fit_LW_sp_UT1.mat'));
    UT_ref = UT_ref.S;
    UT_struct_fields(UT_ref,S);
end;

%=====================
if flag_refit_lw
for ds=1:n_sets            
    
    n_spectra = numel(S{ds});
    
    for iImino=1:n_iminos % a_fit_LW_many
        
        imino_ppm = imino_ppm_list(iImino); % a_fit_LW_many
        fit_half_width = optns.fit_width(iImino)/2;
        
        if optns.adjust_fitted_range_to_center_on_peak_max
            pkBounds_0 = [imino_ppm+fit_half_width imino_ppm-fit_half_width];
            % TODO: spin out below "ints" and "ppms" functions into anonymous / external
            % functions - they are repeating.
            cutIds_0 = cellfun(@(x) find(x.XAxis < pkBounds_0(1) & x.XAxis > pkBounds_0(2)), S{ds}, 'un', 0);
            ints_0 = cellfun(@(source,extract) source.Data(extract), S{ds}, cutIds_0, 'un', 0);
            ppms_0 = cellfun(@(source,extract) source.XAxis(extract), S{ds}, cutIds_0, 'un', 0);
            
            [~,max_pos] = cellfun(@(x) max(x), ints_0, 'un', 0);
            ppm_at_max_pos = cellfun(@(x,y) x(y), ppms_0, max_pos, 'un', 0);
            
            % Below results in individual ppm range for fitting each peak.            
            pkBounds_array = cellfun(@(pos) pos+[fit_half_width, -fit_half_width], ppm_at_max_pos, 'un', 0);
            cutIds = cellfun(@(x,bounds) find(x.XAxis < bounds(1) & x.XAxis > bounds(2)), S{ds}, pkBounds_array, 'un', 0);
            ints = cellfun(@(source,extract) source.Data(extract), S{ds}, cutIds, 'un', 0);
            ppms = cellfun(@(source,extract) source.XAxis(extract), S{ds}, cutIds, 'un', 0);
        else        
            % Cut out the part for fitting:
%             pkBounds = [imino_ppm+fit_half_width imino_ppm-fit_half_width];
%             cutIds = cellfun(@(x) find(x.XAxis < pkBounds(1) & x.XAxis > pkBounds(2)), S{ds}, 'un', 0);
            % Below uses the same fitting range for all peaks.
            pkBounds_array = repmat({[imino_ppm+fit_half_width imino_ppm-fit_half_width]},n_spectra,1);
            cutIds = cellfun(@(x,bounds) find(x.XAxis < bounds(1) & x.XAxis > bounds(2)), S{ds}, pkBounds_array, 'un', 0);
            ints = cellfun(@(source,extract) source.Data(extract), S{ds}, cutIds, 'un', 0);
            ppms = cellfun(@(source,extract) source.XAxis(extract), S{ds}, cutIds, 'un', 0);
        end;


        %% Fit Lorentzian and calc data stats
        %=====================================

        % [YPRIME PARAMS RESNORM RESIDUAL] = LORENTZFIT(X,Y,P0,BOUNDS,NPARAMS,OPTIONS)
        [yprime,params,resnorm,residual] = deal( cell(n_spectra,1) );

        % Pre-allocate for results
        [ymin,ymax,yspread,xATymax,P1,P2,P3,C,FWHM,A,yATxc] = deal( nan(n_spectra,1) );

        % Pre-allocate for post-regression diagnostics
        [R2,R2adj,dof,mse,chi2,chi2red,gofP,gofPred] = deal( nan(n_spectra,1) );
        % % Alternatively - preallocate in one structure to keep fit results:
        % fitRes = repmat(struct('popt',[],'resnorm',[],'residual',[],'exitflag',[],'output',[],'lambda',[],'jacobian',[],'R2',[],'R2adj',[]),...
        %     numOfFits,1);

        man_FWHM = nan; % in TopSpin - at 6.88 intensity - half-way from baseline hump
        P0 = [];

        for i=1:n_spectra
            % Fit
            [yprime{i} params{i} resnorm{i} residual{i}] = lorentzfit2(ppms{i},ints{i},P0,[],'3'); % zero baseline (no constant term)

            % Calc stats    
            ymin(i) = min(ints{i});
            [ymax(i), ymax_idx] = max(ints{i});
            yspread(i) = ymax(i)-ymin(i);
            xATymax(i) = ppms{i}(ymax_idx);

            P1(i) = params{i}(1);
            P2(i) = params{i}(2);
            P3(i) = params{i}(3);
            if numel(params{i}) == 3
                C(i) = 0;
            else
                C(i) = params{i}(4);
            end

            FWHM(i) = 2 * sqrt( P1(i)/((P1(i)/P3(i)+C(i))/2-C(i)) - P3(i));
                        
            A(i) = params{i}(1)/params{i}(3);
            yATxc(i) = C(i)+A(i);

            %%% Post-regression statistics
            ydata = ints{i};
            n = length(ydata);
            p = length(params{i});

            %%% R2
            % Compute the residual values as a vector signed numbers:
            yresid = ydata - yprime{i}; % == residual{i} == fitRes.residual
            % Square the residuals and total them obtain the residual sum of squares (Sum of Squared Residuals = SSR):
            SSResid = sum(yresid.^2); % == resnorm{i} == fitRes.resnorm // residual sum of squares (target optimization parameter!)
            % Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
            SSTotal = (n-1) * var(ydata);
            % Compute coefficient of determination - R2
            R2(i) = 1 - SSResid / SSTotal;
            % Compute "adjusted coeff of determination":
            R2adj(i) = 1 - (SSResid/SSTotal)*(n-1)/(n-p);

            %%% Chi2
            r = yresid;
            dof(i) = n-p;
            mse(i) = SSResid/dof(i);
            rmse = sqrt(mse(i)); % kind-of sigma (error) - since don't have actual sigma(error) for each point.
            chi2(i) = sum((yresid./rmse).^2);
            chi2red(i) = chi2(i)/dof(i); % In our case this == 1, since we don't have experimental errors.
            % gofP(i) = 1 - gammainc(chi2(i)/2,dof(i)/2); % F.Geier code:
            gofP(i) = 1-chi2cdf(chi2(i),dof(i)); % MIT code
            gofPred(i) = 1-chi2cdf(chi2red(i),dof(i));

            % GOF: P(X>=Chi2)
            clear ydata n p r yresid SSResid SSTotal rmse
        end


        %% Derive kex, Kop and dG
        %============================
        BF = S{ds}{1}.Acqus.BF1;
        if isfield(optns, 'temperature') && ~isempty(optns.temperature)
            var_T = optns.temperature; % K    
        else % read temperature from the audita.txt
            var_T = S{ds}{1}.Acqus.TE;
        end;
        
        FWHM(imag(FWHM)~=0) = nan; % replace all complex numbers with NaNs
%         lw{ds} = FWHM.*BF;
        lw{ds,iImino} = FWHM.*BF; % a_fit_LW_many
        fitted_ppm_center{ds,iImino} = P2; % a_fit_LW_many

        %% Plot and display fit results
        %============================
        if optns.displayFits
            rows = 5;
            columns = min([optns.n_spectra_to_disp, n_spectra]);
            sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

            FIG=FIG+1;
            figure(FIG); set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 0 columns*sbWidth rows*sbHeight]);
            SB = 0; % Reset numbering of subplots for new figure

            % xticks = linspace(pkBounds(2),pkBounds(1),4);

            for i=1:columns
            % columns here is a subset of n_spectra -- so columns(i) is
            % correctly mapped to "spectra"(i) (i.e. all values FWHM, etc
            % are correct.

                %%% Disp data and fit
                %========================    
                SB=SB+1;
                sb(SB);
                plot(ppms{i},ints{i},'o-');
                hold on;
                plot(ppms{i},yprime{i},'r','LineWidth',2);
                set(gca,...
                    'XDir','reverse'...
                    );
                axis tight;

            %     titleString = sprintf('%s | %.1f Hz',expt_names{ds}{i},FWHM(i)*S{i}.Procs.SF);
                titleString = {sprintf('%s  /  %.2f ppm  /  %s', dset_id{ds}, fitted_ppm_center{ds,iImino}(i), expt_names{ds}{i}),sprintf('LW=%.1f Hz', FWHM(i)*S{ds}{i}.Procs.SF)}; % TODO: should be BF here?!
                title(titleString,'FontSize',fSize,'Interpreter','None');

                %%% Disp residual
                %========================
                sb(SB+columns);
                plot(ppms{i},residual{i});
                set(gca,...
                    'XDir','reverse'...
                    );
                axis tight;

                %%% Disp summary of calcs
                %========================
                h2 = sb(SB+columns*2);
                plot([],[]);
                axis off;
                hold on;

                if ~flag_display_fit_details
                    results = strcat(...
                        sprintf('fitRange = %.1f %.1f',pkBounds_array{i}(1),pkBounds_array{i}(2)),{char(10)},...
                        sprintf('ymin | C = %.0f | %.0f',ymin(i),C(i)),{char(10)},...
                        sprintf('FWHM = %.3f | %.3f',man_FWHM,FWHM(i)),{char(10)},...
                        sprintf('ppm center = %.3f',fitted_ppm_center{ds,iImino}(i)),{char(10)},...
                        {char(10)},...
                        sprintf('Chi2red = %.3g',chi2red(i)),{char(10)},...
                        sprintf('R2adj = %.3f',R2adj(i)),{char(10)},...
                        {char(10)}...
                    );
                else        
                    results = strcat(...
                        sprintf('fitRange = %.1f %.1f',pkBounds_array{i}(1),pkBounds_array{i}(2)),{char(10)},...
                        sprintf('ymin | C = %.0f | %.0f',ymin(i),C(i)),{char(10)},...
                        sprintf('ymax | yATxc = %.2g | %.2g',ymax(i),yATxc(i)),{char(10)},...
                        sprintf('ymax-ymin | A = %.2g | %.2g',yspread(i),A(i)),{char(10)},...
                        sprintf('xAtymax | P2 = %.2f | %.2f',xATymax(i),P2(i)),{char(10)},...
                        sprintf('FWHM = %.3f | %.3f',man_FWHM,FWHM(i)),{char(10)},...
                        sprintf('ppm center = %.3f',fitted_ppm_center{ds,iImino}(i)),{char(10)},...
                        {char(10)},...
                        sprintf('Chi2 = %.2g',chi2(i)),{char(10)},...
                        sprintf('P(X>=Chi2) = %.2g',gofP(i)),{char(10)},...
                        sprintf('Chi2red = %.2g',chi2red(i)),{char(10)},...
                        sprintf('P(X>=Chi2red) = %.2g',gofPred(i)),{char(10)},...
                        sprintf('R2adj = %.3f',R2adj(i)),{char(10)},...
                        {char(10)}...
                    );
                end

                % Position text relative to the axes:
                xax = get(gca,'Xlim');
                yax = get(gca,'Ylim');
                text(min(xax),max(yax),results,'FontSize',fSize,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
                hold off;

            end
        end; % optns.displayFits
        
    end; % for-loop over iminos

end % fit linewidths for each dataset

d.time = time;
d.lw = lw;

if flag_save_to_scratch; save( fullfile(scratch_dir, sprintf('%s.mat', scratchfile_prefix)), 'd'); end;

end; % flag_refit_lw

if flag_load_from_scratch; load(fullfile(scratch_dir, sprintf('%s.mat', scratchfile_prefix))); end;

if flag_unittest_LW
    fprintf('\nTHIS WAS NOT TESTED IN a_fit_LW_many (worked fine in a_fit_LW)\n');
    UT_ref = load(fullfile(scratch_dir,'a_fit_LW_UT1.mat'));
    UT_ref = UT_ref.d;    
    UT_struct_fields(UT_ref,d);
end;


%% Display
%============================
FIG=FIG+1;
columns = n_iminos; % a_fit_LW_many
rows = 2;

fig_handle = figure(FIG); set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 0 columns*sbSide rows*sbSide]);
sb = @(x) subplot(rows,columns,x,'FontSize',fSize+4);

    
%%% Display Linewidth
%==============
if optns.display_lw
    SB = 0; % Reset numbering of subplots for new figure
    
    for iImino=1:n_iminos  % a_fit_LW_many      
        sb(iImino);
        
        data_to_plot = d.lw(:,iImino); % a_fit_LW_many
    
        plot(d.time{1}(1:end)./60,data_to_plot{1}(1:end), optns.plotSym, 'LineWidth', plotLW_list(1), 'Color', colors(1,:));

        hold on;
        axis tight;

        ylabel('Linewidth [Hz]', 'FontSize', fSize+4);
        xlabel('time [h]', 'FontSize', fSize+4);

        ymax_for_all_sets = max(cell2mat(cellfun(@(x) max(x), data_to_plot, 'un', 0)));    
        ylim([0 ymax_for_all_sets*1.2]);

        xmax_for_all_sets = max(cell2mat(cellfun(@(x) max(x), d.time, 'un', 0)))./60;    
        if xmax_for_all_sets < optns.xaxis_start
            disp('xmax_for_all_sets lower than optns.xaxis_start, shifting axis window by 2h.');
            optns.xaxis_start = xmax_for_all_sets-2;
        end
        xlim([optns.xaxis_start xmax_for_all_sets]);

        if n_sets > 1
            for ds=2:n_sets
                plot(d.time{ds}(1:end)./60,data_to_plot{ds}(1:end), optns.plotSym, 'LineWidth', plotLW_list(ds), 'Color',colors(ds,:));
            end
        else
        end
        hold off;

        if isfield(optns, 'legend') && ~isempty(optns.legend) && iImino==1
            axP = get(gca,'Position'); 
            legend(optns.legend, 'Location', 'SouthOutside', 'FontSize', fSize+1,'Interpreter','None');
            set(gca, 'Position', axP);
        end;
                
        titleString = sprintf('imino at %.2f ppm', imino_ppm_list(iImino));
        title(titleString);
        
    end; % loop over iminos

end; % optns.display_lw

if optns.save_figure
    file_name = sprintf('%s_all_%s.pdf', strjoin(dset_id,'_'), datestr(now,'yymmdd'));
    save_figure(fig_handle,file_name);
end;    

end % main function

function time0 = get_time0(data_dir,dataset_name)
    notes_path = fullfile(data_dir,dataset_name,'notes.txt');
    % <time>2016-08-16 13:08:30</time0>
    filetext = fileread(notes_path);
    pattern = '<time0>(.*?)</time0>';
    time0 = regexp(filetext,pattern,'tokens');
    assert(numel(time0) == 1, 'Abort: more than one time0 matches found in notes.txt.');    
    time0 = char(time0{:}); % convert from cell to a string % 2017-07-19: added char()
    clear filetext;
end


function save_figure(fig_handle,file_name)
%% Save figure as pdf
%=============================================================
    set(fig_handle,'Units','Inches');
    pos = get(fig_handle,'Position');
    set(fig_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    figure_export_dir = 'figure_images';
    if ~exist( figure_export_dir ,'dir'), mkdir(figure_export_dir), end;
    print(fig_handle,fullfile(figure_export_dir, file_name),'-dpdf','-r1200')
end