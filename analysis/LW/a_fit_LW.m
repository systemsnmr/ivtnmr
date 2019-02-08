function a_fit_LW(nmr_data_path, dsets, imino_ppms, optns)
% !!! ITS GOOD TO ANALYZE ONE SET - AND THEN SAVE LW TO IVTNMR !!!
% Function for visualizing - make separate!

%%% Obligatory parameters
if nargin == 0
    clear; close all;
    
    analysis_root = fullfile(pwd, '..');
    nmr_data_path = fullfile(analysis_root, '..', 'spectra');
%     nmr_data_path = '/Volumes/Data/yar/_eth2/data_NMR/spectra/';
                
    dsets = {...
        '170322_IN75a_EV2_co-UP1_303K_600'
        '180314_IN72b_SMN214_co-NUP1_303K_600'
    };
    imino_ppms = [13.18, 14.033, 13.518, 13.197, 11.937, 11.623];        
end;

% Defaults for optional parameters (are set all together at the moment)
if any(nargin == [0 3])
    optns.fitWidth = 0.1;
    
    optns.display_lw = 1;    
    optns.xaxis_start = 5; % h
    optns.plotSym = 'o-';
    
    optns.display_dG = 1; %% -- better exclude?!! -- this should be a separate analysis!
    optns.temperature = 303; % in case if setting in the acqus is not exact.    
    
    optns.save_figure = 1;
%     optns.save_to_global_ivtnmr

    optns.displayFits = 0;
    optns.n_spectra_to_disp = 5; % Spectra shown when showing LW fits.    
end;

%% Fixed Settings
%=====================
iminoPPM = imino_ppms(1); % takes first imino by default

fSize = 9;
scale = 1;
sbSide = 220;
sbWidth = sbSide * scale;
sbHeight = sbSide * scale;

display_fit_details = 0;
global FIG; FIG=0;

addpath( genpath(fullfile(pwd,'lib')) );

n_sets = numel(dsets);
colors = lines(n_sets);
plotLW_list = ones(1, n_sets);

%%% Preallocate to store variables for several sets:
% (use cells, cuz different dsets may h diff amount of data
lw = cell(n_sets,1);
kex = cell(n_sets,1);
kex_unfold = cell(n_sets,1);
kex_unfold_mean = cell(n_sets,1);
Kop = cell(n_sets,1);
dG_folding = cell(n_sets,1);
dG_folding_mean = cell(n_sets,1);
time = cell(n_sets,1);

%% Actual fitting
%=====================

for ds=1:n_sets    
    files = cell2mat(arrayfun(@(x) str2num(x.name), dir(fullfile(nmr_data_path, dsets{ds},'6*')), 'un', 0));
    expnos = files(files>=6000 & files<=6999)';

    time0 = get_time0(nmr_data_path, dsets{ds});
    time{ds} = getNMRTime(dsets{ds},expnos,time0); % time = getNMRTime(datasetName,experiments,time0)

    % TODO:
    % expt_names - derive from TIME
    expt_names = arrayfun(@num2str, expnos, 'unif', 0);

    rbSpectra = cellfun(@(x) fullfile(nmr_data_path, dsets{ds},x,'pdata/1/1r'), ...
        arrayfun(@num2str, expnos, 'un', 0), 'un', 0);

    pkBounds = [iminoPPM+optns.fitWidth iminoPPM-optns.fitWidth];

    %% Import spectra
    %==============================
    n_spectra = numel(expnos);

    S = rbnmr(rbSpectra);

    % Cut out the part for fitting:
    cutIds = cellfun(@(x) find(x.XAxis < pkBounds(1) & x.XAxis > pkBounds(2)), S, 'UniformOutput', false);
    ints = cellfun(@(source,extract) source.Data(extract), S, cutIds, 'UniformOutput', false);
    ppms = cellfun(@(source,extract) source.XAxis(extract), S, cutIds, 'UniformOutput', false);


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
    %     [yprime{i} params{i} resnorm{i} residual{i}] = lorentzfit(ppms{i},intensities{i},[],[]);
    %     [yprime{i} params{i} resnorm{i} residual{i}] = lorentzfit2(ppms{i},ints{i},P0,[]);
        [yprime{i} params{i} resnorm{i} residual{i}] = lorentzfit2(ppms{i},ints{i},P0,[],'3'); % zero baseline (no constant term

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
        %%% Simplified version of stuff from imino_kex_SMN210_WESF_full_v9f_new_SMN29_110_CUG.m
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
        %%% This part was redefined slightly now (compared to
        %%% imino_kex_SMN210_WESF_full_v9f_new_SMN29_110_CUG.m)
        r = yresid;
        dof(i) = n-p;
        mse(i) = SSResid/dof(i);
        rmse = sqrt(mse(i)); % kind-of sigma (error) - since don't have actual sigma(error) for each point.
        chi2(i) = sum((yresid./rmse).^2);
        chi2red(i) = chi2(i)/dof(i); % In our case this == 1, since we don't have experimental errors.
        %%% As in F.Geier's code:
        % gofP(i) = 1 - gammainc(chi2(i)/2,dof(i)/2);
        %%% MIT code
        gofP(i) = 1-chi2cdf(chi2(i),dof(i));
        gofPred(i) = 1-chi2cdf(chi2red(i),dof(i));

        % GOF: P(X>=Chi2)

        clear ydata n p r yresid SSResid SSTotal rmse
    end


    %% Derive kex, Kop and dG
    %============================
    BF = S{1}.Acqus.BF1;
    if isfield(optns, 'temperature') && ~isempty(optns.temperature)
        var_T = optns.temperature; % K    
    else % read temperature from the audita.txt
        var_T = S{1}.Acqus.TE;
    end;
    
    FWHM(imag(FWHM)~=0) = nan; % replace all complex numbers with NaNs
    lw{ds} = FWHM.*BF;
    kex{ds} = lw{ds}*pi; % k = LW*pi = (FWHM*B0 field)*pi
    kex_base_flip_plus_R2 = 61.4; % Hz: STSL1 TTD 28mM Mg 303K (see lorentzian_v04a_SMN.m)
    
    kex_unfold{ds} = kex{ds} - kex_base_flip_plus_R2;

    % % TTRED - with 50mM Arg/Glu
    % kex_int_Ura_TTRED = 2094645.5; 
    % kex_int_Gua_TTRED = 1959861.5;
    % TTD - w/o Arg/Glu
    kex_int_Ura_TTD = 1042983.7;
    kex_int_Gua_TTD = 956480.3;
    kex_int = mean([kex_int_Ura_TTD, kex_int_Gua_TTD]);

    % This should be an under-estimate of dG - because the LW of ST2 in PO4/KCl
    % is narrower than it would have been in TTD
    Kop{ds} = kex_unfold{ds} ./ kex_int;
    const_R = 1.9858775; % cal K-1 mol-1
    dG = -const_R .* var_T .* log(Kop{ds}) ./ 1000; % kcal/mol
    dG_folding{ds} = -dG;
    dG_folding_mean{ds} = mean(dG_folding{ds});
    
    kex_unfold_mean{ds} = mean(kex_unfold{ds});
    
    % Need to include removal of T2?!
    
    %% Plot and display fit results
    %============================
    if optns.displayFits && ds==1
        rows = 5;
        columns = n_spectra_to_disp;
        sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

        FIG=FIG+1;
        figure(FIG); set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 0 columns*sbWidth rows*sbHeight]);
        SB = 0; % Reset numbering of subplots for new figure

        % xticks = linspace(pkBounds(2),pkBounds(1),4);

        for i=1:n_spectra_to_disp

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

        %     titleString = sprintf('%s | %.1f Hz',expt_names{i},FWHM(i)*S{i}.Procs.SF);
            titleString = {expt_names{i},sprintf('LW*pi=%.1f Hz', FWHM(i)*S{i}.Procs.SF*pi)}; % TODO: should be BF here?!
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

            if ~display_fit_details
                results = strcat(...
                    sprintf('fitRange = %.1f %.1f',pkBounds(1),pkBounds(2)),{char(10)},...
                    sprintf('ymin | C = %.0f | %.0f',ymin(i),C(i)),{char(10)},...
                    sprintf('FWHM = %.3f | %.3f',man_FWHM,FWHM(i)),{char(10)},...
                    sprintf('kex (FWHM*BF*pi) = %.3f',kex{ds}(i)),{char(10)},...
                    sprintf('dG_folding = %.2f',dG_folding{ds}(i)),{char(10)},...
                    {char(10)},...
                    sprintf('Chi2red = %.3g',chi2red(i)),{char(10)},...
                    sprintf('R2adj = %.3f',R2adj(i)),{char(10)},...
                    {char(10)}...
                );
            else        
                results = strcat(...
                    sprintf('fitRange = %.1f %.1f',pkBounds(1),pkBounds(2)),{char(10)},...
                    sprintf('ymin | C = %.0f | %.0f',ymin(i),C(i)),{char(10)},...
                    sprintf('ymax | yATxc = %.2g | %.2g',ymax(i),yATxc(i)),{char(10)},...
                    sprintf('ymax-ymin | A = %.2g | %.2g',yspread(i),A(i)),{char(10)},...
                    sprintf('xAtymax | P2 = %.2f | %.2f',xATymax(i),P2(i)),{char(10)},...
                    sprintf('FWHM = %.3f | %.3f',man_FWHM,FWHM(i)),{char(10)},...
                    sprintf('kex (FWHM*BF*pi) = %.3f',kex{ds}(i)),{char(10)},...
                    sprintf('dG (RTln(kunf/kint)) = %.2f',dG(i)),{char(10)},...
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
    end % optns.displayFits
    
    
end % fit datasets    


%% Display
%============================
FIG=FIG+1;
columns = 2;
rows = 2;

fig_handle = figure(FIG); set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 0 columns*sbSide rows*sbSide]);
sb = @(x) subplot(rows,columns,x,'FontSize',fSize+4);

    
%%% Display kex
%==============
if optns.display_lw
    SB = 1; % Reset numbering of subplots for new figure
    sb(SB);

    data_to_plot = lw;
    
    plot(time{1}(1:end)./60,data_to_plot{1}(1:end), optns.plotSym, 'LineWidth', plotLW_list(1), 'Color', colors(1,:));

    hold on;
    axis tight;

    ylabel('Linewidth [Hz]', 'FontSize', fSize+4);
    xlabel('time [h]', 'FontSize', fSize+4);
    
    ymax_for_all_sets = max(cell2mat(cellfun(@(x) max(x), lw, 'un', 0)));    
    ylim([0 ymax_for_all_sets*1.2]);

    xmax_for_all_sets = max(cell2mat(cellfun(@(x) max(x), time, 'un', 0)))./60;    
    xlim([optns.xaxis_start xmax_for_all_sets]);

    if n_sets > 1
        for ds=2:n_sets
            plot(time{ds}(1:end)./60,data_to_plot{ds}(1:end), optns.plotSym, 'LineWidth', plotLW_list(ds), 'Color',colors(ds,:));
        end
    else
    end
    hold off;

end % 


if optns.display_dG
    %%% Display dG
    %==============
    %%% Used to plot 2:end in all! (changed 2017-04-06)
    
    SB = 2; % Reset numbering of subplots for new figure
    sb(SB);

    plot(time{1}(1:end)./60,dG_folding{1}(1:end), optns.plotSym,'Color',colors(1,:));
    hold on;
    axis tight;

    yax = get(gca,'Ylim');
    ylim([mean(yax)-0.5 mean(yax)+0.5]);
    
    xax = get(gca,'XLim');
    xlim([optns.xaxis_start xax(2)]);
    
    ylabel('RNA stability, dG [kcal/mol]', 'FontSize', fSize+4);
    xlabel('time [h]', 'FontSize', fSize+4);

    if n_sets > 1
        for ds=2:n_sets
            plot(time{ds}(1:end)./60,dG_folding{ds}(1:end), optns.plotSym,'Color', colors(ds,:));
        end
    else
        hold off;
    end
    
%         title_text = sprintf('%s (%.3g ppm)\n dG mean=%.3g kcal/mol', expt_title, iminoPPM, dG_folding_mean);
    title_text = sprintf('dG (imino %.2f ppm) for', iminoPPM);
    title({title_text, '(GU)UA(UA), TTD-Mg pH7.7, 30nt RNA', 'see code: dG deriv depends', 'on buffer, RNA size,', 'basetype and neighbors'}, 'FontSize', fSize+1);

    if isfield(optns, 'legend') && ~isempty(optns.legend)
        legend(optns.legend, 'Location', 'Best', 'FontSize', fSize+3);
    end
    
%     disp(dG_folding_mean);
         
end % display_dG

if optns.save_figure
    file_name = sprintf('%s_%.3f_ppm_%.2f_fit.pdf', dsets{ds}, iminoPPM, optns.fitWidth);
    save_figure(fig_handle,file_name);
end;    

end % main function


function time0 = get_time0(data_dir,dataset_name)
    notes_path = fullfile(data_dir,dataset_name,'notes.txt');
    % <time>2016-08-16 13:08:30</time0>
    filetext = fileread(notes_path);
    pattern = '<time0>(.*?)</time0>';
    time0 = regexp(filetext,pattern,'tokens','once');
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


