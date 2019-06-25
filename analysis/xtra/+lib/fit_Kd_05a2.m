% Based on nmrkdCombined
% 2016-11
%  - Added lsq fitting of the function and output of params
%  - Adjusted simulation to display the output Kd from the fits.
%  - Varying normalization of BOTH L and P to see if back-calc Kd is affected
% Based on nmrkdCombined_sim_UP1_conc_norm.m
% 2017-03-20
%  - Modified to use both P and L as fixed vectors during optimization!
%   (sorta NMR/ITC titration when sample gets diluted)
% 2017-04-02 - Based on nmrkdCombined_170320_var_P_and_L_NMR_fit.m
%  - Adapted to fit data from IVTNMR series
% 2017-08-17 - fit_Kd_01 - based on nmrkdCombined_170405_fit_all
%  - Optimized code to make it an externally-called function

% 2017-08-18 fit_Kd_02
%  - added initial guesses for d_sat (max_y) and Kd(max_x*0.5) params (p n
%  find what PRISM does - just assumed "good" data).
%  - added ERROR calculation

% 2017-08-24 fit_Kd_03 (changes input - doing version up):
%  - Added ability to use fixed_d_sat (checking dsat input too)
%   (expects either 'auto_SCALEVALUE' e.g. 'auto_2.x', or numerical == # of peaks)

% - fit_Kd_04b_no_dsat - added flag to not add csp_fit_long_obs if dSat is
% not fixed

%   Based on fit_Kd_04b_no_dsat.m:
% fit_Kd_05a - generalizing for other use:
%   - removed fixed labels of axes
% fit_Kd_05a2 (for github ivtnmr deposition):
%   - adds optns structure
%   - can define units in the optns

% TODO
% v0X (changes input args - version up)
%  - Include options argument
%  - apart from ligand_conc, target_csp, target_conc - put other params
%  into OPTIONS structure
%  - add input checks


% Consider setting Kd as FIRST parameter in the fitting - then can simplify
% downstream stuff when dsat is fixed - popt(1) will always be Kd!

function o = fit_Kd_05a2(ligand_conc, target_csp, target_conc, plot_results, fix_dsat)
% Can fit multiple peaks with same ligand_conc -- pass them as rows in target_csp argument.
% For linearly-behaving CSP - set fix_dsat >= 2*max_csp

global FIG;    
global SB;

if nargin == 0
    %%% Import all integrated data:
    % - PO4/NTP/RNA concentrations, as converted from 31P integrals by create_data_str.m
    % - HSQC chemical shifts - extrapolated for common time-vector
    data_struct_file = '/Volumes/Data/yar/Dropbox/Science/Teaching/2017-03_NMR_BlockKurs/BK2017/data_modeling/170331_IN60_IN70/v6a_nD_oM_HSQC/data_nmr_selected_G6.mat'; 

    load(data_struct_file);

    rna_conc = sort( dsel.y(3,:) ); % TMP SORT here is a hack - couple of values are going down in concentration, which causes ODE complaints.

    target_csp = dsel.hsqc.mean'; % Convert to row-vector to match analyses below

    rna_conc = rna_conc-rna_conc(1); % Normalize to have 0uM RNA at the start of reaction.
    rna_conc = rna_conc.*1000./4; % Convert mM-to-uM (1000), and convert 31P atoms to 4nt-long RNA.

    ligand_conc = rna_conc;

    target_conc = 150; % IVTNMR fixed
    
    plot_results = 1;
    
    fix_dsat = 'auto_2';    
%     fix_dsat = ones( size(target_csp,1) ,1); % testing numeric version
    fix_dsat = [];
    
end

dsatIsFixed = exist('fix_dsat','var') && ~isempty(fix_dsat);

DEBUG = 0;

%% Input checks
%=========================
n_points = numel(ligand_conc);
n_csp_peaks = size(target_csp,1);

%%% Checking input dsat value
if dsatIsFixed
    if ischar(fix_dsat) && strcmp(fix_dsat(1:4),'auto')
        fprintf(1,'indeed auto\n');
        fixed_dsat_arr = max(target_csp,[],2) .* str2num(fix_dsat(6:end)); % potentially for input of several peaks
        
    elseif isnumeric(fix_dsat)
        assert(n_csp_peaks == numel(fix_dsat), 'Error: # of fix    _d_sat values does not match # of peaks!');
        fixed_dsat_arr = fix_dsat;
        
    else
        error('Expecting fix_dsat either auto_x.x or numeric matching # of peaks');
    end
end % dsatIsFixed


if numel(target_conc) == 1
    fprintf(1,'Input target concentration is a scalar, not vector.\n');
    fprintf(1,'Assuming its constant (%d) during titration.\n',target_conc);
    repmat(target_conc,1,n_points);
end


%% Main script
%=========================

% Display settings
%=================================================
if plot_results
    sbSide = 180; sbWidth = sbSide; sbHeight = sbSide*1.5; fSize = 12;
    rows = 2; columns = n_csp_peaks;
    SB=0; % global defined at the top of the function
    sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

    % global defined at the top of the function
    if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
    FIG = FIG+1;
    set(figure(FIG), 'Color', repmat(1,1,3), 'Position', [0 0 columns*sbWidth rows*sbHeight]) %#ok<RPMT1> % [x y width height]
end

% Params for REFERENCE. Initial GUESS params are more below
%=================================================
ref_dsat = 0.1; % delta_sat
ref_kd = 100; % Kd
x0 = 0; % Initial PL complex

%%% Function for simulating results
%=================================================
syms PL delta_sat2 Kd2 L2 Etot2;
PL = delta_sat2 * ((Kd2 + L2 + Etot2) - sqrt((Kd2 + L2 + Etot2)^2 - 4*L2*Etot2)) / (2*Etot2);

o(n_csp_peaks,1) = struct();

for pk=1:n_csp_peaks    

    %%% Least-squares (minimization of residuals) fit of Kd model (equation) to data
    %================================================
    ydata = target_csp(pk,:); % protein chemical shift
    if dsatIsFixed
        fixed_dsat_value = fixed_dsat_arr(pk);
    end
    xdata = ligand_conc; % L - ligand conc
    hasSameSizeXYpoints = numel(ydata)==numel(xdata);
    assert(hasSameSizeXYpoints, 'Error: number of X and Y data points must match');

    % Initial parameter guesses
%     p0 = [1 1]; % [dsat Kd] init gues - either set close, or increase MaxFunEvals, MaxIter, TolFun, ...
    if dsatIsFixed
        p0 = [max(xdata)*0.5]; % [Kd] init gues - either set close, or increase MaxFunEvals, MaxIter, TolFun, ...
    else
        p0 = [max(ydata) max(xdata)*0.5]; % [dsat Kd] init gues - either set close, or increase MaxFunEvals, MaxIter, TolFun, ...
    end

    % For comparison of lsqcurvefit and lsqnonlin implementations see
    % /Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra/140815_SMN210_pH_IN_refs_303K_600/SMN210_imino_intens_vs_pH_lsqnonlin_fit.m
    % anonymous Kd func
    % using ANALYTICAL solution
    if dsatIsFixed
        analytical_kd_func = @(p,xdata) fixed_dsat_value * ((p(1) + xdata + target_conc) - sqrt((p(1) + xdata + target_conc).^2 - 4 .* xdata .* target_conc)) ./ (2*target_conc);
    else
        analytical_kd_func = @(p,xdata) p(1) * ((p(2) + xdata + target_conc) - sqrt((p(2) + xdata + target_conc).^2 - 4 .* xdata .* target_conc)) ./ (2*target_conc);
    end

    options = optimset('MaxFunEvals',1e5,'MaxIter',1e4,'TolFun',1e-7);

    [popt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(analytical_kd_func,p0,xdata,ydata,[],[],options);

    % Residual values as a vector signed numbers:
    yresid = ydata - analytical_kd_func(popt,xdata); % == 'residual' in the output
    % Square the residuals and total them obtain the residual sum of squares:
    SSresid = sum(yresid.^2); % == 'resnorm' - residual sum of squares (was optimized)
    % Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
    SStotal = (length(ydata)-1) * var(ydata);
    % Compute coefficient of determination - R2
    R2 = 1 - SSresid / SStotal;
    % Compute "adjusted coeff of determination" (now correct, no braket missing)
    % 1-(SSR/SST)*(n-1)/(n-p) // if p includes constant term. if it doesn't - use (n-p-1) in denominator)
    R2adj = 1 - (SSresid / SStotal) * (length(ydata)-1)/(length(ydata)-length(popt));

    if nargin==0 || DEBUG
        if dsatIsFixed
            disp('----- Initial(p) vs fitted(popt) parameters');
            disp({'d_sat(fix)', 'Kd', 'Etot', 'R2adj';...
                ref_dsat, ref_kd,'var', '';...
                fixed_dsat_value, popt(1), 'var', R2adj});
        else
            disp('----- Initial(p) vs fitted(popt) parameters');
            disp({'d_sat', 'Kd', 'Etot', 'R2adj';...
                ref_dsat, ref_kd,'var', '';...
                popt(1), popt(2), 'var', R2adj});
        end
    end;
            
    %%% Post-regression statistics
    %==================================
    %=======================================================================%
    %%%% SPIN THIS OUT INTO A SEPARTE FUNCTION WHEN REUSING NEXT TIME!!! %%%%
    %=======================================================================%
    
    %%% Simplified / combined versions of:
    % lorentzian_v04e_SMN.m (that was simplified version of imino_kex_SMN210_WESF_full_v9f_new_SMN29_110_CUG.m)
    % and
    % 2015-02_SE_CI_SD_from_NLS_fits/Kd_fit_v2_using_normalyzed_Cov_matrix_for_prism_calc
    
    if 1
    n_data = length(ydata);
    n_par = length(popt);
    n_obs = 1; % @Y: here have only single observable
        
    %%% Quality of the fit
    %==
    %%% This part was redefined slightly in lorentzian_v04e_SMN compared to imino_kex_SMN210_WESF_full_v9f_new_SMN29_110_CUG
    dof = n_data-n_par;
    mse = SSresid ./ dof;
    % A complicated way to do the same - from MathWorks:
	% mse = sum(abs(r).^2)/(size(full(J),1)-size(full(J),2)); % Mean Squared Error = Sum of Squares / N
    rmse = sqrt(mse); % kind-of sigma (error) - since don't have actual sigma(error) for each point.
    chi2 = sum((yresid./rmse).^2); % sum of residuals normalized by mean error
    chi2red = chi2 ./ dof; % In our case this == 1, since we don't have experimental errors.
        
    if ~isreal(chi2)
        fprintf(1,'\n pk=%i \n',pk);
        warning('Fit returns complex optimized parameters. Taking only real part for gofP & gofPred calc.');
        chi2 = real(chi2);
        chi2red = real(chi2red);
    end
    gofP = 1-chi2cdf(chi2,dof); %%% MIT code
    gofPred = 1-chi2cdf(chi2red,dof); %%% MIT code

    %%% Error calc - MathWorks-MSE (also matches results of internal PRISM fits). See
    %==
    %%% MathWorks-MSE (also matches results of internal PRISM fits). See
    %%% "Parameter confidence interval estimation / err"
    % http://ch.mathworks.com/matlabcentral/newsreader/view_thread/86198
    [Q,R] = qr(full(jacobian),0);
    Rinv = inv(R);
    Sigma = Rinv*Rinv'*mse; % ?C. Sorta "Covariance matrix" - scaled by MSE.
    
    % Assuming sqrt(diag(C)) = SD
    sd = sqrt(diag(Sigma));
    sem_www1 = sd' ./ sqrt(n_data*n_obs);
    sd_www1 = sem_www1 .* sqrt(n_data*n_obs);

    level = 0.95;    
    % sd*fac matches PRISM-based output closer than sem*fac !
    % In principle - degrees of freedom should be taken into account here !!
    fac = sqrt(2)*erfinv(2*(1-(1-level)/2)-1); % is 1.96 for level = 0.95
    confi_www1 = sd_www1*fac; % was sem_www1*fac, but seems underestimating
    
    end % error calc

    
    %%% Simulate fitted solution
    %===================================
    xscale = ligand_conc;
    if dsatIsFixed
        fit_result = subs(PL,[delta_sat2 Kd2 Etot2 L2],[fixed_dsat_value popt(1) target_conc {xscale}]);
%         fit_long_obs = linspace(0,10,30);
        fit_long_obs = [xscale, linspace(xscale(end),3,200)];
        fit_long = subs(PL,[delta_sat2 Kd2 Etot2 L2],[fixed_dsat_value popt(1) target_conc {fit_long_obs}]);
    else
        fit_result = subs(PL,[delta_sat2 Kd2 Etot2 L2],[popt(1) popt(2) target_conc {xscale}]);
    end
    

    %%% Plot data and final solution
    %================================================
    if plot_results
        SB=SB+1; sb(SB);

        % Data
        plot(xdata,ydata,'o','MarkerSize',9);
        hold on;

        % Plot fitted solution
        plot(xscale,fit_result,'--r','LineWidth',2);
        axis tight;
        
        if nargin == 0
            fit_name = dsel.hsqc.names{pk};
        else
            fit_name = num2str(pk);
        end
%         title_string = sprintf('%s\nKd = %.0f ?mM', fit_name, popt(2));
        title_string = fit_name;
        if dsatIsFixed
            title_string = strcat(title_string, sprintf('\nfix_dsat=%.3f', fixed_dsat_value));
            kd_val = popt(1);
            kd_sd = sd_www1(1);
            kd_ci = confi_www1(1);
        else
            title_string = strcat(title_string, sprintf('\ndsat=%.3f (%.3f)', popt(1), sd_www1(1)));
            kd_val = popt(2);
            kd_sd = sd_www1(2);
            kd_ci = confi_www1(2);
        end
        title_string = strcat(title_string, {char(10)});
        title_string = strcat(title_string, sprintf('\nKd=%.2g (%.2g) ?mM', kd_val, kd_sd));
        title_string = strcat(title_string, sprintf('\nKd CI = %.1f to %.1f', kd_val-kd_ci, kd_val+kd_ci));
        title_string = strcat(title_string, {char(10)});
        title_string = strcat(title_string, sprintf('\nR2=%.3f',R2));
        title_string = strcat(title_string, sprintf('\nR2adj=%.3f',R2adj));
        
        title(title_string,'FontSize',fSize, 'Interpreter', 'None');

        if pk==1
            xlab_string = {'ligand [?mM]'};
            legend({'data','fit'},'Location','SE');
        else
            xlab_string = 'ligand [?mM]';
        end

        xlabel(xlab_string);
        ylabel('\deltaCS (H+N)');
    end % plot_results
    
    %%% Assemble structure for output
    o(pk).popt = popt;
    if dsatIsFixed
        o(pk).pnames = {'Kd'};
    else
        o(pk).pnames = {'d_sat', 'Kd'};
    end
    o(pk).R2 = R2;
    o(pk).R2adj = R2adj;
    o(pk).csp_fit = fit_result;
    o(pk).sd = sd_www1;
    o(pk).sem = sem_www1;
    o(pk).confi = confi_www1;
    
    o(pk).resnorm = resnorm;
    o(pk).residual = residual;
    o(pk).exitflag = exitflag;
    o(pk).output = output;
    o(pk).lambda = lambda;
    o(pk).jacobian = jacobian;
    o(pk).SSresid = SSresid;
    o(pk).SStotal = SStotal;

    o(pk).dof = dof;
    o(pk).mse = mse; % SSresid/dof
    o(pk).rmse = rmse; % sqrt(mse)
    o(pk).chi2 = chi2; % sum((yresid./rmse).^2)
    o(pk).chi2red = chi2red; % chi2/dof
    o(pk).gofP = gofP;
    o(pk).gofPred = gofPred;        

    if nargin ~= 0 && dsatIsFixed % currently the fit_long is only generated for whed dsat is fixed (as in external call)
        o(pk).csp_fit_long_obs = fit_long_obs;
        o(pk).csp_fit_long = fit_long;
    end
    
end % n_csp_peaks

if nargin == 0
    disp('o(1) example');
    disp(o(1));
end;

end % function