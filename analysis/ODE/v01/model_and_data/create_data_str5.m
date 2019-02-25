function create_data_str5(dset,use_average_PO4,integr_ref,export_files, hn_optns, optns)
% getNMRData4 - now can estimate integration error (if provided with noise_region)

%% CHANGELOG
%===========================
%%% 2017-06-16. @YAR. - create_data_str2
%   - Removed all old / unused parts - pH-determination from 1D1H!!
%   - removed import of iminos completely
%   - set hsqc import as optional (if not imported - is placed as NAN into
%   matrix)
%   - somewhat automated determination of available expnos

%%% 2017-07-07. @YAR. - create_data_str3
%   - Enabled running from external call:
%       - optional arguments
%       - check if running internally (nargin).
%   - Check that AUTO-detected EXPNOS are thousands i.e. in range a(a >= 5000 & a <= 5999)


%%% 2018-05-03. @YAR. - create_data_str4
%   - Adds import of HSQC data from specified ivtnmr_path
%   - Interpolates CSP to the 31P vector
%   - CSP NORMALIZATION to max_CSP is done in the MODEL FILE!
%   >>> SEE EVN "IVTNMR IN48a and other analyses / model" - A 12


%%% 20190204 @YN - create_data_str5
%   - Adds 'optns' input parameter - e.g. to set datasave dir.
%   - !! TODO: drop 'pure' versions - just do PNM.
%   20190220 - Removes hard-coded CURDIR paths - uses optns.nmr_data_path now
%   20190225
%   - makes use 31P 'PNM' file when doing 'intern' (any(strcmp(integr_ref,{'PNM', 'intern'}))


% TODO:
% - Automatic detection of expt range:
%    - Move into helper function
%    

% - Automatic decision on minimal data vector length
% - Make such that create data str can work with single or multiple
% datasets (for multiple then use error calculation, for single - skip)
% - check that dataVectorLength is not crucial for iminos and etc (as of v5a just
% dropping the end-points which have NAN after interpolation.

%% Script Settings
%==================================

addpath(genpath('lib'));

fSize = 12;
scale = 1;
sbWidth = 130 * scale;
sbHeight = 130 * scale;
% global SB; SB=0;
global FIG; FIG=0;

plotSym = 'o-'; % symbol in the graph

disp_PO4_Mg_final_conc = 1;

% selectDataForDisplay = [1 2 3 4 5 6]; % Comment out to plot all.

% % Options to save data to mat file (to speed up re-analyses):
% saveImportedDataToMatFile = 1;
% reimportData = 1;

if nargin == 0
    fprintf(1,'%s requires at least 3 input arguments. Running test example.\n\n', mfilename);
    % Settings to be varied by user if running directly from script
    dset = {...
%         '170404_IN72i_SMN2_post-UP1_303K_600';
%         '170308_IN60a_pR02_co-A1R1_303K_600';
        '170310_IN70a_pR02_co-UP1_303K_600';
    };

    integr_ref = 'PNM'; % PNM or pure
    integr_ref = 'pure'; % PNM or pure
    integr_ref = 'intern'; % PNM or pure or intern (uses internal then)
    
    use_average_PO4 = 1; % overrides use_external_PO4ref
    
    export_files = 1;
    
    hn_optns.ivtnmr_path = '/Volumes/Data/yar/Dropbox/_eth2/data_MathModeling/170821_IN6x_IN7x_CSP_n_Kd_fits/datasave_03c/ivtnmr_IN70a.mat';    
    hn_optns.hn_peak_names = {'H33', 'R75'};
    hn_optns.p_conc = 0.15; % mM
    
    % Visualization
    plot_hsqc = 1;
    plot_p31 = 1;
    plot_final = 1;
    
else
    plot_hsqc = optns.plot_hsqc;
    plot_p31 = optns.plot_p31;
    plot_final = optns.plot_final;        
end

if ~exist('export_files','var') || isempty(export_files)
    export_files = 0;
end

dsetname = regexp(dset{1}, '_', 'split');
dsetname = dsetname{2};

%% Data settings
%==================================
refPeak = 7; % aNTP - peak used to reference concentration from integral.
refConc = 20; % mM. concentration to which the ref peak should be set
useT1Corrections = 1; % to have more accurate concentration estimates (T1 of PO4 is longer than of RNA).
tot_phosph_in_NTPs = 3*20; % mM of total PO4 - from NTPs

set_RNA_T1_to_1 = 0; % not using - was just for tests
use_external_PO4ref = 1; % If this is '1' - then PO4 is not re-normalized additionally agains aNTP
% .. i.e. keeps normalization to external reference.

if strcmp(integr_ref, 'intern')
    use_average_PO4 = 0;
    use_external_PO4ref = 0;
end

noise_31P = [-24,-28];

if exist('hn_optns', 'var') && ~isempty(hn_optns)
    include_hsqc = 1;
else
    include_hsqc = 0;
end

err_calc = 1; % 1: std of 1st and 2nd pt. 2: uniform (max) std from (1). 3: noise-based error
% was: use_uniform_std

if any(strcmp(integr_ref,{'PNM', 'intern'}))
    p31_integr_file = 'integr_results_31P_PNM.txt'; % integr_results_31P_pure_PO4.txt. default '31P'
else
    p31_integr_file = 'integr_results_31P_pure_PO4.txt'; % integr_results_31P_pure_PO4.txt. default '31P'
end


%% Data variables
%=======================================================
CURDIR = optns.nmr_data_path;
% CURDIR = '~/Desktop/BK2017/NMR';

dset_t0 = cellfun(@(x) get_time0(optns.nmr_data_path,x), dset); % uses helper func (see EOF)

dset_exp_31P = {...
    cell2mat(arrayfun(@(x) str2num(x.name), dir(fullfile(CURDIR,dset{1},'5*')), 'unifo', 0))';
    };

% Filter out only experiments > 5000 - i.e. drop any 5, 50, 500
dset_exp_31P{1} = dset_exp_31P{1}(dset_exp_31P{1} >= 5000);

numOfSets = numel(dset);

% ====== Define shortest time vector to use as common time-base ======
% Vectors passed to interpl need to be of same length, so have to set:
% 1) Expt set with the smallest number of full expt series:
refSet = 1;
% 2) Which experiment was running first in [expt series] - for interp1 to work.
refExpt = dset_exp_31P;

% Manually set fixed data vector length
% v5a - decided to omit this, now just using NAN-check before saving dsel at the end
% dataVectorLength = 33;

% Get common time-vector based on above parameters:
targetTime = getNMRTime( fullfile(optns.nmr_data_path, dset{refSet}), refExpt{refSet},dset_t0{refSet});
% v5a - disabled
% if exist('dataVectorLength','var')
%     targetTime = targetTime(1:dataVectorLength);
% end


%% HSQC
%==========================================
if include_hsqc
%%% In old-style CDS - this could have been done for several datasets at
%%% once - to get STD & etc from replicates. But in create_data_str4 -
%%% assumes only one dataset.
assert(numOfSets==1, 'Error: create_data_str4 currently can only parse one HSQC dataset at a time.')

load(hn_optns.ivtnmr_path);

hn_peak_numbers = cellfun(@(x) x(2:end), hn_optns.hn_peak_names, 'un', 0);
[~, hn_peak_idx_in_ivtnmr_struct] = ismember( hn_peak_numbers, s.hn.names);
    
% % TODO: preallocate
% % 'integr_regions',[],'err',[] - only when auto-calc of error from noise
% hsqc_raw(1,numOfSets) = struct('data',[],'time',[]);
% % in create_data_str4 'time' in the above structure is not used.
% 
% hsqc_norm = hsqc_raw;
% 
% for i=1:numOfSets
%     hsqc_raw(i).data = s.hn.csp(hn_peak_idx_in_ivtnmr_struct,:);    
%     % Scaling by protein conc and RNA size
%     % RNA - in nts. need to adj Prot conc to RNA expressed in conc of 31P atoms.
%     %%% DECIDED TO DO THIS NOW IN THE MODEL FILE?!
% %     hsqc_norm(i).data = hsqc_raw(i).data .* hn_optns.p_conc .* s.rna_length;
%     hsqc_norm(i).data = hsqc_raw(i).data;
%     hsqc_norm(i).time = s.hn.time;
% 
%     % Transpose to match old notation
%     hsqc_raw(i).data = hsqc_raw(i).data';
%     hsqc_norm(i).data = hsqc_norm(i).data';
% end


% create_data_str4: since we might cut the HSQC vector - need to have the
% targetTime for extrapolation matched to that length.
max_time_in_CSP_data = max(s.hn.time);
[~, last_31P_time_point_near_max_CSP] = min(abs(targetTime-max_time_in_CSP_data));

% Interpolate to common time vector
hsqc = struct('time',[],'mean',[],'std',[]); % Preallocate structure to hold mean & std values after interpolation.
%%% create_data_str4 - not using this here.
% v2 - expects one structure with several datasets inside.
% if d n give a numeric time vector (first argument) - will auto-interpolate for the shortest one.
% [hsqc.time, hsqc.mean, hsqc.std] = interpolateGetMeanGetSTD2(targetTime,hsqc_norm);

hsqc.time = targetTime(1:last_31P_time_point_near_max_CSP);
n_peaks = numel(hn_optns.hn_peak_names);
hsqc.mean = nan(n_peaks,numel(hsqc.time));
for i=1:n_peaks
    %%% REMOVED EXTRAPOLATION FROM HERE (initially added, but that could
    %%% give very high error from the last time-point!
    hsqc.mean(i,:) = interp1(s.hn.time, s.hn.csp(hn_peak_idx_in_ivtnmr_struct(i),:), hsqc.time);
end; clear i;

% Estimating STD as an average of the lowest quartile(or half) of pair-wise STDs
% between consequtive points. Reasoning:
% 1. Assuming there is some saturation (equilibration) at the end of the process - the
% std between two consecutive points in this region will equal the error.
% 2. Using one or two lowest quartiles - seems a conservative estimate where the changes
% are few.

fraction_for_std = 1/2;
data_for_std = hsqc.mean(:,ceil(end-end*fraction_for_std):end);
% assemble points into 3D array - putting consequtive values in 3rd dim
pair_wise_in_3D = cat(3, data_for_std(:,1:end-1), data_for_std(:,2:end));
% calculate std across 3rd dimension:
pair_wise_STD_list = std(pair_wise_in_3D,0,3);
% get MEAN (average. MODE  was too small!) STD for the selected region:
average_std_for_peaks = nanmean(pair_wise_STD_list,2);
% repeat this for whole set:
hsqc.std = repmat(average_std_for_peaks,1,numel(hsqc.time));

% Transpose to match old notation
hsqc.mean = hsqc.mean';
hsqc.std = hsqc.std';

hsqc.names = hn_optns.hn_peak_names;

% STD - as the max difference between interpolated and non-interpolated
% -- This was an old-way for calculating the STD.
% -- It does not make too much logic, cuz if time-vectors are close, STD
% will be lower. and if one of the time-points is far in the beginning of
% the trace (e.g. recorded many HSQCs, but few 31P) - the STD will be high.
% % if numOfSets == 1
% %     hsqc_std = std(cat(3, hsqc_norm(1).data, hsqc.mean), 0, 3);
% % %     % each peak own std:
% % %     maxstd = max(hsqc_std); % default - across rows
% % %     hsqc.std = repmat(maxstd, numel(targetTime));
% %     % same std for all
% %     maxstd = max(hsqc_std(:)); % default - across rows
% %     hsqc.std = repmat(maxstd, size(hsqc.mean));
% % end

%%% Plot HSQC to check
%==========================================
if plot_hsqc
subplots = numel(hsqc.names);
FIG = FIG+1; columns = subplots; rows = 2;
figure(FIG); set(figure(FIG), 'Color', [1,1,1], 'Position', [0 500 columns*sbWidth*1.4 rows*sbHeight*1.4]);        

toPlot = hsqc;

for i=1:subplots
    subplot(rows,columns,i);
%     plot(hsqc_norm(1).time,hsqc_norm(1).data(:,i), 'o-');
%     hold on;
    
    % plot interpolated
%     plot(hsqc.time,hsqc.mean(:,i), 'o-r'); 
%     plot(toPlot.time,toPlot.mean(:,i));
    errorbar(toPlot.time, toPlot.mean(:,i), toPlot.std(:,i),'r');
    axis tight;
    if i==1
        ylabel('mM (31P atoms)');
        legend(dsetname,'Location','Best');
    end;
    title(strcat(hsqc.names(i)),'FontSize',11);    
end

end

end % include hsqc


%% P31 import, normalization & scaling
%==========================================
% Vars for normalization and calc of the total RNA integrals:
%%% v1 ("OLD") (can overwrite the automatic T1 corr generation above)
T1_corr_old   = [2.8020    1.1840    1.2180    1.1520    1.2560    1.3170    1.1100    1.0660    1.0940    1.109];
%%% v2 ("NEW")
T1_corr_new   = [2.7719    1.3584    1.1576    1.7731    1.3164    1.1008    1.0702    1.9033    1.0604    0.9753];
%%% v3 (IN46c - 2016 NEW)
T1_corr_IN46c = [3.0929    1.1980    1.1661    1.1403    1.3203    1.3471    1.0509    1.1193    1.0418    1.0610];
% p31.names =          {'PO4',   'RNA',   'gNTP',   'gRNA?!',  'bNDP?',  'aNDP?',  'aNTP',   'aRNA?!',  'bNTP',  'bRNA?'};
% T1_corr_factor = mean([T1_corr_old; T1_corr_new]);
T1_corr_factor = T1_corr_IN46c;

% Probably this one is not even needed here, since free PO4 is NOT
% normalized to T1, but referenced directly to own calibration?!
PO4_free_T1_corr = 3.8917; % from 161017_20mM_PO4_TTD_integral_ref_303K_600/analysis_T1.m

RNApeaks = [2 4 8 10];

if set_RNA_T1_to_1
    T1_corr_factor(RNApeaks) = 1; % TMP. v3e_change_T1_corr_f_RNA
end
% T1_corr_factor(1:end) = 1; % TMP. v3e_change_T1_corr_f_RNA

% TODO: could n make work w/o rigid preallocation of array fields.
% p31_raw(1,numOfSets) = struct(); % this d n ra!
p31_raw(1,numOfSets) = struct('data',[],'time',[],'integr_regions',[],'err',[]); % Preallocate memory for the target array:

p31_norm = p31_raw;

for i=1:numOfSets
    % NTP_RNA_border
%     p31_raw(i) = getNMRData3(dset{i},dset_exp_31P{i},dset_t0{i},'31P',2,2,'',noise_31P); % 1st row - ref spectrum. 1st column - TopSpin spectr index. 
    p31_raw(i) = getNMRData4(dset{i},dset_exp_31P{i},dset_t0{i}, p31_integr_file, 2,2, CURDIR, noise_31P); % 1st row - ref spectrum. 1st column - TopSpin spectr index. 
        
    % Referencing / scaling.
    % PO4 is already referenced during integration. 
    % TODO: check in the 31P_datasets - the first spectrum should NOT BE
    % 5000. If it is - means the data WAS NOT REFERENCED!
    if use_external_PO4ref && ~use_average_PO4
        p31_norm(i).data(:,1) = p31_raw(i).data(:,1);
        p31_norm(i).err(:,1) = p31_raw(i).err(:,1);
        norm_start = 2;
    else
        norm_start = 1;
    end    
    
    % T1-normalize integrals against aNTP, and then scale to 20mM (starting conc of alpha PO4 in NTPs).
    refValue = p31_raw(i).data(1,refPeak); % get the integral value of aNTP peak.
    for k = norm_start:size( p31_raw(i).data, 2 ) % cycle through data
        p31_norm(i).data(:,k) = p31_raw(i).data(:,k)...
            ./ (refValue .* T1_corr_factor( k )/T1_corr_factor(refPeak))... % T1-normalize
            .* refConc; % scale to (usually) 20mM
        
        % normalize errors too!
        p31_norm(i).err(:,k) = p31_raw(i).err(:,k)...
            ./ (refValue .* T1_corr_factor( k )/T1_corr_factor(refPeak))... % T1-normalize
            .* refConc; % scale to (usually) 20mM
    end            

    if use_average_PO4
        p31_norm(i).data(:,1) = mean([p31_raw(i).data(:,1) p31_norm(i).data(:,1)],2);
        p31_norm(i).err(:,1) = mean([p31_raw(i).err(:,1) p31_norm(i).err(:,1)],2);
    end

    n_of_31P_peaks = size(p31_norm(i).data,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Extend the source array to add MgPO4 observable and RNA sum value:
    p31_norm(i).data = padarray(p31_norm(i).data, [0 2], 0, 'post');        
    p31_norm(i).err = padarray(p31_norm(i).err, [0 2], 0, 'post');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Generate MgPO4 observable:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Total integral over all 31P peaks
    tot31p = sum( p31_norm(i).data(:,1:n_of_31P_peaks) ,2);
    tot31p_err = sum( p31_norm(i).err(:,1:n_of_31P_peaks) ,2);
        
%     % Normalize total integral to total PO4, and same for error
%     % This is not needed anymore - normalization to refPeak above r e enough
%     tot31p_norm = tot31p ./ (tot31p(1,1) ./ tot_phosph_in_NTPs); 
%     tot31p_err_norm = tot31p_err ./ (tot31p(1,1) ./ tot_phosph_in_NTPs);
    % To e sure: check x above normalization y n enough:
    assert(max(tot31p) <= tot_phosph_in_NTPs, ...
        'Normzalization error: data says there is more 31P in sample than possible.');
        
    MgPO4 = tot31p(1,1) - tot31p;
    MgPO4_err = tot31p_err; % error does/shall not change here!!!
    % Append MgPO4 to the output matrix:
    p31_norm(i).data(:,end-1) = MgPO4;
    p31_norm(i).err(:,end-1) = MgPO4_err;
    clear norm_31P_T1 norm_31P_T1_err tot31p tot31p_norm tot31p_err tot31p_err_norm MgPO4 MgPO4_err;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate total RNA integral and error - and add to the data vector.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Make a sum of RNA peaks, and write it into the source array:
    p31_norm(i).data(:,end) = sum( p31_norm(i).data(:,RNApeaks) ,2);
    p31_norm(i).err(:,end) = sum( p31_norm(i).err(:,RNApeaks) ,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    p31_norm(i).time = p31_raw(i).time;
end

p31 = struct('time',[],'mean',[],'std',[]); % Preallocate structure to hold mean & std values after interpolation.

% v2 - expects one structure with several datasets inside.
% if d n give a numeric time vector (first argument) - will auto-interpolate for the shortest one.
[p31.time, p31.mean, p31.std] = interpolateGetMeanGetSTD2(targetTime,p31_norm);
p31.names = {'PO4','RNA','gNTP','gRNA?!','bNDP?','aNDP?','aNTP','aRNA?!','bNTP','bRNA?','MgPO4','RNAtotal'};

%%% Error estimate:
%%% 1 - from std of first two points in spectrum.
%%% This is OK for peaks which are low at the start and/or change slowly.
%%% But otherwise this is not ideal, since observables get diff errors
%%% (e.g. PO4 which is low at start, and aNTP which is high at start).
%%% Ideally we want to quantify only contribution of noise.
%%% 2 - uniform std
%%% 3 - noise from spectrum

p31_std = std(p31.mean(1:2,:));

switch err_calc
    case 1 % for each peak separately
        p31.std = repmat(p31_std, size(p31.mean,1), 1);
    case 2 % uniform - max of all
        p31.std = repmat(max(p31_std), size(p31.mean));
    case 3 % from noise and peak width
        if numOfSets > 1
            fprintf(1,strcat('WARNING: Error calc from noise for several dsets not implemented yet!\n\n',...
            'Using noise from the first set.\n'));
        end
        % Just using error from first dataset, assuming its the only one!
        p31.std = repmat(p31_norm(1).err, size(p31.mean,1), 1);
end

clear p31_std;

% return;


%% Plot p31 to check
%========================
if plot_p31
toPlot = p31;
% plot(toPlot.time,toPlot.mean); errorbar(toPlot.mean,toPlot.std); % one plot
% Separate plots:
subplots = size(toPlot.mean,2); % or numel(plot_names)
columns = 12;
rows = max(2, ceil( subplots/columns )); % auto-calc how many rows needed
FIG = FIG+1;
figure(FIG); set(figure(FIG), 'Color', [1,1,1], 'Position', [0 500 columns*sbWidth rows*sbHeight]);        

for i=1:subplots
    subplot(rows,columns,i);
    plot(toPlot.time,toPlot.mean(:,i)); errorbar(toPlot.time, toPlot.mean(:,i),toPlot.std(:,i),'b');
    axis tight;
    title(strcat(toPlot.names(i)),'FontSize',11);
    if i==1
        legend(dsetname,'Location','Best');
    end
    
%     if i == 11 % plot individual Mg data
%         hold on
%         plot(p31_norm(1).time,p31_norm(1).data(:,11),'r');
%         plot(p31_norm(2).time,p31_norm(2).data(:,11),'g');
%         disp(toPlot.mean(end,i));
%         hold off
%     end
    
end

clear p31_raw p31_norm;
end % plotting p31

% return;


%% Selecting specific data for observables and (if necessary) scaling to correct concentration. 
%==============================================================================================

% ===== PO4 =====  (scaled during integration in TopSpin).
PO4     = p31.mean(:,1);
PO4_std = p31.std(:,1);

% ===== NTP =====  (scaled during data import).
NTP     = p31.mean(:,7);
NTP_std = p31.std(:,7);

% ===== NDP =====  (scaled during data import).
% Before made an average of bNDP & aNDP.
% Now - using aNDP signal (as aNTP for NTPs).
NDP     = p31.mean(:,6);
NDP_std = p31.std(:,6);

% ===== MgPO4 =====  (scaled during data import).
MgP     = p31.mean(:,end-1);
MgP_std = p31.std(:,end-1);

% ===== RNA =====  (scaled during data import. uses a combination of RNA peaks as observable).
RNA     = p31.mean(:,end);
RNA_std = p31.std(:,end);

% RNA as fraction of total RNA - now moved to the main script!
RNAtotal = RNA; % was: RNA.*smn_fraction;
RNAtotal_std = RNA_std; % was: RNA_std.*smn_fraction;
Abortives = RNA; % was: RNA-RNAtotal;
Abortives_std = RNA_std; % was: RNA_std-RNAtotal_std;

% ===== Protein =====
if include_hsqc
    HSQC_peak = 1; % include first peak into the main struct, others - below separately
    PROT     = hsqc.mean(:,HSQC_peak);
    PROT_std = hsqc.std(:,HSQC_peak);
    PROT_name = sprintf('%s(%s)', s.prot, hsqc.names{HSQC_peak});
else
    PROT = zeros(size(NTP));
    PROT_std = zeros(size(NTP));
    PROT_name = 'PROT not def.';
end

%% Assembly of data structures // used to do with diff imino observables.
%======================================================================
    
dsel.sname = dsetname;
dsel.name = dset{:};
dsel.t = targetTime;

if include_hsqc
    hsqc_points = numel(hsqc.time);
    dsel.names = {'PO4', 'NTP', 'RNAtotal', 'Abortives', 'NDP', 'MgPO4', PROT_name};
    if hsqc_points == numel(targetTime)
        dsel.y = [PO4 NTP RNAtotal Abortives NDP MgP PROT];
        dsel.s = [PO4_std NTP_std RNAtotal_std Abortives_std NDP_std MgP_std PROT_std];
    else
        % make full vector w/o protein
        dsel.y = [PO4 NTP RNAtotal Abortives NDP MgP];
        dsel.s = [PO4_std NTP_std RNAtotal_std Abortives_std NDP_std MgP_std];    
        % cut it down to hsqc size, and add the protein:    
        dsel.y = [dsel.y(1:hsqc_points,:) PROT];
        dsel.s = [dsel.s(1:hsqc_points,:) PROT_std];    
    end
else
    dsel.names = {'PO4', 'NTP', 'RNAtotal', 'Abortives', 'NDP', 'MgPO4'};
    dsel.y = [PO4 NTP RNAtotal Abortives NDP MgP];
    dsel.s = [PO4_std NTP_std RNAtotal_std Abortives_std NDP_std MgP_std];    
end; % include_hsqc

% Have to transpose the DATA vectors for the main modeling script;
dsel.y = dsel.y';
dsel.s = dsel.s';

%%% Cut vectors to remove end-points with NaN (could get due to
%%% interpolation w/o EXTRApolation).
real_data = ~sum(isnan(dsel.y),1);
dsel.y = dsel.y(:,real_data);
dsel.s = dsel.s(:,real_data);
dsel.t = dsel.t(real_data);

if include_hsqc
    % Save all HSQC peaks too
    dsel.hn.note = 'CSP NOT scaled by protein conc and RNA size.';
    dsel.hn.mean = hsqc.mean(real_data,:);
    dsel.hn.std = hsqc.std(real_data,:);
    dsel.hn.names = hsqc.names;
    dsel.hn.p_conc = hn_optns.p_conc;
    dsel.hn.rna_length = s.rna_length;
end

    if export_files % Export files?
        
        if include_hsqc
            export_suffix = '_HN';
        else
            export_suffix = '';
        end
        
        if use_average_PO4
            save_mat_name = sprintf('data_nmr_selected_%s_avg_%s%s.mat', dsetname, integr_ref, export_suffix);
        else
            save_mat_name = sprintf('data_nmr_selected_%s_%s%s.mat', dsetname, integr_ref, export_suffix);
        end
        
        fprintf(1,'Data saved into %s\n', save_mat_name);
%         save(strcat('./data/',save_mat_name),'dsel'); % cds4
        save( fullfile(optns.data_for_fit, save_mat_name), 'dsel'); % cds5
        
        
    end % if "Export files?"

% Plot of results (was used before the above SWITCH option)
if plot_final
% Separate plots:
subplots = numel(dsel.names);
columns = subplots; % modify if want specific num of rows
rows = max(2, ceil( subplots/columns )); % auto-calc how many rows needed, but minimum two
FIG = FIG+1;
figure(FIG); set(figure(FIG), 'Color', [1,1,1], 'Position', [0 500 columns*sbWidth rows*sbHeight]);
for i=1:subplots
    subplot(rows,columns,i);
    plot(dsel.t,dsel.y(i,:)); errorbar(dsel.t, dsel.y(i,:),dsel.s(i,:));
    axis tight;
    title(strcat(dsel.names(i)),'FontSize',fSize);
    if i==1
        legend(dsetname,'Location','Best');
    end;
end

end % if "Plot of results"

if disp_PO4_Mg_final_conc
    fprintf(1,'== Components and sum of PO4:\n');
    tot = sum(dsel.y([1 3 6],end)) + dsel.y(2,end)*3 + dsel.y(5,end)*2; % NTP has 3 phosphates
    fprintf(1,'%.1f(PO4) + %.1f(MgPO4) + %.1f(NTP)*3 + %.1f(NDP)*2 + %.1f(RNA) = %.1f\n', dsel.y(1,end), dsel.y(6,end), dsel.y(2,end), dsel.y(5,end), dsel.y(3,end), tot);
end

end

%% Helper functions
%==============================
function time0 = get_time0(nmr_data_path, dset)
    CURDIR = nmr_data_path;

    notes_path = fullfile(CURDIR,dset,'notes.txt');

    filetext = fileread(notes_path);
    pattern = '<time0>(.*?)</time0>';
    time0 = regexp(filetext,pattern,'tokens');
    assert(numel(time0) == 1, 'Abort: more than one time0 matches found in notes.txt.');    
    time0 = time0{:}; % convert from cell to a string
    clear filetext;
end
