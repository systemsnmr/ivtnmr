function d = getNMRData4(datasetName,experiments,time0,dataType,startingRow,startingColumn,dataPath,noise_region)

%%% 2015-08-12. @YAR.
%   - Added dataPath - optional parameter to indicate the path to the folder with data ("bruker").
%   - Preallocated array for "time".

%%% 2016-09-28. @YAR. - getNMRData3
%   - Added read-out and passing back "integr_regions"
%   - Added integration error estimation - using noise-region

%%% 2017-06-16. @YAR. - getNMRData4
%   - Fixed noise-based integration error estimation (added power2 to
%   rmse_norm)
%   - Added ';' delimeter in the importdata statement

%% Help / annotation
% d = getNMRData3(datasetName,experiments,time0,dataType[,sRow,sCol,dataPath,noise_region])
% Example call:
% getNMRData3('130520_IN17_pl13_TB_31P_600',[500:526],'2013-05-20 18:42:10','PO4')
% 
% Retrieves time and integral values from NMR data.
% 
% Input arguments:
% 
% datasetName = name of NMR directory with data.
% experiments = a range of experiments representing the data. 
% time0 = data and time at which experiment was started (for calculation of
% time vector).
% dataType = type of data: [PO4, 31P, imino, imino_fast, 1H, tocsy, hsqc, hsqc_CS].
%            PO4, 31P, 1H, imino, imino_fast are expected to be exported from TopSpin.
%            tocsy, hsqc, hsqc_CS - batch integration file from Cara.
% sRow = number of rows to skip in the output (e.g. if integration used a
% reference spectrum which is not in the dataset.
% sCol = number of columns to skip (definitly for Bruker spectrum ID)

% dataPath = if data is not in /Volumes/Data/yar/....

% noise_region = a region w/o pekas to be used for error estimation


%% Error checks: TODO (add more of these)
%% Params for testing (are set only when function is run "locally" (when no parameters are passed in)
if nargin == 0
    % datasetName = '130520_IN17_pl13_TB_31P_600';
    % experiments = [500:526];
    % time0 = datevec('2013-05-20 18:42:10','yyyy-mm-dd HH:MM:SS');
    % dataType = 'tocsy';

    % datasetName = '130724_IN26_pl13_TB_RE_600';
    % time0 = '2013-07-24 18:54:59';
    % experiments = [2000:2058];
    % dataType = 'imino';
    % startingRow = 1;
    % startingColumn = 2;

    datasetName = '170308_IN60a_pR02_co-A1R1_303K_600';
    time0 = '2017-03-08 17:34:30';
    experiments = [5000:5034];
    dataType = 'integr_results_31P_pure_PO4.txt';
    startingRow = 2; % remove reference row!
    startingColumn = 2; % remove spectra index
    noise_31P = [-24,-28];
    noise_region = noise_31P;
end

% Checking dataPath variable
if nargin>6 && ~isempty(dataPath)
    if ~ischar(dataPath)
        error('dataPath parameter should be a String');
    else
        datasetFolder = strcat(fullfile(dataPath,datasetName),'/');
    end
else
    datasetFolder = strcat('/Volumes/Data/yar/_eth2/data_NMR/spectra/',datasetName,'/');
end

% Checking noise range - done at the end

%% Actual function

time0 = datevec(time0,'yyyy-mm-dd HH:MM:SS');

%% Get time
numOfExp = numel(experiments);
time = nan(numOfExp,1);
for i=1:numOfExp
    exp = experiments(i);
    filetext = fileread( strcat(datasetFolder, num2str(exp), '/audita.txt') );

    sFormat = 'started at (?<year>\d+)-(?<month>\w+)-(?<day>\d+) (?<hour>\d+):(?<minute>\d+):(?<second>\d+).(?<ms>\d+)';
    sArray = regexp(filetext,sFormat,'names');
    sString = strcat(sArray.year,'-',sArray.month,'-',sArray.day,{' '},sArray.hour,':',sArray.minute,':',sArray.second,'.',sArray.ms);
    timeStart = datevec(sString,'yyyy-mm-dd HH:MM:SS.FFF');
%     disp(sVector);
    eFormat = '(   1,<(?<year>\d+)-(?<month>\w+)-(?<day>\d+) (?<hour>\d+):(?<minute>\d+):(?<second>\d+).(?<ms>\d+)';
    eArray = regexp(filetext,eFormat,'names');
    eString = strcat(eArray.year,'-',eArray.month,'-',eArray.day,{' '},eArray.hour,':',eArray.minute,':',eArray.second,'.',eArray.ms);
    timeEnd = datevec(eString,'yyyy-mm-dd HH:MM:SS.FFF');
%     disp(eVector);

    t = (etime(timeStart,time0) + etime(timeEnd,timeStart)/2)/60; % (time from start of series + half exp time)/normalize to minutes
    if t < 0
        t = 0;
    end
    time(i) = t;
end

% Get integrals
%================================
switch dataType
    case 'PO4'
        integralFile = 'integr_results_31P_PO4.txt';
    case '31P'
        integralFile = 'integr_results_31P.txt';
    case 'imino'
        integralFile = 'integr_results_1H_imino.txt';
    case 'imino_fast'
        integralFile = 'integr_results_1H_imino_fast.txt';
    case 'imino_WESF_3pt'
        integralFile = 'integr_results_1H_imino_WESF_3pt.txt';
    case '1H'
        integralFile = 'integr_results_1H.txt';
    case 'tocsy'
        integralFile = 'integr_results_tocsy.txt';
    case 'hsqc'
        integralFile = 'integr_results_hsqc.txt';
    case 'hsqc_CS'
        integralFile = 'chemshift_results_hsqc.txt';
    otherwise
        % User can provide the name of the file instead
        integralFile = dataType;
end

dataFile = strcat(datasetFolder,integralFile);
if ~exist(dataFile,'file')
    error('Unexpected dataType. Allowed values are: PO4, 31P, imino, imino_fast, imino_WESF_3pt, 1H, tocsy, hsqc_CS(chemical shifts), hsqc(integrals), or need to specify an exact "filename.txt"');
else
    raw = importdata(strcat(datasetFolder,integralFile),';');
end
    
% Omit rows/columns as indicated in the input
data = raw.data(startingRow:end,startingColumn:end); 
% data = raw.data(:,:); 
% plot(time,data,'o');


% Get integration regions (added on 2016-09-28 @YN)
%====================================================
% integr_regions = raw.textdata; 
filepath = fullfile(datasetFolder,integralFile);
header_lines = 10;
fileID = fopen(filepath, 'r');
integr_regions = textscan(fileID, '#   %f  %f  %*f  %*f  # for region %*u', 'HeaderLines', header_lines);
integr_regions = cell2mat(integr_regions);
fclose(fileID);


% Calculate integration error for each peak - from noise (added on 2016-09-28 @YN)
% TODO - perhaps move to a separate function!
%====================================================
E = [];
if (nargin>7 || nargin == 0) && exist('noise_region','var') && ~isempty(noise_region)
    %%% Check of input
    if ~isnumeric(noise_region) || numel(noise_region) ~= 2
        error('noise_region parameter is expected to pass two numbers');
    else
        [noise_l, noise_r] = deal(noise_region(1), noise_region(2));
    end
    
    %%% Actual code
    % Helper: gives indexes of points for ppm range in rbnmr data structure
    f_find_ids = @(l,r,rbnmr_set)  find(rbnmr_set.XAxis < l & rbnmr_set.XAxis > r);

    % specify and read first spectrum
    ref_expt = experiments(1);
    spec = fullfile(datasetFolder,num2str(ref_expt),'pdata/1/1r');
    s = rbnmr(spec);
    
    % calc rmse in the given range
    ints = s.Data( f_find_ids(noise_l, noise_r, s) );
    rmse = sqrt(mean((ints-mean(ints)).^2)); % Root Mean Squared Error
    
    % get normalized and non-normalized integral of the first peak
    first_region = integr_regions(1,:);
    first_region_intensities = s.Data( f_find_ids(first_region(1), first_region(2), s) );
%     disp('not sure trapezoid is correct here (@Y assumed this is what TopSpin uses');
%     disp('but @FD says its only for special cases - e.g. smth moving.');
    first_region_integral = trapz(first_region_intensities); % Read notes above!
    first_region_integral_norm = data(1,1);
    rmse_norm = rmse * first_region_integral_norm/first_region_integral;

    % Number of points in each peak
    % diff between left and right edge of the peak, scaled by PTpPPM
    N = abs(diff(integr_regions,1,2)) ./ ((s.Procs.SW_p / s.Procs.SF) / s.Procs.SI);
    % Error of integration for each peak.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% In below POWER-of-2 WAS MISSING in getNMRData3. < E = sqrt( rmse_norm^2 .* N ) >
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = sqrt( rmse_norm.^2 .* N )'; % Transposed to match orientation of data vector

end

% Export
%====================================================
d = struct;
d.data = data;
d.time = time;
if nargin>7 && exist('noise_region','var') && ~isempty(noise_region)
    d.integr_regions = integr_regions;
    d.err = E;
end


% Visualize if just testing
%====================================================
if nargin == 0
    rows = 2;
    columns = size(d.data,2);
    n_points = numel(d.time);
    
    sbSide = 150; fSize = 11;
    figure(1);
    set(figure(1), 'Color', repmat(0.95,1,3), 'Position', [0 300 columns*sbSide rows*sbSide]);
    global SB; SB=0;
    sb = @(x) subplot(rows,columns,x,'FontSize',fSize);
    for peak=1:columns
        SB=SB+1; sb(SB);
        plot(d.time,d.data(:,peak));
        if exist('noise_region','var') && ~isempty(noise_region)
            errorbar(d.time, d.data(:,peak), repmat(E(peak),n_points,1) ,'b')
        end
        axis tight;
    end
end

end % whole function

