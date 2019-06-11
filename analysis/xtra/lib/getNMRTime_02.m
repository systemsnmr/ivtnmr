function time = getNMRTime_02(dataset_full_path,experiments,time0)

% Versions
% getNMRTime_02
% - introduces timezone-awareness
% - removes the option to just use dataset name!

%% Params for testing (are set only when function is run "locally" (when no parameters are passed in)
if nargin == 0
    datasetName_or_full_path = '160816_IN47a_p213-B3E_noRE_NUP1_303K_600';
    experiments = [5032, 2032, 6000+32*3];
    time0 = '2016-08-16 13:08:30';
    
    datasetName_or_full_path = '160816_IN47a_p213-B3E_noRE_NUP1_303K_600';
    experiments = [5032, 2032, 6000+32*3];
    time0 = '2016-08-16 13:08:30 +0200';
    
    dataset_full_path = fullfile(pwd,'dataset');
    experiments = 5000:5005;
    time0 = '2017-03-22 22:11:50 +0000';
%     time0 = '2017-03-22 22:11:50';
end

% fprintf('before\n %s \n', time0);

%% Actual script
%===========================
% Correct time0 to UTC if needed.
time0_parts = regexp(time0, ' ', 'split'); % strsplit(time0,' '); % strsplit - newer Matlab
n_parts = numel(time0_parts);
time0_field_length = numel(time0);

if n_parts==2 && time0_field_length==19
    fprintf('Time0 defined w/o timezone. \nAssuming notes.txt and audita.txt are in same timezone (e.g. measured on TopSpin < 4).\n')
    time0 = datevec(time0);
elseif n_parts==3 && numel(time0)==25 && numel(time0_parts{3})==5
    fprintf('Time0 defined with timezone. \nAssuming >=TopSpin4 and >=NEO console, and \nconverting time0 to UTC to match audita format.\n')
    time0 = datevec(time0(1:end-6));
    HH = str2double(time0_parts{3}(1:3));
    MM = str2double(strcat(time0_parts{3}(1),time0_parts{3}(4:5)));
    time0(4) = time0(4)-HH; % this should be minus if timezone has +
    time0(5) = time0(5)-MM; % this should be minus if timezone has +
    time0 = datevec(datenum(time0)); % convert to numeric and then back
else
    error('Unrecognized time0 format (expecting CHAR optionally with timezone: YYYY-MM-DD HH:MM:SS[ x0000]). Aborting..\n');
end;

n_expnos = numel(experiments);
time = nan(n_expnos,1);
for i = 1:n_expnos
    iExp = experiments(i);
    filetext = fileread( fullfile(dataset_full_path, num2str(iExp), 'audita.txt') );

%     sFormat = 'started at (?<year>\d+)-(?<month>\w+)-(?<day>\d+) (?<hour>\d+):(?<minute>\d+):(?<second>\d+).(?<ms>\d+)';
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

end