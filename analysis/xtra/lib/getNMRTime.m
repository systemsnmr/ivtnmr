function time = getNMRTime(datasetName_or_full_path,experiments,time0)

%% Params for testing (are set only when function is run "locally" (when no parameters are passed in)
if nargin == 0
    % datasetName = '130627_IN22_pl23_TB_600';
    % experiments = [500:528];
    % time0 = datevec('2013-06-28 00:12:00','yyyy-mm-dd HH:MM:SS');

    % datasetName = '130724_IN26_pl13_TB_RE_600';
    % experiments = [2000:2058];
    % time0 = '2013-07-24 18:54:59';

    % datasetName = '130902_tdh_FBP_T2_LOGSY_STD_750';
    % experiments = [312 325 330:347 352];
    % time0 = '2013-09-03 15:23:32';

    % datasetName = '130903_tdh_G6P_T2_LOGSY_STD_750';
    % experiments = [312 325 330:335];
    % time0 = '2013-09-04 00:50:40';

%     datasetName = '130815_IN32_NTP_degr_TBRE_600';
%     experiments = [5000:5050];
%     time0 = '2013-08-15 23:21:00';
    
%     datasetName = '160412_Metab_AroG_catalysis_600';
%     experiments = [120:140];
%     time0 = '2016-04-12 13:47:35';        

    datasetName_or_full_path = '160816_IN47a_p213-B3E_noRE_NUP1_303K_600';
    experiments = [5032, 2032, 6000+32*3];
    time0 = '2016-08-16 13:08:30';

end


%% Actual script
time0 = datevec(time0);
%%% Hack added later (2017-08-15) to allow non-hard-fixed nmr_data_path (by
%%% passing full path in the input)
if strcmp(datasetName_or_full_path(1),'/')
    datasetFolder = strcat(datasetName_or_full_path,'/');
else
    datasetFolder = strcat('/Volumes/Data/yar/_eth2/data_NMR/spectra/',datasetName_or_full_path,'/');
end

time = [];
% var(:)' forces row vector (var(:) - column)
for i = experiments(:)'
    filetext = fileread( strcat(datasetFolder, num2str(i), '/audita.txt') );

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
    time = [time; t];
end

% disp(time);

end