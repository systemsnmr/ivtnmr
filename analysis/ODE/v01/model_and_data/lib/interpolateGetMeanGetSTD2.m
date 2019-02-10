function [target_time, dataMean, dataSTD] = interpolateGetMeanGetSTD2(time,datasets)

% interpolateGetMeanGetSTD2 (version2)
% instead of separately defined datasets, requires as an input one
% structure with several datasets inside.
% Initially implemented here: /Volumes/Data/yar/Dropbox/_eth2/data_MathModeling/140721_v5m_copy_IN36_IN37_v1/create_data_str.m

%% Get Mean and STD
numOfdatasets = length(datasets);

target_time = [];
% If time vector for interpolation is not provided - interpolate to the
% shortest time vector

if ~isnumeric(time)
    target_time = datasets(1).time; % set time to time of first dataset
    for i=2:numOfdatasets
        if length(datasets(i).time) < length(target_time)
            target_time = datasets(i).time;
        end
    end
else
    target_time = time;
end

for i=1:numOfdatasets
    % Interpolate and add to dataArray
    dataArray(:,:,i) = interp1(datasets(i).time,datasets(i).data,target_time);    
end

dataMean = mean(dataArray,3); % get mean across 3rd dimension
dataSTD = std(dataArray,0,3); % get mean across 3rd dimension

% plot(target_time,dataMean,'om');
% errorbar(target_time,dataMean,dataSTD,'om');

% % If want to return data structure, use this:
% d = struct;
% d.mean = dataMean;
% d.time = target_time;
% d.std = dataSTD;

end