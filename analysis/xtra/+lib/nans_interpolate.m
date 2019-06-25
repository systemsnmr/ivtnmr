function [xaxis, traj_select] = nans_interpolate(xaxis, traj_select)
% Interpolates NaN-data in the Kd fits. Assumes peaks in rows in traj_select.

%% Params for testing (when function is run "locally")
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    xaxis = {[1:10]'; [2:20]'};
    traj_select = {[1:10; 1:10]; [2:20; 2:20]};
    
    traj_select{2}(1,8) = NaN;
    traj_select{1}(1,4) = NaN;
end

%% Actual script

n_sets = numel(traj_select);
set_ids = 1:n_sets;
peak_ids = cellfun(@(x) 1:size(x,1), traj_select, 'un', 0);

nan_idx = cellfun(@(x) isnan(x), traj_select, 'un', 0);
peaks_with_nans = cellfun(@(x) logical(sum(x,2)), nan_idx, 'un', 0);
arrays_with_nans = logical(cell2mat(cellfun(@(x) sum(x), peaks_with_nans, 'un', 0)));

% traj_select{:}

if any(arrays_with_nans)
    warning('== Interpolating some missing CSP data (NaNs) for Kd fits.');
    for iSet=set_ids(arrays_with_nans)
        for iPeak=peak_ids{iSet}(peaks_with_nans{iSet})
            % Some inspiration https://ch.mathworks.com/matlabcentral/answers/34346-interpolating-nan-s
            y = traj_select{iSet}(iPeak,:);
            nan_y = nan_idx{iSet}(iPeak,:);
            t = xaxis{iSet};
            y(nan_y) = interp1(t(~nan_y), y(~nan_y), t(nan_y), 'linear', 'extrap');
            traj_select{iSet}(iPeak,:) = y;
        end;        
%         traj_select{iSet}(:,nan_idx{iSet}) = [];
%         xaxis{iSet}(nan_idx{iSet}) = [];
    end;
end;

% fprintf(1,'\n\n\n== After fix');
% traj_select{:}

end