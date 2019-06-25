function [xaxis, traj_select] = nans_drop(xaxis, traj_select)
% Removes NaN-data for the Kd fits. Assumes peaks in rows in traj_select.

%% Params for testing (when function is run "locally")
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    xaxis = {[1:10]'; [2:20]'};
    traj_select = {[1:10; 1:10]; [2:20; 2:20]};
    
    traj_select{2}(1,8) = NaN;
end

%% Actual script

n_sets = numel(traj_select);
set_ids = 1:n_sets;
nan_idx = cellfun(@(x) logical(sum(isnan(x), 1)), traj_select, 'un', 0);
arrays_with_nans = logical(cell2mat(cellfun(@(x) sum(x), nan_idx, 'un', 0)));

% xaxis
% traj_select

if any(arrays_with_nans)
    for iSet=set_ids(arrays_with_nans)
        traj_select{iSet}(:,nan_idx{iSet}) = [];
        xaxis{iSet}(nan_idx{iSet}) = [];
    end
end;

% xaxis
% traj_select

end