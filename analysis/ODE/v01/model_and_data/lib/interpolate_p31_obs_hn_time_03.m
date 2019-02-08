function [obs_out, err_out] = interpolate_p31_obs_hn_time_03(nmr_data_path, dset_names, optns)

global FIG;
global SB;

% Interpolates P31 observables to HN time vector.
% 31P data is smoother than HN CSP - more accurate for interpolation.

% _03 - 20190204 @YN:
%   - updates local ivtnmr by default, but has a toggle to update global or
%   not.

%% Params to run "locally"
%===========================
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    nmr_data_path = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra';
    dset_names = {...
'170308_IN60a_pR02_co-A1R1_303K_600'
'170503_IN61a_SMN1_co-A1R1_303K_600'
'170309_IN62a_p214_co-A1R1_303K_600'
'170607_IN65b_EV2_co-A1R1_303K_600'
'170310_IN70a_pR02_co-UP1_303K_600'
'170313_IN71a_p114_co-UP1_303K_600'
'170426_IN72a_SMN2_co-UP1_303K_600'
'170515_IN75b_EV2_co-UP1_303K_600'
        };
        
    optns.save_to_global_ivtnmr = 0;
    optns.plot_results = 0;
end;

% Generate dataset ID, assuming its the second position in the dset_name
dset_id_names = cellfun(@(x) regexp(x,'_','split'), dset_names, 'un', 0);
dset_id_names = cellfun(@(x) x{2}, dset_id_names, 'un', 0);

global DATA_IVTNMR;

%% Actual script
%===========================
n_sets = numel(dset_names);

obs_out = cell(n_sets,1);
err_out = cell(n_sets,1);

if optns.plot_results
    rows = 2;
    columns = n_sets;
    sbSide = 200;
    fSize = 11;
    
    % global defined at the top of the function
    if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
    FIG=FIG+1;
    
    figure(FIG); set(figure(FIG), 'Color', repmat(0.95,1,3), 'Position', [0 300 columns*sbSide rows*sbSide]);
	SB=0; % global defined at the top of the function
    sb = @(x) subplot(rows,columns,x,'FontSize',fSize);
end

clear i;
for i=1:n_sets
%     fprintf('== i=%i\n', i);
    % this was reading from global:
%     data_file = fullfile(nmr_data_path,dset_names{i},'ivtnmr.mat'); %

    data_file = fullfile(DATA_IVTNMR, sprintf('ivtnmr_%s.mat', dset_id_names{i}));
    load(data_file);
    
    obs = s.p31.obs;
    err = s.p31.obs_std;
    time = s.p31.time;
    hn_time = s.hn.time;
    % interp1 - can also extrapolate - either by default if using
    % 'spline', 'pchip', 'cubic'. Or if forcing by 'extrap' option.
    
    %%%%% Its DANGEROUS TO BLINDLY use INTERPOLATION here:
    %%%%% in case if HN vector >> longer than 31P one!
    hn_end = hn_time(end);
    t_end = time(end);
    difference_threshold = 0.02; % 0.05 - 5%
    if (hn_end-t_end)/hn_end > difference_threshold
        error('%s - Aborting interpolation: HN time >> P31 time.', dset_names{i});
    end
        
%     (hn_time(end)-time(end))/hn_time(end)
    % Trasponsing to colum so can operate on a matrix column-wise
    obs_hn_time = interp1(time', obs', hn_time')';
    err_hn_time = interp1(time', err', hn_time')';
    obs_out{i} = obs_hn_time;
    err_out{i} = err_hn_time;
    
    % Update the structure
    s.p31.obs_hn_time = obs_hn_time;
    s.p31.obs_std_hn_time = err_hn_time;
    
    if optns.plot_results
        SB=SB+1; sb(SB);
        plot(time, obs,'-');
        hold on;
        plot(hn_time, obs_hn_time, 'o');
        axis tight;
    
        sb(SB+n_sets);
        plot(s.hn.time, s.p31.obs_hn_time);
        axis tight;
    end    
        
    %%% Save
    save( data_file, 's'); % updates local copy    
    if optns.save_to_global_ivtnmr
        save( fullfile(nmr_data_path, dset_names{i}, 'ivtnmr.mat') ,'s');
    end

end


end
