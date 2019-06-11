function integr_results = integrate1D(nmr_data,intrng)
% Integrate NMR spectra in specified regions.
% v001 @YN 2017-08-08
%   - input is automatically parsed with rbnmr (expno or procno)

% TODO
% - 

%% Params for testing (when function is run "locally")
if nargin == 0
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);
    intrng = getIntegrRegions('/Users/yar/Dropbox/_eth2/project_Primase/bruker/140128_PRM2_323K_600/2000');
    nmr_data = '/Users/yar/Dropbox/_eth2/project_Primase/bruker/140128_PRM2_323K_600/2000/pdata/1/1r';
    nmr_data = '/Users/yar/Dropbox/_eth2/project_Primase/bruker/140128_PRM2_323K_600/2000';
    nmr_data = {...
        '/Users/yar/Dropbox/_eth2/project_Primase/bruker/140128_PRM2_323K_600/2000'
        '/Users/yar/Dropbox/_eth2/project_Primase/bruker/140128_PRM2_323K_600/2001'
        };
end

%% Actual script
%=========================
% Helper: gives indexes of points for ppm range in rbnmr data structure
f_find_ids = @(l,r,rbnmr_set)  find(rbnmr_set.XAxis < l & rbnmr_set.XAxis > r);

% read nmr spectra
% ref_expt = experiments(1);
% spec = fullfile(datasetFolder,num2str(ref_expt),'pdata/1/1r');
s = rbnmr(nmr_data);

n_spec = numel(s);
n_int = size(intrng,1);

% Check that all spectra have same size
if size( unique( cell2mat( cellfun(@(x) size(x.Data), s, 'unif', 0) ), 'rows' ), 1) ~= 1
    error('NMR data to be integrated must be same size.');
end

% % from project_Metab_CCM/code/160705_PeakPicking/peakDetectionRef.m
% % calculate piece-wise numerical integral, step size 50 data points
% y_int = [];
% x_int = [];
% for i=1:50:length(x)-50
%    y_val=y(i:i+49);
%    y_int(end+1,1)=trapz(y_val);
%    x_int(end+1,1)=x(i+25);
% end

integr_results = nan(n_spec,n_int);

% TODO: efficiency should be improveable! (i.e. w/o for-loops)? E.g. on
% combined matrix instead of cell array:
% intensities = cell2mat( cellfun(@(x) x.Data', s, 'unif', 0) );
% trapz(intensities(:,start:finish),dim);

for i=1:n_spec
    for j=1:n_int % could do cellfun?
        %%% TODO: better make min(intrng(j,:)) and max(intrng(j,:)) -- so
        %%% that can specify ppm bounds in any orientation (increasing or
        %%% decreasing). In the below notation - have to use NMR-style
        %%% orientation - lower ppm on the right.
        y_val = s{i}.Data( f_find_ids(intrng(j,1), intrng(j,2), s{i}) );
        integr_results(i,j) = trapz( y_val );
    end
end

end % main func