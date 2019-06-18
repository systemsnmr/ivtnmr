function full_peaklist = get_cara_peakshifts_04(cara_repo_path, target_project, target_peaklist, optns)

% > IN CURRENT IMPLEMENTATION - NEED TO HAVE ANALYZED SPECTRA IN LINEAR SEQUENCE for CARA IDs!

% @Y: seems this func has two ways of getting nodes - with or w/o XPATH
% expressions.

% 2017-07-17 - @YN:
% This is a quick-n-dirty implementation to get peak shifts from cara
% repositories. Many improvements can be done :-) E.g.
% - Use XPATH to directly select specific branches of XML tree - w/o
% needing to traverse whole thing. E.g.
%       expression = xpath.compile('//project[@name="IN70a"]/peaklist[@name="170425_IN70a"]/spec');
%       sel_nodes = expression.evaluate(xDoc, XPathConstants.NODESET);
% > see file:///Volumes/Data/yar/Dropbox/Programming/Matlab/tests/test_xml3.m

% - ...

% get_cara_peakshifts.m - used for Kd analysis
% get_cara_peakshifts_02.m - modified version, to not depend on TIME0 /
% IVTNMR stuff
% get_cara_peakshifts_03 - adds reading of peak AMPLITUDE!
% get_cara_peakshifts_04 - uses optns structure to get info about nmr data
% path (with SPECTRA) -- not just read it from cara file)


%% If called locally
%===================================
DEBUG = 1;

if nargin == 0
    fprintf(1,'%s requires at least X input arguments. Running test example.\n\n', mfilename);    
%     file_name = '/Volumes/Data/yar/Dropbox/Programming/Matlab/testdata/170425_UP1_A1R1_R0_trace.cara'; 
%     target_project = 'IN70a';
%     target_peaklist = '170425_IN70a';
    file_name = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/170915_5Pz_PO4_Mg_NTP.cara';
    target_project = '5Pz';
    target_peaklist = 'hnco';
    
    file_name = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra/190415_IN115a_pA20_co-P50N525_100uM_298K_600/IN115a_v2_trace.cara';
    target_project = 'v1';
    target_peaklist = 'IN115a';
    optns.nmr_data_path = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra/190415_IN115a_pA20_co-P50N525_100uM_298K_600';
else
    file_name = cara_repo_path;
end

% Test variables for existence
% if ~exist('export_files','var') || isempty(export_files)
%     export_files = 0;
% end

%% Code
%===================================
xDoc = xmlread(file_name);

import javax.xml.xpath.*;
factory = XPathFactory.newInstance;
xpath = factory.newXPath;

% 2. Find all projects
% The getElementsByTagName method returns a deep list that contains information about the child nodes:
all_proj = xDoc.getElementsByTagName('project');
n_proj = all_proj.getLength;
fprintf(1,'\n== CARA file : %i projects\n', n_proj);

% % attr = proj.item(0).getAttributes();
% attr = char(all_proj.item(0).getAttribute('name'));
% all_attr = all_proj.item(0).getAttributes;

for k = 0:n_proj-1 % because this is JAVA object - indexing goes from ZERO
    proj = all_proj.item(k);
    proj_name = proj.getAttribute('name'); % traverse names
    if strcmp(proj_name, target_project)
            
        % Get all peaklists
        all_peaklists = proj.getElementsByTagName('peaklist');
        n_peaklists = all_peaklists.getLength;
        fprintf(1,'\n  == Project %s : %i peaklists\n', char(proj_name), n_peaklists);

        for L = 0:n_peaklists-1
            peaklist = all_peaklists.item(L);
            peaklist_name = peaklist.getAttribute('name'); % traverse names
            if strcmp(peaklist_name, target_peaklist)
                
                peaks = peaklist.getElementsByTagName('peak');
                n_peaks = peaks.getLength;

                spec = peaklist.getElementsByTagName('spec');
                n_spec = spec.getLength;
                                                
                fprintf(1,'\n  == Peaklist %s : %i peaks, %i spec\n', char(peaklist_name), n_peaks, n_spec);
                
                fl.H = nan(n_peaks, n_spec);
                fl.N = nan(n_peaks, n_spec);
                fl.amp = nan(n_peaks, n_spec);
                fl.peak_tags = cell(n_peaks,1);
                                
                last_spec_id = str2num(spec.item(n_spec-1).getAttribute('id'));
%                 fl.spec_ids = str2num(spec.item(0).getAttribute('id')) : last_spec_id;                
                
                fl.spec_ids2 = nan(1, n_spec);
                
                for iSpec = 0 : (n_spec-1)
                    fl.spec_ids2(iSpec+1) = str2num(spec.item(iSpec).getAttribute('id'));
                end
                
                fl.spec_ids = fl.spec_ids2;
                                                
                if numel(fl.spec_ids) ~= n_spec
                    fprintf(1,'=====================================================\n');
                    fprintf(1,'WARNING - spectra may have non-linear numbering (ids).\n');
                    fprintf(1,'=====================================================.\n');
                end
                
                %%% To get the TIME vector.
                % Search for the spectrum with ID matching ID of last spectrum 
                % in the peaklist.
                spectrum = proj.getElementsByTagName('spectrum'); % get all spectrum nodes
                n_spectrum = spectrum.getLength;
                
                item_idx = n_spectrum-1; % search from the end of the list
                while str2num(spectrum.item(item_idx).getAttribute('id')) ~= last_spec_id
                    item_idx = item_idx-1; % search from the end of the list
                end
                
                % Get data path of the matched spectrum. Extract dset name
                % and expno numbers.
                fl.last_path = char(spectrum.item(item_idx).getAttribute('path'));
                path_parts = regexp(fl.last_path, filesep, 'split');
                
%                 fl.nmr_data_path = fullfile(filesep,path_parts{1:end-5});
                fl.nmr_data_path = optns.nmr_data_path;
                
                fl.dset_name = path_parts{end-4};
                fl.last_expno = str2num(path_parts{end-3});
                try 
                    fl.time0 = getTime0(fl.nmr_data_path,fl.dset_name);
                    % In calculating first expno - we assume expnos go
                    % linearly.                             
                    fl.time = getNMRTime( fullfile(optns.nmr_data_path, fl.dset_name),...
                        [fl.last_expno-n_spec+1:fl.last_expno], fl.time0);
                catch
                    warning('Some error in processing time0 in notes.txt. Using 1:n_spec index instead of TIME.');
                    fl.time = 1:n_spec;
                end                

                % Get shifted positions for each peak
                fprintf(1,'Parsing peaklist ...\n');
                for P = 0:n_peaks-1
                    peak = peaks.item(P);
                    peak_id = peak.getAttribute('id');
                    peak_tag = peak.getAttribute('tag');
                    fl.peak_tags{P+1} = char(peak_tag); %%%%% FILL MATRIX %%%%%
                    
                    if DEBUG && nargin == 0
                        fprintf(1,'\t\tid=%s tag=%s\n', char(peak_id), char(fl.peak_tags(P+1)));
                    end

                    % Get all POS elements for specific peak:
%                     if k == n_proj-1 && L == n_peaklists-1 && P == n_peaks-1
                        positions = peak.getElementsByTagName('pos');
                        n_positions = positions.getLength;

                        format long;
                        for pos_idx = 0:n_positions-1
                            position = positions.item(pos_idx);
                            pos_spec = position.getAttribute('spec');
                            pos_amp = position.getAttribute('amp');
                            spec_idx = int16(find(eq(fl.spec_ids, str2num(pos_spec))));
                            
                            if DEBUG && nargin == 0
                                fprintf(1,'spec_idx = %i\n', spec_idx);
                            end

                            dim = position.getElementsByTagName('dim');
                            n_dim = dim.getLength;
%                             H_shift = dim.item(0).getAttribute('pos');
%                             N_shift = dim.item(1).getAttribute('pos');

                            %%%%% FILL THE MATRIX %%%%%
                            % There are two non-normal cases here:
                            % 1. NO <POS> element at all - then spec_idx is
                            % absent
                            % 2. NO <DIM> sub-elements in <POS> element
                            % (for this one there is if condition below - to avoid errors in loop)
                            
                            % Both should be treated to set peak position
                            % to value of reference peak (first spectrum).
                            % Done below.
                            
                            if n_dim == 2
                                fl.H(P+1,spec_idx) = str2double(dim.item(0).getAttribute('pos'));  
                                fl.N(P+1,spec_idx) = str2double(dim.item(1).getAttribute('pos'));
                            elseif n_dim == 1
                                fprintf(1,'WARNING - peak has only one alias (H or N) - need to adjust the code for this!.\n');
                            end

%                             fprintf(1,'\t\t\tspec=%s H=%s N=%s\n', char(pos_spec), char(fl.H(P,spec_idx)), char(fl.N(P,spec_idx)));

                            % get_cara_peakshifts_03
                            fl.amp(P+1,spec_idx) = str2double(pos_amp);
                        end

                end
            end % select peaklist

        end % loop over peaklists
    end  % select project
end % loop over projects

%% Fill the NAN parts of matrix - with values from reference spectrum
%======================================================================
% Matrix with 1st columns (ref spectra) instead of NaNs:
H_nan_replace = bsxfun(@times, isnan(fl.H), fl.H(:,1));
N_nan_replace = bsxfun(@times, isnan(fl.N), fl.N(:,1));
% Convert NaN to zeros (cuz NaN+1 = NaN):
fl.H(isnan(fl.H)) = 0;
fl.N(isnan(fl.N)) = 0;
% add to zeros the values from 
fl.H = fl.H + H_nan_replace;
fl.N = fl.N + N_nan_replace;


%% Calc CSP
%======================================================================
fl.H_csp = bsxfun(@minus, fl.H, fl.H(:,1));
fl.N_csp = bsxfun(@minus, fl.N, fl.N(:,1));

fl.HN_csp = sqrt( ( fl.H_csp.^2+(fl.N_csp./5).^2 )./2 );

full_peaklist = fl;

end % main function

