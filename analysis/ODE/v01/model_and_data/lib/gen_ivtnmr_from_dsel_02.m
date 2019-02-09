function gen_ivtnmr_from_dsel_02(dsel_path,dsel_suffixes,init_NTP,nmr_data_path,dset_names,plot_results,save_to_global_ivtnmr,rna_lengths)
% gen_ivtnmr_from_dsel(dsel_file,nmr_data_path,dataset_name)
% This is a converter script, but shall be easier to generate same directly from cds.m !

% + read data
% + normalize by RNA length
% + save into external structure:
%   name, sname, rna_length, rna_conc, RNA_TIME_VECTOR, rna_conc_variant_names
% + save into each dataset directory
% + read protein name from dsetname
% + read co/post from dsetname

% TODO - in a separate file:
% - import CSPs (from collect_CSP_data)
% - import 31P (from this script)
% - fit Kd - for each CSP+31P pair

global FIG;
global SB; 

global DATA_IVTNMR;

if exist('./lib','dir')
    datasave_folder = strcat(DATA_IVTNMR, '/'); % if running in root (runall)
else
    datasave_folder = strcat('../', DATA_IVTNMR, '/'); % if running from lib directory
end

%% Variables
%==================
if nargin == 0 
	
    % Just to make sure we're in the right directory (was complaining)
    [fld,~,~] = fileparts( mfilename('fullpath') );
    cd(fld);

    nmr_data_path = '/Volumes/Data/yar/Dropbox/_eth2/data_NMR/spectra';

    dset_names = {...
    '170308_IN60a_pR02_co-A1R1_303K_600'
    '170310_IN70a_pR02_co-UP1_303K_600'
    '170503_IN61a_SMN1_co-A1R1_303K_600'
    '170313_IN71a_p114_co-UP1_303K_600'
    % '170427_IN71i_SMN1_post-UP1_303K_600'
    '170309_IN62a_p214_co-A1R1_303K_600'
    '170426_IN72a_SMN2_co-UP1_303K_600'
    % '170404_IN72i_SMN2_post-UP1_303K_600'
    '170607_IN65b_EV2_co-A1R1_303K_600'
    '170515_IN75b_EV2_co-UP1_303K_600'
    % '170516_IN75k_EV2_post-UP1_303K_600'
    % '170314_IN60q_pR02_free_303K_600'
    % '170320_IN63a_pRHV3b_co-A1R1_303K_600'
    % '170608_IN65i_EV2_post-A1R1_303K_600'
    % '170627_IN73q_HV3_free_303K_600'
    };

    dsel_path = '/Volumes/Data/yar/Dropbox/_eth2/data_MathModeling/170614_IN60_IN70/v02_kcat_NTP_many_sets/data';
    init_NTP = 20; % mM
    dsel_suffixes = {'avg_PNM','avg_pure','PNM','pure'};
    plot_results = 1;
    
%     datasave_folder = '../datasave/';
    
    save_to_global_ivtnmr = 0;
    
    %%% NOW IS PASSED AS PARAMETER!! HERE ONLY LOCAL TEST!
    rna_lengths = containers.Map({...
        'IN60','IN61','IN62','IN63','IN65','IN70','IN71','IN72','IN73','IN75'},...
        [5.8, 28, 28, 33, 49, 5.8, 28, 28, 33, 49]);
    
else
    
%     datasave_folder = 'datasave/';
end

%%% Independent of local or global run
n_sets = numel(dset_names);


%% Code
%===================
n_subs = numel(dsel_suffixes); % subsets for each dataset
% Preventively open an invis figure, cuz colormap would make if none is open.
figure('Visible','off'); 
colors = colormap(hsv(n_sets));

name_split = cellfun(@(x) regexp(x, '_', 'split'), dset_names, 'unif', 0);
dset_snames = cellfun(@(x) x{2}, name_split, 'unif', 0); % extract short name
rna_names = cellfun(@(x) x{3}, name_split, 'unif', 0); % extract short name
rna_names = regexprep(rna_names, {'pR02', 'p114', 'p214'}, {'RNA0', 'SMN1', 'SMN2'}); % replace some names to more standard

cotr_n_prot = cellfun(@(x) regexp(x, '[\w+](co|post)-([^_]+)', 'tokens', 'once'), dset_names, 'unif', 0);
cotr = cellfun(@(x) x{1}, cotr_n_prot, 'unif', 0);
prot_names = cellfun(@(x) x{2}, cotr_n_prot, 'unif', 0);

for i=1:n_sets    
    fprintf(1,'== Importing and converting data from %s...\n', dset_snames{i});
        
    %%% Import 31P data and calc average of them
    %==============================================
    dsel_files = cellfun(@(x) fullfile(dsel_path,...
        sprintf('data_nmr_selected_%s_%s.mat', dset_snames{i}, x)), dsel_suffixes, 'unif', 0); % gen paths    
    ds = cellfun(@(x) load(x), dsel_files); % import
    
    [n_obs,n_points] = size(ds(1).dsel.y);
    avg = mean(... % mean of the observable
            reshape(... % reshape to same array
                cell2mat(... % convert to one matrix form cells
                    arrayfun(@(x) x.dsel.y, ds, 'unif', 0) ),... % get all data
            [n_obs,n_points,n_subs]),...
          3);
      
    err = max(... % max value of std
            reshape(... % reshape to same array
                cell2mat(... % convert to one matrix form cells
                    arrayfun(@(x) x.dsel.s, ds, 'unif', 0) ),... % get all data
            [n_obs,n_points,n_subs]),...
          [],3);
                     
%     for subs=1:n_subs
%         dsel_file = fullfile(dsel_path,...
%         sprintf('data_nmr_selected_%s_%s.mat', dset_snames{i}, dsel_suffixes{subs}));
        
%         d = load(dsel_file); % for data generation see create_data_str.m
%         data = d.dsel; % Get data matrix from the structure: {'PO4', 'NTP', 'RNAtotal', 'RNAfolded', 'Abortives', 'NDP', 'MgPO4', 'HSQC';}
        
%     end

    %%% Calc RNA conc from NTPs. Adjust obs names
    %==============================================
    obs_name = 'NTP'; % what to extract - need to be an argument
    obs_idx = strcmp(ds(1).dsel.names, obs_name); % get index of the target observable
%         rna_length = rna_lengths( sprintf('%sx%s', dset_snames{rna_idx}(1:2), dset_snames{idx}(4)) ); % can generate INx# by using sprintf
    rna_length = rna_lengths( dset_snames{i}(1:4) );
    rna_from_NTP = (init_NTP - avg(obs_idx,:)) ./ rna_length; % calc mM RNA
    rna_from_NTP_err = err(obs_idx,:) ./ rna_length;
    
    target_obs_idx = 7;
    avg(target_obs_idx,:) = rna_from_NTP;
    err(target_obs_idx,:) = rna_from_NTP_err;
    
    names = ds(1).dsel.names;
    names{3} = 'RNA';
    names{7} = 'RNA_from_NTP';

    %%% Plot to check results
    %==============================================
    
    if plot_results
        c = i; %(i-1)*n_subs+subs;
        rows = 1; columns = 1;
        fSize = 12;
        sbSide = 250;
        SB=0;
        sb = @(x) subplot(rows,columns,x,'FontSize',fSize);

        time = ds(1).dsel.t;

        if i==1            
            % global defined at the top of the function
            if isempty(FIG), FIG=0; end; % if this is the first place to initialize - set to zero
            FIG=FIG+1;
            figure(FIG); set(figure(FIG), 'Color', repmat(0.95,1,3), 'Position', [0 300 sbSide*columns sbSide*rows]);
            
            SB = SB+1; sb(SB);
            plot(time, rna_from_NTP,'Color',colors(c, :));
            errorbar(time, rna_from_NTP, rna_from_NTP_err, 'Color',colors(c, :));
            hold on;
        else
            plot(time, rna_from_NTP,'Color',colors(c, :));
            errorbar(time, rna_from_NTP,rna_from_NTP_err, 'Color',colors(c, :));
        end
        if i==n_sets
            axis tight;
            legend(dset_snames,'Location','best');
        end
    end
        
    
    %%% Generate IVTNMR data-structure
    %===================================
    %%% Common
    s.sname = dset_snames{i};
    s.name = dset_names{i};
    s.nmr_data_path = nmr_data_path;
    s.rna_name = rna_names{i};
    s.rna_length = rna_length;
    s.prot = prot_names{i};
    s.cotr = cotr{i};
    s.time0 = getTime0(nmr_data_path,dset_names{i});
    
    %%% 31P
    s.p31.notes = '31P [obs,std] are [MEAN,MAX] from combining 4 CDS files. RNA0 length set to 5nt!!';
    s.p31.time = ds(1).dsel.t;
    s.p31.expnos = getExpnos(nmr_data_path,dset_names{i},'P31');
%     s.p31.peak_names = [];
%     s.p31.peak_integrals = [];
%     s.p31.peak_integrals_std = [];
    
    s.p31.obs_names = names;
    s.p31.obs = avg;
    s.p31.obs_std = err;
    s.p31.obs_units = 'mM';
    
    %%% HN
    s.hn.expnos = getExpnos(nmr_data_path,dset_names{i},'HN');
    s.hn.time = getNMRTime(fullfile(nmr_data_path,dset_names{i}), s.hn.expnos, s.time0);
    
    %%% Save
    save(strcat(datasave_folder,'ivtnmr_',dset_snames{i},'.mat'),'s');
    if save_to_global_ivtnmr
        save( fullfile(nmr_data_path, dset_names{i}, 'ivtnmr.mat') ,'s');
    end
end


end