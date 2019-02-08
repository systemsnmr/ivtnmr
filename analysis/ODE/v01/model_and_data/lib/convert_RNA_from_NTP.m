function data = convert_RNA_from_NTP(data, source_name, target_name, target_new_name, init_NTP_conc)

%%% Notes:
% - not using scaling by RNA length here!!

%% Params to run "locally"
%===========================
if nargin == 0    
    fprintf(1,'No args in the input, Running %s with test example.\n\n', mfilename);

    flag_plot_before_after = 1;

    %%% DATA
    modelSubFolder = 'v07_NTP_PO4_MgPO4_RNA';    
    dataset_names = {'avg_PNM', 'avg_pure', 'PNM', 'pure'};
    dataset = 2;

    data_files = cellfun(@(x) sprintf('data_nmr_selected_%s.mat',x), dataset_names,...
        'UniformOutput', false);
    
    dataFile = fullfile('..', modelSubFolder, data_files{dataset} );
    
    d = load(dataFile); % for data generation see create_data_str.m
    data = d.dsel; % Get data matrix from the structure: {'PO4','NTP','SMNtotal','SMNfolded','Abortives','NDP','MgPO4','HSQC';}
    
    %%% SETTINGS    
    source_name = 'NTP'; % what to extract - need to be an argument
    target_name = 'SMNtotal';
    target_new_name = 'RNA_from_NTP';
    init_NTP_conc = 20; % mM
else
    flag_plot_before_after = 0;
end

if flag_plot_before_after
    plot(data.t, data.y(3,:),'bo');
    hold on;
end

%% Actual script
%===========================

    %%% Calc RNA conc from NTPs. Adjust obs names
    %==============================================
    % Get data matrix from the structure: {'PO4','NTP','SMNtotal','SMNfolded','Abortives','NDP','MgPO4','HSQC';}    
    NTP_idx = strcmp(data.names, source_name);
    RNA_idx = strcmp(data.names, target_name);
    
    rna_from_NTP = (init_NTP_conc - data.y(NTP_idx,:)); % calc mM RNA
    rna_from_NTP_err = data.s(NTP_idx,:);
    
    data.y(RNA_idx,:) = rna_from_NTP;
    data.s(RNA_idx,:) = rna_from_NTP_err;
    
    data.names{RNA_idx} = target_new_name;

if flag_plot_before_after
    plot(data.t, data.y(3,:),'ro');
    legend({target_name, target_new_name}, 'Location', 'Best', 'Interpreter', 'None');
    axis tight; 
end

end

