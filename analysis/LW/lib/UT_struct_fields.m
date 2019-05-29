function o = UT_struct_fields(s1,s2)

    if nargin==0
        analysis_root = fullfile(fileparts(which(mfilename)),'..');
        scratchdir = 'datasave';
        s1 = load(fullfile(analysis_root,scratchdir,'a_fit_LW_UT1.mat'));
        s1 = s1.d;
        s2 = s1;
%         class(s1)
%         isstruct(s1)
%         iscell(s1)
    end
    
    if isstruct(s1)
        s1_cells = {s1}; % encapsulate into a cell - to generalize the below script.
        s2_cells = {s2};
    elseif iscell(s1) && numel(s1)>1 % flatten cell array
        % TODO: while condition - to continue flattening until reaching
        % isstruct(s1_cells{1})
        s1_cells = vertcat(s1{:});
        s2_cells = vertcat(s2{:});
    end
    
    n_cells = numel(s1_cells);
    
    for iCell=1:n_cells 
        s1 = s1_cells{iCell};
        s2 = s2_cells{iCell};                        

        d_field_names = fieldnames(s1);
        n_fields = length(d_field_names);

%         s1.time{1} = 735.9429;
        comp = @(s1,s2,field) isequal(s1.(field), s2.(field)); % comparison test

        if ~isequal(s1,s2)
            warning('=== UnitTest: structures (in cell %i) dont match', iCell); % overall check
        end;
        for iField = 1:n_fields
            iField_name = d_field_names{iField};
            if ~comp(s1, s2, iField_name)
                if strcmp(iField_name, 'Info')
                    disp('Only Info field mismatch.');
                else
                    warning('=== UnitTest: %s field unequal.\n', iField_name);
                end;
            end;
        end;

    end;
    
end