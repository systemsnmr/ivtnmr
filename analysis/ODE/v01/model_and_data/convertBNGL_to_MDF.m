function convertBNGL_to_MDF(fname)

[p,name,e] = fileparts(fname); 

%%%%%%%%%%%%%%
% read net file
%%%%%%%%%%%%%%

fid = fopen(sprintf('%s.net',name),'r');
if fid == -1
  error('Can''t read net-file \n')
end
body = fscanf(fid,'%c');
fclose(fid);

% get species, indices and compartments 
var = regexp(body,'begin species\n(.*)\nend species','tokens');
if isempty(var)
  error('Break !')
end
var = regexp(var{1},'\n','split');
model.nvar = length(var{1});
for v=1:model.nvar
  varf = regexp(var{1}{v},'\s+','split');
  model.all_var_names{v} = varf{3};
  m = regexp(varf{3},'Loc~([\w]+)','tokens');
  if ~isempty(m)
    model.var_comp{v} = m{1}{1};
  end
end

% IF THERE IS PRODUCTION (Prd) OR DEGRADATION (Deg) , REMOVE THESE STATES
% FROM THE VARIABLES, AND MARK THEIR ROWS:

prd_idx = find(strcmp(model.all_var_names,'Prd()'));
deg_idx = find(strcmp(model.all_var_names,'Deg()'));
prd_deg_idx = [0 0];

% IF BOTH ARE EMPTY, JUST DO NOTHING
if isempty(prd_idx)&&isempty(deg_idx)
% ONE IS EMPTY, REMOVE THE OTHER ONE
elseif isempty(prd_idx) || isempty(deg_idx)
    if ~isempty(prd_idx)
       model.all_var_names(prd_idx) = [];
       prd_deg_idx(1) = 1;
       model.nvar = model.nvar-1;
    elseif ~isempty(deg_idx)
       model.all_var_names(deg_idx) = [];
       prd_deg_idx(2) = 1;
       model.nvar = model.nvar-1;
    end
% IF BOTH ARE NOT EMPTY, ORDER THEM AND REMOVE THEM
elseif (~isempty(prd_idx)) && (~isempty(deg_idx))
    if prd_idx > deg_idx
        model.all_var_names(prd_idx) = [];
        model.all_var_names(deg_idx) = [];
        prd_deg_idx(1:2) = 1;
        model.nvar = model.nvar-2;
    elseif deg_idx > prd_idx
        model.all_var_names(deg_idx) = [];
        model.all_var_names(prd_idx) = [];
        prd_deg_idx(1:2) = 1;
        model.nvar = model.nvar-2;
    end
end


% BUT!! WHAT IF there are no compartments??
% then simply : size(fieldnames(model),1) = 2
% or 
% find(strcmp(fieldnames(model),'var_comp'))
% ans =
%    Empty matrix: 0-by-1

% get molecule types
mols = regexp(body,'begin molecule types\n(.*)\nend molecule types','tokens');
if isempty(mols)
  error('Break !')
end
mols = regexp(mols{1},'\n','split');
model.nmol = length(mols{1});
for v=1:model.nmol
  molsf = regexp(strtrim(mols{1}{v}),'\s+','split');
  model.molecule_types{v} = molsf{2};
end

% Remove the Prd() Deg() from the molecule types:
if prd_deg_idx(1)==1
    model.molecule_types = setdiff(model.molecule_types,'Prd()');
    model.nmol = model.nmol-1;
end

if prd_deg_idx(2)==1
    model.molecule_types = setdiff(model.molecule_types,'Deg()');
    model.nmol = model.nmol-1;
end


% get observations
obs = regexp(body,'begin observables\n(.*)\nend observables','tokens');
if isempty(obs)
  error('Break !')
end
obs = regexp(obs{1},'\n','split');
model.nobs = length(obs{1});
for v=1:model.nobs
  obsf = regexp(strtrim(obs{1}{v}),'\s+','split');
  model.obs_types{v} = obsf{2};
  model.obs_names{v} = obsf{3};
  if length(obsf)>3
    model.obs_pattern{v} = obsf{4};
  end
end
% get parameters
pars = regexp(body,'begin parameters\n(.*)\nend parameters','tokens');
if isempty(pars)
  error('Break !')
end
pars = regexp(pars{1},'\n','split');
model.nallpar = length(pars{1});




for v=1:model.nallpar
  parsf = regexp(strtrim(pars{1}{v}),'\s+','split');
  % REMOVE THE PARAMETERS THAT ARE INITIAL CONDITIONS FOR PRODUCTION AND
  % DEGRADATION
  if (~strcmp(parsf{2},'Prd0')) && (~strcmp(parsf{2},'Deg0'))
      
  model.all_par_names{v} = parsf{2};
  model.all_par_values(v) = eval(parsf{3});
  if strcmp(model.all_par_names{v},'L0')
    model.L0idx_par = v;
  end
  
  else
      if strcmp(parsf{2},'Prd0') == 1
          Prd0_idx = v;
          model.nallpar = model.nallpar-1;
      end
      if strcmp(parsf{2},'Deg0') == 1
          Deg0_idx = v;
          model.nallpar = model.nallpar-1;
      end
  end    
end
  
% IF there are compartments...
if ~isempty(find(strcmp(fieldnames(model),'var_comp'), 1))
% match volume parameters to compartments
comp = unique(model.var_comp(~cellfun(@isempty,model.var_comp)));
model.var_comp_par = zeros(model.nvar,1);
    for c=1:length(comp)
      idx = find(strcmp(sprintf('vol_%s',comp{c}),model.all_par_names));
      model.var_comp_par(strcmp(model.var_comp,comp{c})) = idx;
    end

end

% get fluxes: which flux is generated by which rule?
flux = regexp(body,'begin reactions\n(.*)\nend reactions','tokens');
if isempty(flux)
  error('Break !')
end
flux = regexp(flux{1},'\n','split');
model.nflux = length(flux{1});
for v=1:model.nflux
  fluxf = regexp(strtrim(flux{1}{v}),'\s+','split');
  rule = regexp(fluxf{5},'#([\w]+)','tokens');
  model.flux_rule{v} = rule{1}{1};
end


%%%%%%%%%%%%%%
% read m-file
%%%%%%%%%%%%%%

fid = fopen(sprintf('%s.m',name),'r');
if fid == -1
  error('Can''t read m-file \n')
end
body = fscanf(fid,'%c');
fclose(fid);
% get rateLaw expressions
% IF any!!
rl = regexp(body,'rateLaw[\d]+=[\w\d\/\*]+;\n','match');
if ~isempty(rl)
model.rateLaws = []; 
    for r=1:length(rl)
      model.rateLaws{r} = strtrim(rl{r});
    end

end

% get initial conditions
init = regexp(body,'x0= \[(.*)\];\n\n% Reaction flux function','tokens');
if isempty(init)
  error('Break !')
end
init = regexp(init{1},';','split');
model.idx_x0p = []; 
for v=1:length(init{1})
  s = strtrim(init{1}{v});
  if ~isempty(s)
    if (~strcmp(s,'Prd0')) && (~strcmp(s,'Deg0'))
        match = strmatch(s,model.all_par_names);
        if ~isempty(match)
          model.idx_x0p = [model.idx_x0p;v match];
          if strcmp(s,'L0')
            model.L0idx_var = v;
          end
        end
    end
  end
end

% get fluxes 
flux = regexp(body,'function f= flux\(t,x\)\n(.*)\nend\n\n% Derivative function','tokens');
if isempty(flux)
  error('Break !')
end
flux = regexp(flux{1},'\n','split');
trans_flux = [];
for v=1:length(flux{1})
  model.fluxes{v} = strtrim(flux{1}{v});
  ok = 0;
  
%   if ~isempty(trans_par)
%   
%       for p=1:length(trans_par)
%         fluxf = regexp(model.fluxes{v},trans_par{p},'match');
%         if ~isempty(fluxf)
%           ok = 1;
%         end
%       end
%       
%   
%   
%       if ok
%         trans_flux = [trans_flux v];
%       end
% 
%   end
end
model.trans_flux = trans_flux;

% FIND THE FLUX/FLUXES THAT REFER TO THE PRODUCTION VARIABLE AND REMOVE IT
% BY DELETING * x(?) JUST LETTING EFFECTIVELY "1"
model.fluxes = strrep(model.fluxes,['*x(',num2str(prd_idx),')'],'');

% Sometimes the added species for production/degradation are somewhere
% in the middle of the x's. Correct the fluxes by substituting all the x's
% that are larger than deg_idx/prd_idx with x-1 sequentially

% IF BOTH ARE EMPTY, JUST DO NOTHING
if isempty(prd_idx)&&isempty(deg_idx)
% ONE IS EMPTY, MAKE ALL THE X'S FROM THAT ON TO x-1 where nexessary
elseif isempty(prd_idx) || isempty(deg_idx)
    if ~isempty(prd_idx)
       for i=1:model.nflux
           for j=prd_idx+1:model.nvar+1
           model.fluxes{i} = strrep(model.fluxes{i},['*x(',num2str(j),')'],['*x(',num2str(j-1),')']);
           end
       end
    elseif ~isempty(deg_idx)
       for i=1:model.nflux
           for j=deg_idx+1:model.nvar+1
           model.fluxes{i} = strrep(model.fluxes{i},['*x(',num2str(j),')'],['*x(',num2str(j-1),')']);
           end
       end
    end
% IF BOTH ARE NOT EMPTY, ORDER THEM AND REMOVE THEM
elseif (~isempty(prd_idx)) && (~isempty(deg_idx))
    if prd_idx > deg_idx
       for i=1:model.nflux
           for j=prd_idx+1:model.nvar+1
           model.fluxes{i} = strrep(model.fluxes{i},['*x(',num2str(j),')'],['*x(',num2str(j-1),')']);
           end
       end
       for i=1:model.nflux
           for j=deg_idx+1:model.nvar+1
           model.fluxes{i} = strrep(model.fluxes{i},['*x(',num2str(j),')'],['*x(',num2str(j-1),')']);
           end
       end
    elseif deg_idx > prd_idx
       for i=1:model.nflux
           for j=deg_idx+1:model.nvar+1
           model.fluxes{i} = strrep(model.fluxes{i},['*x(',num2str(j),')'],['*x(',num2str(j-1),')']);
           end
       end
       for i=1:model.nflux
           for j=prd_idx+1:model.nvar+1
           model.fluxes{i} = strrep(model.fluxes{i},['*x(',num2str(j),')'],['*x(',num2str(j-1),')']);
           end
       end
    end
end

% get RHS and stoichiometric matrix
rhs = regexp(body,'f=flux\(t,x\);\n(.*)\nend\n\n% Observables','tokens');
if isempty(rhs)
  error('Break !')
end
rhs = regexp(rhs{1},'\n','split');
model.stoich_mat = sparse(model.nvar,model.nflux);

% for v=1:length(rhs{1})
% WE HAVE TO EXCLUDE THE PRODUCTION DEGRADATION RHSides
for v=setdiff([1:length(rhs{1})], [prd_idx deg_idx])
  rhsl = regexp(strtrim(rhs{1}{v}),'=','split');
  rhsl = strtrim(rhsl{2});
  % modify compartment exchange fluxes
  for f=1:length(trans_flux)
    e = regexp(rhsl,sprintf('f\\(%d\\)',trans_flux(f)),'end');
    if ~isempty(e)
      rhsl = strcat(rhsl(1:e),sprintf('/vol_%s',model.var_comp{v}),rhsl(e+1:end));
    end
  end
  model.rhs{v} = rhsl;
  % build stoichiometric matrix + volume matrix (i.e. stoichiometric
  % matrix where the transport fluxes wille be scaled by the
  % corresponding volume)
  %f = regexp(model.rhs{v},'([\+-])f\(([\d]+)\)','tokens');
  %for n=1:length(f)
  %  fidx = eval(f{n}{2});
  %  model.stoich_mat(v,fidx) = eval(sprintf('%s1',f{n}{1}));
  %end
  %% 
  f = regexp(model.rhs{v},'\s','split');
  
  % in the case where this is just ";", the elements on the stoichiometric
  % matrix are zero
  
  if ~strcmp(f(1),';')==1
      for n=1:length(f)
         r = regexp(f{n},'f\(([\d]+)\)','tokens');
         fidx = eval(r{1}{1});
         r = regexp(f{n},'f','split');
         r = strtok(r{1},'*');
         if strcmp(r,'+')
            model.stoich_mat(v,fidx) = 1;
         elseif strcmp(r,'-')
            model.stoich_mat(v,fidx) = -1;
         else
            model.stoich_mat(v,fidx) = eval(r);
         end
      end
  else
      
  end
end

% remove the empty lines if any (because of production/dengradation)

model.stoich_mat(find(cellfun(@isempty,model.rhs)),:) = [];
model.rhs(find(cellfun(@isempty,model.rhs))) = [];


% get Obs
obs = regexp(body,'w= \[\n(.*)\n\s+\];','tokens');
if isempty(obs)
  error('Break !')
end
model.obs = eval(sprintf('[%s]',obs{1}{1}));
% EXCLUDE COLUMNS OF SPECIES PRODUCTION/DEGRADATION
model.obs(:,[prd_idx deg_idx]) = [];



% *********************************************************


% EXTRACT ALL THE Right Hand SIDES INTERMS OF INDIVIDUAL FLUXES
RHS_ind_fluxes = cell(length(model.rhs),1);
for i=1:length(model.rhs)
regexp(model.rhs{i},' ','split');
RHS_ind_fluxes{i} = regexp(model.rhs{i},' ','split');
% remove the ";" from the last flux:
RHS_ind_fluxes{i}(end) = regexprep(RHS_ind_fluxes{i}(end),';','');
end


% EXTRACT THE EXPRESSION OF EACH FLUX
flux_expression = cell(model.nflux,1);
for i=1:model.nflux
[flux{i},flux_expression{i}] = strtok(model.fluxes{i},'=');
flux_expression{i} = regexprep(flux_expression{i},'[=;]','');
end


% IN RHS, SUBSTITUTE EVERY FLUX COMPONENT BY THE CORRESPONDING EXPRESSION
for i=1:length(model.rhs)
   for j=1:length(RHS_ind_fluxes{i})
   % we want for every flux component to be substituted by the expression 
   % So, e.g.in RHS_ind_fluxes{6}{1} which is "+(-2)*f(1)", we want
   % that the f(1) is detected, and substituted with:
   % flux_expression(find(strcmp(flux,'f(1)')))
   f_expr = RHS_ind_fluxes{i}{j}(cell2mat(regexp(RHS_ind_fluxes{i}(j),'f')):end);
   f_expr_mod = regexprep(f_expr,'\(','\\(');
   f_expr_mod = regexprep(f_expr,'\)','\\)');
   RHS_ind_fluxes{i}{j} = regexprep(RHS_ind_fluxes{i}{j},f_expr_mod,flux_expression(find(strcmp(flux,f_expr)))); 
   end   
end

% MERGE RHSides in one cell for each variable:
RHS_expressions = cell(length(model.rhs),1);
for i=1:length(model.rhs)
    RHS_expr = [];
   for j=1:length(RHS_ind_fluxes{i})
      RHS_expr = cat(2,RHS_expr,RHS_ind_fluxes{i}{j});
   end
   RHS_expressions{i} = RHS_expr;
end


% REMOVE THE RIGHT AND LEFT PARENTHESIS IN THESE EXPRESSIONS
% SO THAT WE ONLY HAVE THE x1... xn as variables
for i=1:length(model.rhs)
  RHS_expressions{i} = regexprep(RHS_expressions{i},'\(','');
  RHS_expressions{i} = regexprep(RHS_expressions{i},'\)','');
end


% Save the expressions in the structure
model.rhs_exp = RHS_expressions;

%% ADD SOME LIMITS FOR MIN AND MAX VALUES:
% Define npar!!!
model.npar = size(model.all_par_names,2)-size(model.idx_x0p,1);

idx_kinpar = setdiff(1:model.npar,model.idx_x0p(:,2));
P = model.all_par_values(idx_kinpar);
model.par_names = model.all_par_names(idx_kinpar);
model.par_init  = P;
model.par_min   = zeros(size(P));
model.par_max   = 100*P;

% Get the initial conditions
x0 = zeros(size(model.all_var_names,2),1)';
x0(model.idx_x0p(:,1)) = model.all_par_values(model.idx_x0p(:,2));

model.var_names = cellfun(@(x) ['x',num2str(x)],num2cell(1:model.nvar),'UniformOutput',false);
model.var_init  = x0;
model.var_min   = zeros(size(x0));
model.var_max   = 100*x0;

%% WRITE MDF FILE

% VARIABLES
txt_var = [];
for n=1:model.nvar
  txt_var = sprintf('%s%s  %s  %s  %s    # %s\n',txt_var,['x',num2str(n)],...
                                            num2str(model.var_init(n)),...
                                            num2str(model.var_min(n)),...
                                            num2str(model.var_max(n)),...
                                            model.all_var_names{n});
end

% PARAMETER NAMES
txt_par = [];
for n=1:length(model.par_names)
  txt_par = sprintf('%s%s  %s  %s  %s    \n',txt_par,model.all_par_names{n},...
                                        num2str(model.par_init(n)),...
                                        num2str(model.par_min(n)),...
                                        num2str(model.par_max(n)));
end

if ~isempty(find(strcmp(fieldnames(model),'rateLaws'), 1))
    for r=1:length(model.rateLaws)
      txt_par = sprintf('%s%s\n',txt_par,model.rateLaws{r});
    end
end


% FLUXES
txt_flux = [];
for n=1:size(model.fluxes,2)
  txt_flux = sprintf('%s%s\n',txt_flux,...
      regexprep(regexprep(model.fluxes{n},'[\(\);]',''),'=',' = '));
end

outfile = sprintf('%s.mdf',name);

try
  fid = fopen(outfile,'w');
catch
  error('Can''t open output file!')
end

% % RHS (in terms of fluxes)
% txt_rhs = [];
% for n=1:size(RHS_ind_fluxes,1)
%     nth_flux = RHS_ind_fluxes{n};
%   txt_rhs = sprintf('%s%s\n',txt_rhs,regexprep(strcat(nth_flux{:}),'[\(\)]',''));
% end

%  RHS (in terms of fluxes)
txt_rhs = [];
for n=1:size(model.rhs,2)
%     nth_flux = model.rhs{n};
  txt_rhs = sprintf('%s%s\n',txt_rhs,regexprep(model.rhs{n},'[\(\);]',''));
end

% OBSERVABLES
txt_obs = [];
for n=1:size(model.obs,1)
    idx_obs = find(model.obs(n,:));
    txt_ind_obs = [];
    for j=1:size(idx_obs,2)
    txt_ind_obs = sprintf('%s+%s',txt_ind_obs,['x',num2str(idx_obs(j))]);
    end
  txt_obs = sprintf('%s%s    # %s\n',txt_obs,txt_ind_obs,model.obs_pattern{n});
end



fwrite(fid,sprintf('<DOC>\n'));
fwrite(fid,sprintf('...add text...\n'));
fwrite(fid,sprintf('</DOC>\n\n'));

fwrite(fid,sprintf('<VAR>\n'));
fwrite(fid,txt_var);
fwrite(fid,sprintf('</VAR>\n\n'));

fwrite(fid,sprintf('<PAR>\n'));
fwrite(fid,txt_par);
fwrite(fid,sprintf('</PAR>\n\n'));

fwrite(fid,sprintf('<DEF>\n'));
fwrite(fid,txt_flux);
fwrite(fid,sprintf('</DEF>\n\n'));

fwrite(fid,sprintf('<RHS>\n'));
fwrite(fid,txt_rhs);
fwrite(fid,sprintf('</RHS>\n\n'));

fwrite(fid,sprintf('<OBS>\n'));
fwrite(fid,txt_obs);
fwrite(fid,sprintf('</OBS>\n'));

fclose(fid);



end