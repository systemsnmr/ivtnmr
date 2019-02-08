function model = model4fit(mdfile)
% model = model4fit(mdfile)
%
% @YarN 2016-10-26 - modified to also save BNGL-based proper variable names
%
% This function reads a model definition file and calculates the derivatives
% of the right hand side (rhs) of the ODE and the observation function with
% respect to the state variables and the parameters. This step requires
% access to the symbolic toolbox of MATLAB. It generates several output
% files for the ODE rhs, the observation function and the respective
% Jacobian and sensitivities. It returns a model structure which can be
% passed to the other fitting routines (fitModel,odeM).
%
% Note that you need write access to the current working directory!
%
%
% The general structure of a model definition file follows below:
%
% Variables are declared as
%
% <VAR>
% var_name1  initial_value  lower_bound  upper_bound
% var_name2  initial_value  lower_bound  upper_bound
% ...
% </VAR>
% 
% Parameters are declared as
%
% <PAR>
% par_name1  initial_value  lower_bound  upper_bound
% par_name2  initial_value  lower_bound  upper_bound
% ...
% </PAR>
% 
% Variable and parameter names should conform to MATLABs naming convention
% concerning symbols. Floating point numbers for the initial value, lower
% and upper bounds can be written such as 0.01, 10^-2, 1e-2 or 1E-2.
%
% The rhs of the ODE is defined as 
%
% <RHS>
% rhs_eqn1
% rhs_eqn2
% ...
% </RHS>
%
% Note that the order of equations in the rhs should correspond to the order
% of variable declations in order to match each variable with its
% corresponding equation. An rhs equation may contain any state variable and
% any parameter which was defined. Additionally, it can contain any
% algebraic expression the symbolic toolbox understands. Use 't' in order to
% indicate an explicite dependance on time, e.g., exp(-p1*t). The same
% freedom applies to the observation function which is defined as
%
% <OBS>
% obs_eqn1
% obs_eqn2
% ...
% </OBS>
%
% A model definition file may contain a definition. Any occurance of the
% left hand side (lhs) in either the ODE rhs or the observation functions is
% replaced by the rhs of the definition. Definitions can make use of
% preceding definitions.
%
% <DEF>
% lhs1 = rhs1
% lhs2 = rhs2
% lhs3 = lhs1 + lhs2
% ...
% </DEF>
%
% Documentation can be included as
% 
% <DOC>
% Documentation here ...
% </DOC>
%
% To comment a line (or the remainder of a line) use either matlab's
% comment symbol '%' or the hash symbol '#' 

doc_start = '<DOC>';
doc_end = '</DOC>';

par_start = '<PAR>';
par_end = '</PAR>';

var_start = '<VAR>';
var_end = '</VAR>';

def_start = '<DEF>';
def_end = '</DEF>';

rhs_start = '<RHS>';
rhs_end = '</RHS>';

obs_start = '<OBS>';
obs_end = '</OBS>';

if ~license('test','symbolic_toolbox')
  error('ERROR: this code depends on the symbolic toolbox!')
  return
end

fid = fopen(mdfile,'r');
if fid == -1
  error(sprintf('Can''t read file: %s \n',mdfile))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read file contents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'Reading input file ...\n')
body = fscanf(fid,'%c');
fclose(fid);
model = struct;

%% documentation
doc = regexp(body,sprintf('%s(.*)%s',doc_start,doc_end),'tokens');
model.doc = strtrim(doc{1}{1});

%% model parameter
par = regexp(body,sprintf('%s\n(.*)\n%s',par_start,par_end),'tokens');
if isempty(par)
  error('No parameters given!')
end
par = regexp(par{1},'\n','split');
NPar = length(par{1});
model.npar = NPar;
for p=1:NPar
  parf = regexp(par{1}{p},'\s+','split');
  model.par_names{p} = parf{1};
  model.par_init(p) = estr2num(parf{2});
  model.par_min(p) = estr2num(parf{3});
  model.par_max(p) = estr2num(parf{4});
end

%% variable names
var = regexp(body,sprintf('%s\n(.*)\n%s',var_start,var_end),'tokens');
if isempty(var)
  error('No variables given!')
end
var = regexp(var{1},'\n','split');
NVar = length(var{1});
model.nvar = NVar;
for v=1:NVar
  varf = regexp(var{1}{v},'\s+','split');
  model.var_names{v} = varf{1};
  model.var_init(v) = estr2num(varf{2});
  model.var_min(v) = estr2num(varf{3});
  model.var_max(v) = estr2num(varf{4});
  model.var_names_full{v} = varf{6}; % @YarN 2016-10-26.
end

%% definitions
NDef = 0;
def = regexp(body,sprintf('%s\n(.*)\n%s',def_start,def_end),'tokens');
if ~isempty(def)
  def = regexp(def{1},'\n','split');
  NDef = length(def{1});
  model.def_lhs = {[]};%cell(NDef,1);
  model.def_rhs = {[]};%cell(NDef,1);
  for a=1:NDef
    deff = regexp(def{1}{a},'=','split');
    lhs = strtrim(deff{1});
    deff = regexp(deff{2},'[%#]','split');
    rhs = strtrim(deff{1});
    rhs = subs(rhs,model.def_lhs,model.def_rhs,0);
    model.def_lhs{a} = lhs;
    model.def_rhs{a} = rhs;
  end
end

%% Right hand side of ode
rhs = regexp(body,sprintf('%s\n(.*)\n%s',rhs_start,rhs_end),'tokens');
rhs = regexp(rhs{1},'\n','split');
NRHS = length(rhs{1});
model.rhs = cell(NRHS,1);
for r=1:NRHS
  s = regexp(rhs{1}{r}, '[%#]', 'split');
  if ~isempty(s{1})
    if NDef>0
      res = subs(s{1},model.def_lhs,model.def_rhs,0);
    else
      res = s{1};
    end
    model.rhs{r} = char(res);
  end
end

%% observations
NObs = 0;
obs = regexp(body,sprintf('%s\n(.*)\n%s',obs_start,obs_end),'tokens');
if ~isempty(obs)
  obs = regexp(obs{1},'\n','split');
  NObs = length(obs{1});
  for o=1:NObs
    s = regexp(obs{1}{o}, '[%#]', 'split');
    if ~isempty(s{1})
      if NDef>0
        res = subs(s{1},model.def_lhs,model.def_rhs,0);
      else
        res = s{1};
      end
      model.obs{o} = char(res);
    end
  end
else
  % in case no observations are given, assume that all variables are observed
  NObs = NVar;
  for v=1:NObs
    model.obs{v} = model.var_names{v};
  end
end
model.nobs = NObs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% error check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(NRHS==NVar)
  error('Number of equations and number of variables does not match!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% symbolic calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'Calculating Jacobian and sensitivities ...\n')

if ~exist('t','var')
  syms t % always add time variable 't' to symbols
end

%% ODE Jacobian
model.dfdx = cell(NRHS,NVar);
for r=1:NRHS
  for v=1:NVar
%     model.dfdx{r,v} = char(diff(model.rhs{r},model.var_names{v},1));
      % @YaroslavN. 2014-09-12. Modified this after the warnings:
      % Warning: The method char/diff will be removed in a future release. Use
      % sym/diff instead. For example diff(sym('x^2')). After removal diff('x^2')
      % will return diff(double('x^2')).
    model.dfdx{r,v} = diff(sym(model.rhs{r}),sym(model.var_names{v}),1);
  end
end

%% ODE sensitivities
model.dfdp = cell(NRHS,NPar);
for r=1:NRHS
  for p=1:NPar
      % @YaroslavN. 2014-09-12. See note above.
%     model.dfdp{r,p} = char(diff(model.rhs{r},model.par_names{p},1));
    model.dfdp{r,p} = diff(sym(model.rhs{r}),sym(model.par_names{p}),1);
  end
end

%% Obs Jacobian
model.dgdx = cell(NObs,NVar);
for r=1:NObs
  for v=1:NVar
      % @YaroslavN. 2014-09-12. See note above.
%     model.dgdx{r,v} = char(diff(model.obs{r},model.var_names{v},1));
    model.dgdx{r,v} = diff(sym(model.obs{r}),sym(model.var_names{v}),1);
  end
end

%% Obs sensitivities
model.dgdp = cell(NObs,NPar);
for r=1:NObs
  for p=1:NPar
      % @YaroslavN. 2014-09-12. See note above.
%     model.dgdp{r,p} = char(diff(model.obs{r},model.par_names{p},1));
    model.dgdp{r,p} = diff(sym(model.obs{r}),sym(model.par_names{p}),1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'Creating output files ...\n')
%% chop extension
[p,fname] = fileparts(mdfile);
model.name = fname;

%% create preamble 
pre = sprintf('\n%% Variables \n');
for n=1:NVar
  pre = sprintf('%s%s = x(%d,:);\n',pre,model.var_names{n},n);
end
pre = sprintf('%s\n%% Parameters\n',pre);
for n=1:NPar
  pre = sprintf('%s%s = p(%d);\n',pre,model.par_names{n},n);
end
pre = sprintf('%s\n',pre);

%% RHS file
model.hrhs = eval(sprintf('@RHS%s',fname));
rhsfile = sprintf('RHS%s.m',fname);
fid = fopen(rhsfile,'w');
if fid == -1
  error(sprintf('Can''t write file: %s \n',rhsfile))
end
fwrite(fid,sprintf('function dxdt = RHS%s(t,x,p)\n',fname));
fwrite(fid,pre);
fwrite(fid,sprintf('dxdt = zeros(size(x));\n'));
for l=1:NRHS
  fwrite(fid,sprintf('dxdt(%d,:) = %s; %% RHS %s\n',l,vectorize(model.rhs{l}),model.var_names{l}));
end
fclose(fid);

%% RHS jacobian file
model.hjac = eval(sprintf('@JAC%s',fname));
jacfile = sprintf('JAC%s.m',fname);
fid = fopen(jacfile,'w');
if fid == -1
  error(sprintf('Can''t write file: %s \n',jacfile))
end
fwrite(fid,sprintf('function dfdx = JAC%s(t,x,p)\n',fname));
fwrite(fid,pre);
fwrite(fid,sprintf('dfdx = zeros(%d,size(x,2));\n',NVar*NVar));
foo = sparse(NRHS,NVar);
for n=1:NRHS*NVar
  if ~strcmp(model.dfdx{n},'0')
    foo(n) = 1;
    fwrite(fid,sprintf('dfdx(%d,:) = %s;\n',n,vectorize(model.dfdx{n})));
  end
end
fwrite(fid,sprintf('dfdx = reshape(dfdx,%d,%d,size(x,2));',NVar,NVar));
fclose(fid);
model.adjMatrix = sparse(foo'); % transpose of Jacobian !

%% RHS sensitivity file
model.hsens = eval(sprintf('@SENS%s',fname));
sensifile = sprintf('SENS%s.m',fname);
fid = fopen(sensifile,'w');
if fid == -1
  error(sprintf('Can''t write file: %s \n',sensifile))
end
fwrite(fid,sprintf('function dfdp = SENS%s(t,x,p)\n',fname));
fwrite(fid,pre);
fwrite(fid,sprintf('dfdp = zeros(%d,size(x,2));\n',NVar*NPar));
for n=1:NRHS*NPar
  if ~strcmp(model.dfdp{n},'0')
    fwrite(fid,sprintf('dfdp(%d,:) = %s;\n',n,vectorize(model.dfdp{n})));
  end
end
fwrite(fid,sprintf('dfdp = reshape(dfdp,%d,%d,size(x,2));',NVar,NPar));
fclose(fid);


%% observation file
model.hobs = eval(sprintf('@OBS%s',fname));
obsfile = sprintf('OBS%s.m',fname);
fid = fopen(obsfile,'w');
if fid == -1
  error(sprintf('Can''t write file: %s \n',obsfile))
end
fwrite(fid,sprintf('function y = OBS%s(t,x,p)\n',fname));
fwrite(fid,pre);
fwrite(fid,sprintf('y = [\n'));
for l=1:NObs
  fwrite(fid,sprintf('%s\n',vectorize(model.obs{l})));
end
fwrite(fid,sprintf('];\n'));
fclose(fid);

%% observation jacobian file
model.hobsjac = eval(sprintf('@OBSJAC%s',fname));
obsjacfile = sprintf('OBSJAC%s.m',fname);
fid = fopen(obsjacfile,'w');
if fid == -1
  error(sprintf('Can''t write file: %s \n',obsjacfile))
end
fwrite(fid,sprintf('function dgdx = OBSJAC%s(t,x,p)\n',fname));
fwrite(fid,pre);
fwrite(fid,sprintf('dgdx = zeros(%d,size(x,2));\n',NObs*NVar));
for n=1:NObs*NVar
  if ~strcmp(model.dgdx{n},'0')
    fwrite(fid,sprintf('dgdx(%d,:) = %s;\n',n,vectorize(model.dgdx{n})));
  end
end
fwrite(fid,sprintf('dgdx = reshape(dgdx,%d,%d,size(x,2));',NObs,NVar));
fclose(fid);

%% observation sensitivity file
model.hobssens = eval(sprintf('@OBSSENS%s',fname));
obssensfile = sprintf('OBSSENS%s.m',fname);
fid = fopen(obssensfile,'w');
if fid == -1
  error(sprintf('Can''t write file: %s \n',obssensfile))
end
fwrite(fid,sprintf('function dgdp = OBSSENS%s(t,x,p)\n',fname));
fwrite(fid,pre);
fwrite(fid,sprintf('dgdp = zeros(%d,size(x,2));\n',NObs*NPar));
for n=1:NObs*NPar
  if ~strcmp(model.dgdp{n},'0')
    fwrite(fid,sprintf('dgdp(%d,:) = %s;\n',n,vectorize(model.dgdp{n})));
  end
end
fwrite(fid,sprintf('dgdp = reshape(dgdp,%d,%d,size(x,2));',NObs,NPar));
fclose(fid);

foo = ls; % force silent re-reading the contents of the current directory


function n = estr2num(str)
[n,ok] = str2num(str);
if ~ok
  error(sprintf('Cannot read number: %s',str))
end
