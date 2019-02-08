function res = odeM(model,tspan,p,options)
% res = odeM(model,tspan,p,[options])
%
% Solves an ODE system and the corresponding sensitivity equations. Uses
% the matlabs ode15s stiff integrator. Additionally calculates the
% observations and corresponding sensitivities.
%
% Input arguments:
%
% model = model structure as returned by model4fit
% tspan = vector with integration time points
% p = parameter vector. The parameter vector includes the initial values
%     for integrating the ODE system as well as the parameters of the ODE
%     system and the observation function first. It has the following
%     structure p = [k1 ... km x01 ... x0n]. The order of x0 and k1
%     should correspond to the order defined in the model structure
% options (optionally) = options passed to the ode15s integrator (only
%                        for integration of the ODE system)
%
% odeM returns a structure with the following fields:
%
% t = integration time points (same as tspan)
% x = state values at each time point
% y = observations at each time point
% dxdx0 = state sensitivities with respect to initial values
% dxdp = state sensitivities with respect to parameters
% dydx0 = observation sensitivities with respect to initial values
% dydp = observation sensitivities with respect to parameters
%
% Additionally the following fields are generated which can be used for
% plotting a smooth trajectory
%
% tr = integration time points
% xr = state values at each time point in tr
% yr = observations at each time point in tr

%========================================================
% TODO - flag to fit from log-fit to normal
%========================================================

nvar = model.nvar;
npar = model.npar;
nobs = model.nobs;

if isfield(model,'logparam')
    p(1:npar) = 10.^p(1:npar);
end

% initial conditions
x0 = p(model.npar+[1:model.nvar]);
S0 = zeros(nvar,npar);
J0 = diag(diag(ones(nvar)));
JS0 = [S0(:); J0(:)];

% integrate ODE
if nargin == 4
  options = odeset(options,'Jacobian',model.hjac);
else
  options = odeset('Jacobian',model.hjac);
end
solODE = ode15s(model.hrhs,tspan,x0,options,p);
res.t = tspan;
res.x = deval(solODE,tspan);
res.y = feval(model.hobs,res.t,res.x,p);

% resolved trajectory (just for plotting reasons)
res.tr = solODE.x;
res.xr = solODE.y;
res.yr = feval(model.hobs,res.tr,res.xr,p);

% integrate jacobian and sensitivities => dxdx0, dxdp
options = odeset('Jacobian',@helper2);
solS = ode15s(@helper1,tspan,JS0,options,p,model,solODE);
S = deval(solS,tspan);
res.dxdp = reshape(S(1:nvar*npar,:),nvar,npar,length(tspan));
res.dxdx0 = reshape(S(nvar*npar+1:end,:),nvar,nvar,length(tspan));

% calculate sensitivities of the observation => dydx0, dydp
obsjac = feval(model.hobsjac,res.t,res.x,p); 
obssens = feval(model.hobssens,res.t,res.x,p); 
res.dydp = zeros(nobs,npar,length(tspan));
res.dydx0 = zeros(nobs,nvar,length(tspan));
for t=1:length(tspan) 
  res.dydp(:,:,t) = obssens(:,:,t) + obsjac(:,:,t)*res.dxdp(:,:,t);
  res.dydx0(:,:,t) = obsjac(:,:,t)*res.dxdx0(:,:,t);
end  

function dSdt = helper1(t,S,p,model,sol)
% sensitivity equations
x = deval(sol,t);
dfdp = feval(model.hsens,t,x,p); 
dfdx = feval(model.hjac,t,x,p); 
S = reshape(S,model.nvar,model.nvar+model.npar);
dSdt = [dfdp zeros(model.nvar)] + dfdx*S;
dSdt = dSdt(:);

function J = helper2(t,S,p,model,sol)
% Jacobian of the sensitivity equations
x = deval(sol,t);
J = sparse(kron(diag(diag(ones(model.nvar+model.npar))),feval(model.hjac,t,x,p)));