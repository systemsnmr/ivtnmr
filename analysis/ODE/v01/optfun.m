function [resid,J] = optfun(popt,model,data,fit,options)
% [resid,J] = optfun(popt,model,data,fit,options)
%
% Helper function passed to lsqnonlin. Integrates ODE and sensitivity
% equations and calulates residuals and gradient of residuals.

resid = []; 
J = [];

p = fit.ptot;
p(fit.actp) = popt;

tdat = data.t;
ydat = data.y;
err = 1./data.s;
err(isinf(err)) = 1; % exclude zero std values
 
% integrate model  
int = odeM(model,tdat,p);

% calculate residuals
resid = (int.y-ydat).*err;
resid = resid(:);
  
% collect Jacobian of residuals
J = horzcat(int.dydp,int.dydx0);
J = J(:,fit.actp,:); % exclude inactive parameters and arrange to 
                     % match the structure of the residual vector
J = permute(J,[1 3 2]);
J = reshape(J,model.nobs*length(int.t),length(fit.actp));
J = repmat(err(:),1,length(fit.actp)).*J;

if options.output == 1
  % print sum of squared residuals only
  fprintf(1,'(res2 = %.2f) \r',resid'*resid)
elseif options.output == 2
  % print current parameter set
  fprintf(1,'%.4f\t',popt)
  fprintf(1,'(res2 = %.2f) \r',resid'*resid)
elseif options.output == 3
  % print current parameter set and plot solution
  fprintf(1,'%.4f\t',popt)
  fprintf(1,'(res2 = %.2f) \r',resid'*resid)
  plotSolution(int.tr,int.yr,tdat,ydat)
else
  % do nothing
end  

