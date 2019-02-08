% M-file for model SMN created by BioNetGen 2.1.7
function [t_out,obs_out,x_out]= xIVT(tend)

Nspecies=9;
Nreactions=8;
obs_out=zeros(1);
% Parameters
kcat=1;
kcatPPi=100;
kprecip=1;
kdissolve=1;
kdegNTP=0.0001;
kabortive=1;
kon=1;
koff=0.1;
NTP_0=20;
NDP_0=1;
Prot_0=0.15;
PO4=0;
NTP=0;
SMN=0;
Aborts=0;
NDP=0;
MgPO4=0;
Prot_bound_H33=0;
Prot_bound_R75=0;

% Intial concentrations
x0= [ NTP_0; 0; 0; 0; 0; NDP_0; 0; Prot_0; 0;];

% Reaction flux function
function f= flux(t,x)
    f(1)=kcat*x(1);
    f(2)=kabortive*x(1);
    f(3)=kcatPPi*x(4);
    f(4)=kprecip*x(3);
    f(5)=kdissolve*x(5);
    f(6)=kdegNTP*x(1);
    f(7)=kon*x(2)*x(8);
    f(8)=koff*x(9);
end

% Derivative function
function d= xdot(t,x)
    f=flux(t,x);
    d(1,1)=  -f(1) -f(2) -f(6);
    d(2,1)=  +f(1) -f(7) +f(8);
    d(3,1)=  +2*f(3) -f(4) +f(5) +f(6);
    d(4,1)=  +f(1) +f(2) -f(3);
    d(5,1)=  +f(4) -f(5);
    d(6,1)=  +f(6);
    d(7,1)=  +f(2);
    d(8,1)=  -f(7) +f(8);
    d(9,1)=  +f(7) -f(8);
end

% Observables
function o= obs(x)
  w= [
      0 0 1 0 0 0 0 0 0;
      1 0 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0 1;
      0 0 0 0 0 0 1 0 0;
      0 0 0 0 0 1 0 0 0;
      0 0 0 0 1 0 0 0 0;
      0 0 0 0 0 0 0 0 1;
      0 0 0 0 0 0 0 0 1;
     ];
  o= x*w';
end

onames={'PO4','NTP','SMN','Aborts','NDP','MgPO4','Prot\_bound\_H33','Prot\_bound\_R75'};

% Integrate ODEs
[t_out,x_out]= ode15s(@xdot, [0 tend], x0);
% compute observables
obs_out= obs(x_out);
% plot observables
plot(t_out,obs_out);
legend(onames);
end
