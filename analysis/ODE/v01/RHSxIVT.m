function dxdt = RHSxIVT(t,x,p)

% Variables 
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);
x5 = x(5,:);
x6 = x(6,:);
x7 = x(7,:);
x8 = x(8,:);
x9 = x(9,:);

% Parameters
kcat = p(1);
kcatPPi = p(2);
kprecip = p(3);
kdissolve = p(4);
kdegNTP = p(5);
kabortive = p(6);
kon = p(7);
koff = p(8);

dxdt = zeros(size(x));
dxdt(1,:) = - kabortive.*x1 - kcat.*x1 - kdegNTP.*x1; % RHS x1
dxdt(2,:) = 0.17241379310344827586206896551724.*kcat.*x1 + koff.*x9 - kon.*x2.*x8; % RHS x2
dxdt(3,:) = kdegNTP.*x1 + 2.*kcatPPi.*x4 + kdissolve.*x5 - kprecip.*x3; % RHS x3
dxdt(4,:) = kabortive.*x1 + 0.82758620689655172413793103448276.*kcat.*x1 - kcatPPi.*x4; % RHS x4
dxdt(5,:) = kprecip.*x3 - kdissolve.*x5; % RHS x5
dxdt(6,:) = kdegNTP.*x1; % RHS x6
dxdt(7,:) = kabortive.*x1; % RHS x7
dxdt(8,:) = koff.*x9 - kon.*x2.*x8; % RHS x8
dxdt(9,:) = kon.*x2.*x8 - koff.*x9; % RHS x9
