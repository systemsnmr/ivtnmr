function y = OBSxIVT(t,x,p)

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

y = [
x3
x1
x2 + x9
x7
x6
x5
x9
x9
];
