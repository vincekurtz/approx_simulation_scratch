% We can formulate the simulation problem as Input-Output stability of the
% joint system, so here we'll consider starting from an ISS perspective

clear all;
clc;

%% ISS of example system
x1 = sdpvar(1);
x2 = sdpvar(1);
u1 = sdpvar(1);
u2 = sdpvar(1);

% System definition
f = [-x2-1.5*x1^2-0.3*x1^3;
     4*(x1+u1)-4*(x2+u2)-4.5*(x1+u1)^2-1.5*(x1+u1)^3];

% Candidate ISS Lyapunov function
[V,vc] = polynomial([x1;x2],3);
Vdot = jacobian(V,[x1;x2])*f;

% Class K/K_inf functions that should bound Vdot
[alpha,ac] = polynomial([x1;x2],3);
[gamma,gc] = polynomial([u1;u2],4);

% SOS constraints
F = [sos(V)];
F = F + [sos(alpha),sos(gamma)];
F = F + [sos(-alpha+gamma-Vdot)];

[sol,u,Q,res] = solvesos(F,[],[],[vc;ac;gc]);

if ~sol.problem
    disp("Solution Found!")
    disp(' ')
    disp("alpha: ")
    sdisplay(replace(alpha,ac,value(ac)),[3])
    disp(' ')
    disp("gamma: ")
    sdisplay(replace(gamma,gc,value(gc)),[3])
    disp(' ')
    disp("V: ")
    sdisplay(replace(V,vc,value(vc)),[3])
else
    disp("No solution found :(")
end
