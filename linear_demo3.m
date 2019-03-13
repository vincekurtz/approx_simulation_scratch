%% An attempt an finding a Barrier Certificate with our own SOS formulation

clear all;
clc;

%% System Definition
% We'll consider 2 linear systems which we know are approximately bisimilar
% with a given interface
u1 = sdpvar(1);
u2 = sdpvar(1);
x1 = sdpvar(2,1);
x2 = sdpvar(1);

A1 = [0 1; 
      0 0 ];
B1 = [0;1];

A2 = 0;
B2 = 1;

% Interface
K = [-1 -1];
P = [1;0];
interface = u2+K*(x1-P*x2);

f1 = A1*x1+B1*interface;  
f2 = A2*x2+B2*u2;

g1 = [1 0]*x1;
g2 = [1]*x2;


%% Finding a solution
% Now we'll try to find
%   Bisimulation function V(x1,x2)
%   Class K function gamma(|u|)
%   SOS function lambda(u,x1,x2)
% Such that
%   V(x1,x2)-|g1(x1)-g2(x2)| is SOS
%   lambda(u2,x1,x2)          is SOS
%  -dV/dx1*f1 - dV/dx2*f2 +lambda(u2,x1,x2)(gamma(u2)-V) is SOS

% Define V
[V,c] = polynomial([x1;x2],2);

% Define lambda
[lambda,clambda] = polynomial([u2;x1;x2],8);

% Define gamma
gamma = 10000*abs(u2);

% Now we define the partial of V w.r.t x1 and x2
nablaV_x1 = jacobian(V,x1);
nablaV_x2 = jacobian(V,x2);

% Add the SOS constraints
%F = [sos(V-(g1-g2)^2)];
%F = F + [sos(lambda)];
F = [sos( -nablaV_x1*f1-nablaV_x2*f2+lambda*(gamma-V))];  % even this constraint alone is infeasible

% Now solve!
[sol, u, Q] = solvesos(F,[],[],[c;clambda]);

if ~sol.problem
    disp("SOS Solution Found!")
    Vn = clean(u{1}'*Q{1}*u{1},1.e-6);
    sdisplay(Vn);
else
    disp("Constraints Failed. :(")
end