%% Finding Bisimulation Functions for a Simple Linear System with SOSP2

clear all;
clc;

%% System Definition
% We'll consider 2 linear systems which we know are approximately bisimilar
% both of which take the same input u.
A1 = [ 0  1  0; 
      -1 -1  1;
       0  0  0];
B1 = [0;1;1];

A2 = 0;
B2 = 1;

% First we need to define the system sybmolically
u = sdpvar(1);
x1 = sdpvar(3,1);
x2 = sdpvar(1);

f1 = A1*x1+B1*u;  
f2 = A2*x2+B2*u;

g1 = [1 0 0]*x1;
g2 = [1]*x2;

% Bounding functions on control and states
u_min = -1;
u_max = 1;
rho = (u-u_min)*(u_max-u);

x_min = -1;   % for now we'll consider the same bounds for all state vars.
x_max = 1;
tao1 = [(x1(1)-x_min)*(x_max-x1(1)) ;
        (x1(2)-x_min)*(x_max-x1(2)) ;
        (x1(3)-x_min)*(x_max-x1(3)) ];
tao2 = (x2-x_min)*(x_max-x2);


%% Finding a solution
% Now we'll try to find bisimulation function V(x1,x2) such that:
%  -V(x1,x2)+(g1(x1)+g2(x2))^2 is SOS
%   sigma1(x1,u) is SOS
%   sigma2(x2,u) is SOS
%   sigma3(x1,u) is SOS
%   sigma4(x2,u) is SOS
%  -dV/dx1*f1 - dV/dx2*f2 - lambda*V - sigma1*rho 
%                 - sigma2*rho - sigma3*tao1 - sigma4*tao2 is SOS

% Define V
[V,c,v] = polynomial([x1;x2],2);

% Define sigma*
[sigma1,cs1,vsigma1] = polynomial([x1;u],2);
[sigma2,cs2,vsigma2] = polynomial([x2;u],2);
[sigma3,cs3,vsigma3] = polynomial([x1;u],2);
[sigma4,cs4,vsigma4] = polynomial([x2;u],2);

% Now we define the partial of V w.r.t x1 and x2
nablaV_x1 = jacobian(V,x1);
nablaV_x2 = jacobian(V,x2);

% Choosing some value of lambda (?)
lambda = 1e-2;

% Add the SOS constraints
F = [sos(V-(g1-g2)^2)];
F = F + [sos(sigma1)] + [sos(sigma2)] + [sos(sigma3)] + [sos(sigma4)];
F = F + [sos( -nablaV_x1*f1-nablaV_x2*f2-lambda*V-sigma1*rho-sigma2*rho-sigma3*tao1-sigma4*tao2)];

% Now solve!
[sol, u, Q] = solvesos(F,[],[],[c;cs1;cs2;cs3;cs4]);

if ~sol.problem
    disp("SOS Solution Found!")
    Vn = clean(u{1}'*Q{1}*u{1},1.e-6);
    sdisplay(Vn);
else
    disp("Constraints Failed. :(")
end
