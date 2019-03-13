%% Finding Bisimulation Functions for a Simple Linear System

clear all;
clc;

% We'll consider 2 linear systems which we know are approximately bisimilar
% both of which take the same input u.
A1 = [ 0  1  0; 
      -1 -1  1;
       0  0  0];
B1 = [0;1;1];

A2 = 0;
B2 = 1;

%% To double check this, we can simulate both forward in time 

% Simulation parameters
dt = 0.1;
Ns = 1000;

% Place to record the results
x1 = zeros(3,Ns);
x2 = zeros(1,Ns);

% Initial conditions
x1(:,1) = [0.5;0;0];
x2(:,1) = 0;
for i = 1:(Ns-1)
    u = 0.1*sin(0.03*i);  % control pre-specified
    
    % Simulate x1 forward in time
    x1_dot = A1*x1(:,i)+B1*u;
    x1(:,i+1) = x1(:,i) + x1_dot*dt;
    
    % Simulate x2 forward in time
    x2_dot = A2*x2(:,i)+B2*u;
    x2(:,i+1) = x2(:,i) + x2_dot*dt;
    
end

% Plot the resulting trajectories: x1 should follow x2
hold on
plot(x1(1,:))
plot(x2(1,:))
legend("x_1$","x_2")

%% Now on to the interesting part: calculating a Bisimulation funciton
% We'll try to find bisimulation function V(x1,x2) such that:
%  -V(x1,x2)+(g1(x1)-g2(x2))^2 is SOS
%  -dV/dx1*f1-dV/dx2*f2-lambda*V is SOS, for every u in the grid

% First we define the gridding over u
u_min = -0.1;
u_max = 0.1;
n_grids = 1;
u = linspace(u_min,u_max,n_grids);

% First we need to define the system sybmolically
x1 = sdpvar(3,1);
x2 = sdpvar(1);

% We generate f(x,u) for all u's in the grid
f1 = {};
f2 = {};
for i = 1:n_grids
    f1{i} = A1*x1+B1*u(i);  
    f2{i} = A2*x2+B2*u(i);
end

g1 = [1 0 0]*x1;
g2 = [1]*x2;

% Define the V sybmolically, where c are the constants used
[V,c,v] = polynomial([x1;x2],2);

% Now we define the partial of V w.r.t x1 and x2
nablaV_x1 = jacobian(V,x1);
nablaV_x2 = jacobian(V,x2);

% Choosing some value of lambda (?)
lambda = 1e-2;

% Finally, here are our SOS constraints
F = sos(V-(g1-g2)^2);
for i = 1:n_grids
    F = F + [sos(-nablaV_x1*f1{i}-nablaV_x2*f2{i})];
end

% Now run the magic!
[sol, u, Q] = solvesos(F,[],[],[c]);

if ~sol.problem
    disp("SOS Solution Found!")
    Vn = clean(u{1}'*Q{1}*u{1},1.e-6);
    sdisplay(Vn);
else
    disp("Constraints Failed. :(")
end
