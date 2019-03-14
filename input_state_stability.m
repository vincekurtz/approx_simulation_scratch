% We can formulate the simulation problem as Input-Output stability of the
% joint system, so here we'll consider starting from an ISS perspective

clear all;
clc;

%% Original two systems and interface
% Systems in question are a double integrator
A1 = [0 1;
      0 0];
B1 = [0;1];
C1 = [1 0];

% And a single integrator
A2 = 0;
B2 = 1;
C2 = 1;

% An interface that maps the input of system 2 to an input of system 1
% is u2 + K*(x1-P*x2), where
K = [-1, -1];
P = [1;0];

%% ISS of joint system system

x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

u = sdpvar(1);

% System definition
f = [A1*[x1;x2] + B1*[u+K*([x1;x2]-P*x3)];
     A2*x3 + B2*u];
 
% System outputs
g1 = x1;
g2 = x3;
g = x1-x3;

% Candidate ISS Lyapunov function
[V,vc] = polynomial([x1;x2;x3],3);
Vdot = jacobian(V,[x1;x2;x3])*f;

% Class K/K_inf functions that should bound Vdot
[alpha,ac] = polynomial([x1;x2;x3],3);
[gamma,gc] = polynomial([u],4);
[alpha_,a_c] = polynomial([g],3);

% SOS constraints
F = [sos(V-alpha_)];
F = F + [sos(alpha),sos(gamma)];
F = F + [sos(-alpha+gamma-Vdot)];

[sol,u,Q,res] = solvesos(F,[],[],[vc;ac;gc;a_c]);

if ~sol.problem
    disp("Solution Found!")
    disp("V: ")
    sdisplay(replace(V,vc,value(vc)),[3])
else
    disp("No solution found :(")
end

%% Simulating both systems to check the results

% Simulation parameters
dt = 0.1;
Ns = 500;

% Place to record the results
x1 = zeros(2,Ns);
x2 = zeros(1,Ns);
u2 = zeros(1,Ns);

% Initial conditions
x1(:,1) = [0.2;0];
x2(:,1) = 0;
for i = 1:(Ns-1)
    u2(:,i) = 0.1*sin(0.03*i);  % control pre-specified
    
    if i>300
        u2(:,i) = 0;
    end
    
    % Simulate x1 forward in time
    interface = u2(:,i) + K*(x1(:,i)-P*x2(:,i));
    x1_dot = A1*x1(:,i)+B1*interface;
    x1(:,i+1) = x1(:,i) + x1_dot*dt;
    
    % Simulate x2 forward in time
    x2_dot = A2*x2(:,i)+B2*u2(:,i);
    x2(:,i+1) = x2(:,i) + x2_dot*dt;
    
end

% Plot x1 vs x2
figure;
hold on;
plot(x1(1,:))
plot(x2(1,:))
legend("x_1","x_2")

% Plot error
figure;
hold on;
error = x1(1,:)-x2(1,:);
plot(error)