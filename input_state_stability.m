% We can formulate the simulation problem as Input-Output stability of the
% joint system, so here we'll consider starting from an ISS perspective

close all;
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
x = [x1;x2;x3];
x0 = [0.2;0;0];   % initial conditions

u = sdpvar(1);

% System definition
f = [A1*[x1;x2] + B1*[u+K*([x1;x2]-P*x3)];
     A2*x3 + B2*u];
 
% System outputs
g = x1-x3;

% Candidate ISS Lyapunov function
[V,vc] = polynomial([x1;x2;x3],2);
Vdot = jacobian(V,[x1;x2;x3])*f;

% We'll try to minimize this
V0 = replace(V,x,x0);

% For ensuring PD-ness
epsilon = 1e-8;

% Class K/K_inf functions that should bound Vdot
[alpha,ac] = polynomial([x1;x2;x3],3);
[gamma,gc] = polynomial([u],3);
[alpha_,a_c] = polynomial([g],2);
[alpha_overbar,ao_c] = polynomial([x1;x2;x3],2);

% SOS constraints
F = [sos(V-alpha_)];
F = F + [sos(alpha_overbar-V)];
F = F + [sos(alpha),sos(gamma),sos(alpha_),sos(alpha_overbar)];
F = F + [sos(-alpha+gamma-Vdot)];

[sol,u,Q,res] = solvesos(F,[],[],[vc;ac;gc;a_c;ao_c]);

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
x1_sim = zeros(2,Ns);
x2_sim = zeros(1,Ns);
u2 = zeros(1,Ns);

% Initial conditions
x1_sim(:,1) = x0(1:2);
x2_sim(:,1) = x0(3);
for i = 1:(Ns-1)
    u2(:,i) = 0.1*sin(0.03*i);  % control pre-specified
    
    if i>300
        u2(:,i) = 0;
    end
    
    % Simulate x1 forward in time
    interface = u2(:,i) + K*(x1_sim(:,i)-P*x2_sim(:,i));
    x1_dot = A1*x1_sim(:,i)+B1*interface;
    x1_sim(:,i+1) = x1_sim(:,i) + x1_dot*dt;
    
    % Simulate x2 forward in time
    x2_dot = A2*x2_sim(:,i)+B2*u2(:,i);
    x2_sim(:,i+1) = x2_sim(:,i) + x2_dot*dt;
    
end

% Plot x1 vs x2
figure;
hold on;
plot(x1_sim(1,:))
plot(x2_sim(1,:))
legend("x_1","x_2")

% Plot V(x1,x2) vs alpha_(|g1-g2|)
figure;
hold on;

V_sim = zeros(1,Ns);
alpha_sim = zeros(1,Ns);
for i = 1:Ns
    V_sim(:,i)= replace(V,[vc;x1;x2;x3],[value(vc);x1_sim(:,i);x2_sim(:,i)]);
    alpha_sim(:,i) = replace(alpha_,[a_c;x1;x2;x3],[value(a_c);x1_sim(:,i);x2_sim(:,i)]);
end

plot(V_sim)
plot(alpha_sim)
error = x1_sim(1,:)-x2_sim(1,:);
plot(error)
legend("V","alpha_","error")
