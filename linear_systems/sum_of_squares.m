% Computing a simulation function for two linear systems based on formulating
% the problem as an Input-to-Output stability problem

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


%% Defining the joint system and SOS constraints
u = sdpvar(1);

x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

x = [x1;x2;x3];

% System definition
[interface,ic] = polynomial([x;u],1);
%interface = u-x1-x2+x3;  % an interface that we know works decently

f = [A1*[x1;x2] + B1*interface;
     A2*x3 + B2*u];
 
% System outputs
g = x1-x3;

% Candidate ISS Lyapunov function
[V,vc] = polynomial(x,2);
Vdot = jacobian(V,x)*f;

% Class K_inf functions should be even polynomials
alpha1 = g^2;     % enforcing this means that V(x) bounds the output error, without losing much generality
[alpha2,a2c,da2] = even_polynomial(x,2);
[alpha3,a3c,da3] = even_polynomial(x,2);
[alpha4,a4c,da4] = even_polynomial(u,2);

% To help with numerical issues
epsilon = 1e-6;

% SOS constraints
F = [sos(V-alpha1-epsilon*x'*x)];
F = F + [sos(alpha2-V-epsilon*x'*x)];
F = F + [sos(-alpha3+alpha4-Vdot)];
F = F + [sos(da2),sos(da3),sos(da4)];

% It helps BMI solvers to have explicit bounds on variables
F = F + [ 100 >= [vc;a2c;a3c;a4c] >= -100];
F = F + [ 1 >= ic >= -1];

% Options
ops = sdpsettings('penbmi.PBM_MAX_ITER',1000, ...
                  'penbmi.NWT_SYS_MODE',0,   ...
                  'penbmi.LS',1,             ...
                  'penbmi.UM_MAX_ITER',200);
              
% ops = sdpsettings('solver','bmibnb',        ...
%                   'debug',1,                ...
%                   'bmibnb.lpreduce',1,      ...
%                   'sos.model',2,            ...
%                   'bmibnb.roottight',0);

[sol] = solvesos(F,[],ops,[vc;a2c;a3c;a4c;ic]);  % minimizing alpha4 makes V a tigher bound

if ~sol.problem
    disp("Solution Found!")
    disp("V: ")
    sdisplay(replace(V,vc,value(vc)),[3])
    disp("interface: ")
    sdisplay(replace(interface,ic,value(ic)),[3])
else
    disp("No solution found :(")
end

%% Simulating both systems to check the results

% Simulation parameters
dt = 0.1;
Ns = 500;
x0 = [0.2;0;0];   % initial conditions

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
    u1 = replace(interface,[u;x1;x2;x3;ic],[u2(:,i);x1_sim(:,i);x2_sim(:,i);value(ic)]);
    x1_dot = A1*x1_sim(:,i)+B1*u1;
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
xlabel('time')
ylabel('position')
legend("x_1: concrete","x_2: abstract")
title("System Trajectories")

% Plot some of the quantities we've calculated to certify simulation
V_sim = zeros(1,Ns);
Vdot_sim = zeros(1,Ns);
alpha1_sim = zeros(1,Ns);
alpha2_sim = zeros(1,Ns);
alpha3_sim = zeros(1,Ns);
alpha4_sim = zeros(1,Ns);
for i = 1:Ns
    V_sim(:,i)= replace(V,[vc;x1;x2;x3;ic],[value(vc);x1_sim(:,i);x2_sim(:,i);value(ic)]);
    Vdot_sim(:,i)= replace(Vdot,[vc;x1;x2;x3;u;ic],[value(vc);x1_sim(:,i);x2_sim(:,i);u2(:,i);value(ic)]);
    alpha1_sim(:,i) = replace(alpha1,[x1;x2;x3],[x1_sim(:,i);x2_sim(:,i)]);
    alpha2_sim(:,i) = replace(alpha2,[a2c;x1;x2;x3],[value(a2c);x1_sim(:,i);x2_sim(:,i)]);
    alpha3_sim(:,i) = replace(alpha3,[a3c;x1;x2;x3],[value(a3c);x1_sim(:,i);x2_sim(:,i)]);
    alpha4_sim(:,i) = replace(alpha4,[a4c;x1;x2;x3;u],[value(a4c);x1_sim(:,i);x2_sim(:,i);u2(:,i)]);
end

figure;
hold on;
plot(V_sim)
plot(alpha1_sim)
legend("V(x)","\alpha_1(|g(x)|) = |g(x)|^2")
title("Error Bound: \alpha_1(|g(x)|) \leq V(x)")
xlabel('time')

function [p, c, sdp_ds] = even_polynomial(s, n)
    % Helper function to construct a symbolic even polynomial
    %
    % Arguments:
    %    s - a YALMIP symbolic variable
    %    n - the number of coefficients
    %
    % Returns:
    %    p      - a symbolic even polynomial of the form p = c1*s^2+c2s^4+...
    %    c      - the associated constants
    %    sdp_ds - the polynomial s*dp(s)/ds, which forms a sos constraint
    %             for p to be a class K_inf function

    c = sdpvar(n,1);
    raw_poly = [s'*s];
    for i=2:n
        raw_poly(i,:) = raw_poly(1,:)^(i);
    end
    p = c'*raw_poly;

    cds = [2:2:2*n]'.*c;    % coefficients of s*dp(s)/ds
    sdp_ds = cds'*raw_poly;

end

