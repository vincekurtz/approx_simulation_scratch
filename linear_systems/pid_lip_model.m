% Computing a simulation function for two linear systems: one is the linear
% inverted pendulum (LIP) model, and the other is a feedback-linearized double pendulum.
%
% Here we'll use a simple PD controller to track the LIP model and use this interface
% to certify approximate simulation.

close all;
clear all;
clc;

%% Original two systems
% Concrete system is a double integrator in x and y position, which
% is obtained from task-space linearization of the double pendulum
A1 = [0 1 0 0;
      0 0 0 0;
      0 0 0 1;
      0 0 0 0];
B1 = [0 0;
      1 0;
      0 0;
      0 1];
C1 = [1 0 0 0;     % the outputs of both systems are the x and y position only
      0 0 1 0];

% The abstract system is the linear inverted pendulum
omega = 1;
A2 = [0       1 0 0;
      omega^2 0 0 0;
      0       0 0 0;
      0       0 0 0];
B2 = [0        ;
      -omega^2 ;
      0        ;
      0        ];
C2 = [1 0 0 0;
      0 0 1 0];

% external control is position of CoP for the lip
u = sdpvar(1);   

% State of the concrete system
x1 = sdpvar(4,1);

% State of the abstract system
x2 = sdpvar(4,1);

%% Defining an interface as a PD controller that tracks the abstract system

p       = [x1(1);x1(3)];                % Actual COM position       
pd      = [x1(2);x1(4)];                % Actual COM velocity
p_des   = [x2(1);x2(3)];                % Desired CoM position
pd_des  = [x2(2);x2(4)];                % Desired CoM velocity
pdd_des = [omega^2*(x2(1)-u); 0];       % Desired Center of mass accelration

Kd = [1.5 0; 0 0.9];
Kp = [10 0; 0 0.9];
%Kd = 0; K=0;

interface = pdd_des + Kd*(pd_des-pd) + Kp*(p_des-p);

% DEBUG
P = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0];  % P and Q chosen such that PA2 = A1P+B1*Q and C2 = C1*P
Q = [omega^2 0 0 0;
     0       0 0 0];
R = [1;0];
K = -[5 2 0 0;
     0 0      5 2;];
interface = R*u + Q*x2 + K*(x1-P*x2);

%% Finding a Bisimulation function

% Joint system definition
x = [x1;x2];

f = [A1*x1 + B1*interface;
     A2*x2 + B2*u];
 
g = C1*x1 - C2*x2;

% Candidate ISS Lyapunov function
[V,vc] = polynomial(x,2);
Vdot = jacobian(V,x)*f;

% Class K_inf functions should be even polynomials
%alpha1 = g'*g;     % enforcing this means that V(x) bounds the output error, without losing much generality
[alpha1,a1c,da1] = even_polynomial(g,2);
[alpha2,a2c,da2] = even_polynomial(x,2);
[alpha3,a3c,da3] = even_polynomial(x,2);
[alpha4,a4c,da4] = even_polynomial(u,2);

% To help with numerical issues
epsilon = 1e-6;
epsilon = 0;

% SOS constraints
F = [sos(V-alpha1-epsilon*x'*x)];
F = F + [sos(alpha2-V-epsilon*x'*x)];
F = F + [sos(-alpha3+alpha4-Vdot-epsilon*x'*x)];
F = F + [sos(da1),sos(da2),sos(da3),sos(da4)];

ops = sdpsettings('sos.model',2);

[sol] = solvesos(F,[],ops,[vc;a1c;a2c;a3c;a4c;]);  % minimizing alpha4 makes V a tigher bound

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
x0 = [0;0;1;0;
      0;0;1;0];   % initial conditions

% Place to record the results
x1_sim = zeros(4,Ns);
x2_sim = zeros(4,Ns);
u2 = zeros(1,Ns);

% Initial conditions
x1_sim(:,1) = x0(1:4);
x2_sim(:,1) = x0(5:8);
for i = 1:(Ns-1)
    u2(:,i) = 0.1*sin(0.03*i);  % control pre-specified
    
    %if i>300
    %    u2(:,i) = 0;
    %end
    
    % Simulate x1 forward in time
    u1 = replace(interface,[u;x1;x2;],[u2(:,i);x1_sim(:,i);x2_sim(:,i)]);
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
ylabel('x position')
legend("x_1: concrete","x_2: abstract")
title("System Trajectories")


% Plot y1 vs y2
figure;
hold on;
plot(x1_sim(3,:))
plot(x2_sim(3,:))
xlabel('time')
ylabel('y position')
legend("y_1: concrete","y_2: abstract")
title("System Trajectories")

return

% Plot the output error and error bound
V_init = replace(V,[x1;x2],[x1_sim(:,1);x2_sim(:,1)]);
gamma_max = gmma*norm(max(u2));
error_bound = max(V_init, gamma_max);   % by Theorem 1

err = vecnorm(C1*x1_sim-C2*x2_sim);
figure;
hold on;
plot(err);
yline(error_bound,'color','red');
legend("Actual Error","Error Bound")
title("Error Bound vs Actual Error")
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
