% Computing a simulation function for two linear systems: one is the linear
% inverted pendulum (LIP) model, and the other is a feedback-linearized double pendulum

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
B2 = [1        ;
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

%% Finding a Bisimulation function and interface, following Girard and Pappas '09

% Feedback gain that stabilizes the concrete system, from LQR or similar
K = -[1 1.7321 0 0;
     0 0      1 1.7321;];

lambda = sdpvar(1);
Mbar = sdpvar(4,4);
Kbar = K*Mbar;

F = [ [Mbar Mbar*C1'; C1*Mbar eye(2)] >= 0];
F = F + [Mbar*A1' + A1*Mbar + Kbar'*B1'+B1*Kbar + 2*lambda*Mbar <= 0 ];
optimize(F,-lambda);  % maximize lambda
lambda = value(lambda)

M = inv(value(Mbar));

% Double check the original conditions
check1 = all(eig(M-C1'*C1) >= 0);   % M >= C'C
check2 = all(eig((A1+B1*K)'*M + M*(A1+B1*K) + 2*lambda*M) <= 0);  % (A+BK)'M+M(A+BK) <= -2lambdaM

if check1 & check2
    disp("Checks passed!")
else
    disp("Checks Failed!")
    return;
end

% Interface definition
P = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0];  % P and Q chosen such that PA2 = A1P+B1*Q and C2 = C1*P
Q = [omega^2 0 0 0;
     0       0 0 0];
R = [1;0];

interface = R*u + Q*x2 + K*(x1-P*x2);

% Simulation function
V = (P*x2-x1)'*M*(P*x2-x1);

% Class K gamma function of control input
gmma = norm(sqrt(M)*(B1*R-P*B2))/lambda;
disp("Gamma: ")
sdisplay(gmma)
disp("")

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
%figure;
%hold on;
%plot(x1_sim(3,:))
%plot(x2_sim(3,:))
%xlabel('time')
%ylabel('y position')
%legend("y_1: concrete","y_2: abstract")
%title("System Trajectories")

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
