% A demonstration of using an approximate simulation-based interface to balance
% a double-pendulum based on the Linear Inverted Pendulum model

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Definitions                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear Inverted Pendulum
h = 0.75;
g = 9.81;
omega = sqrt(g/h);

A_lip = [0 1; omega^2 0];
B_lip = [0; -omega^2];

K_lip = lqr(A_lip, B_lip, eye(2), 1);  % stabilizing feedback control gain

% Double Pendulum
syms theta1 theta2 theta1_dot theta2_dot real;
syms tau1 tau2 real
q_sym = [theta1;theta2];
qd_sym = [theta1_dot;theta2_dot];

x_sym = [q_sym;qd_sym];
u_sym = [tau1;tau2];

dp_model = autoTree(2);   % an unbranched kinematic tree with 2 bodies
dp_model.gravity = [0;-g;0];
qdd_sym = FDab( dp_model, q_sym, qd_sym, u_sym);
f_sym = [qd_sym; qdd_sym];
f_func = matlabFunction(f_sym,'vars',{x_sym, u_sym});

% Feedback Linearized Inverted Pendulum
xc_sym = 3/4*cos(theta1)+1/4*cos(theta1+theta2);   % COM position
yc_sym = 3/4*sin(theta1)+1/4*sin(theta1+theta2);
pc_sym = [xc_sym;yc_sym];
pc_func = matlabFunction(pc_sym, 'vars',{x_sym});

J_sym = jacobian(pc_sym, q_sym);   % Jacobian maps state to COM velocities
J = matlabFunction(J_sym, 'vars', {q_sym});

[H_sym, C_sym] = HandC(dp_model, q_sym, qd_sym);  % Parameters of the full system dynamics for feedback linearization
H = matlabFunction(H_sym, 'vars', {q_sym});
C = matlabFunction(C_sym, 'vars', {q_sym, qd_sym});

% TODO: figure out how to compute this symbolically
Jdot_sym =  [- (3*cos(theta1)*theta1_dot)/4 - (cos(theta1 + theta2)*(theta1_dot + theta2_dot))/4, -(cos(theta1 + theta2)*(theta1_dot + theta2_dot))/4;
            - (3*sin(theta1)*theta1_dot)/4 - (sin(theta1 + theta2)*(theta1_dot + theta2_dot))/4, -(sin(theta1 + theta2)*(theta1_dot + theta2_dot))/4];
Jdot = matlabFunction(Jdot_sym, 'vars', {q_sym, qd_sym});

% Convinience function that maps joint states to center of mass position
% and velocity, ie [x_c,x_c',y_c,y_c']
pc_dot_sym = J_sym*qd_sym;   % center of mass velocity
x_com_sym = [pc_sym(1);pc_dot_sym(1);pc_sym(2);pc_dot_sym(2)];
x_com = matlabFunction(x_com_sym, 'vars', {x_sym});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface Design                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
C1 = eye(4);     % the outputs of both systems are the full states

% The abstract system is the linear inverted pendulum
A2 = [0       1 0 0;
      omega^2 0 0 0;
      0       0 0 1;
      0       0 0 0];
B2 = [0        ;
      -omega^2 ;
      0        ;
      0        ];
C2 = eye(4);

% external control is position of CoP for the lip
u = sdpvar(1);   

% State of the concrete system
x1 = sdpvar(4,1);

% State of the abstract system
x2 = sdpvar(4,1);

% Interface definition
P = eye(4);  % P and Q chosen such that PA2 = A1P+B1*Q and C2 = C1*P
Q = [omega^2 0 0 0;
     0       0 0 0];
R = [-omega^2;0];
K_joint = -lqr(A1,B1,eye(4),eye(2));    % A control gain that stabilizes the concrete system
                                            % Could derive from a PD controller as well

lambda = 0.1;
Mbar = sdpvar(4,4);
Kbar = K_joint*Mbar;

F = [ [Mbar Mbar*C1'; C1*Mbar eye(4)] >= 0];
F = F + [Mbar*A1' + A1*Mbar + Kbar'*B1'+B1*Kbar + 2*lambda*Mbar <= 0 ];
optimize(F);

M = inv(value(Mbar));

% Double check the original conditions
check1 = all(eig(M-C1'*C1) >= 0);   % M >= C'C
check2 = all(eig((A1+B1*K_joint)'*M + M*(A1+B1*K_joint) + 2*lambda*M) <= 0);  % (A+BK)'M+M(A+BK) <= -2lambdaM

if check1 & check2
    disp("Checks passed!")
else
    disp("Checks Failed!")
    return;
end

% We'll replace u, x1, x2, with matlab symbolic variables so we can save the interface as a function
u = sym('u');
x1 = sym('x1',[4,1]);
x2 = sym('x2',[4,1]);
interface = R*u + Q*x2 + K_joint*(x1-P*x2);
interface_fcn = matlabFunction(interface,'vars',{u,x1,x2});

% Simulation function
V = (P*x2-x1)'*M*(P*x2-x1);
disp("V = (x1-x2)'*M*(x1-x2)")
disp("M: ")
M
disp("")

% Class K function of control input
gmma = norm(sqrt(M)*(B1*R-P*B2))/lambda;
disp("Gamma: ")
sdisplay(gmma)
disp("")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 10;   % length of simulation
dt = 1e-3;

% For storing data to plot later
t_sim = [];
lip_state = [];
lip_control = [];
dp_state = [];

% Torque limits
tau_min = [-100;-100];
tau_max = [100;100];

% Initial state of the full system
x = [0.9;0.1;0;0];         % initial joint states and velocities

% Initial state of the LIP model
com_init = x_com(x); 
x_lip = [com_init(1:2);h;0];  % we start the lip with same x value and velocity as the full system

for t = 1:dt:T

    % Compute the LIP control that will let us balance.
    % Only the x position and velocity are considered
    u_lip = -K_lip*x_lip(1:2);

    % Translate the LIP control to the full system control
    x1 = x_com(x);   % state of the full feedback linearized system
    x2 = x_lip;      % state of the abstract (LIP) system
    u_lin = interface_fcn(u_lip, x1, x2);   % control to apply to full feedback linearized system

    % Translate the linearized control to joint torques
    q = x(1:2);
    qd = x(3:4);
    tau = H(q)*inv(J(q))*(u_lin-Jdot(q,qd)*qd)+C(q,qd);
   
    % Apply the controls to the full system
    dx = BalancerDynamics(x, tau);
    x = x + dx*dt;
    
    % Apply the LIP control to the LIP model
    dx_lip = A2*x_lip + B2*u_lip;
    x_lip = x_lip + dx_lip*dt;

    % Record the corresponding values
    t_sim(end+1,:) = t;
    lip_state(end+1,:) = x_lip;
    lip_control(end+1,:) = u_lip;
    dp_state(end+1,:) = x;

end

% Calculate the actual CoM position for the recorded trajectory
com_state = x_com(dp_state')';

% Plot the errors along the trajectories
subplot(2,2,1);
hold on
plot(t_sim, lip_state(:,1)-com_state(:,1));
title("CoM X position");
xlabel("time")
ylabel("error")

subplot(2,2,2);
hold on
plot(t_sim, lip_state(:,3)-com_state(:,3));
title("CoM Y position");
xlabel("time")
ylabel("error")

subplot(2,2,3);
hold on
plot(t_sim, lip_state(:,2)-com_state(:,2));
title("CoM X velocity");
xlabel("time")
ylabel("error")

subplot(2,2,4);
hold on
plot(t_sim, lip_state(:,4)-com_state(:,4));
title("CoM Y velocity");
xlabel("time")
ylabel("error")

% Plot the net error and Bisimulation function
err = zeros(size(t_sim));
V_sim = zeros(size(t_sim));
for t = 1:length(t_sim)
    err(t) = (lip_state(t,:)-com_state(t,:))*(lip_state(t,:)-com_state(t,:))';
    V_sim(t) = (lip_state(t,:)-com_state(t,:))*M*(lip_state(t,:)-com_state(t,:))';
end

figure;
hold on
plot(t_sim, err);
plot(t_sim, V_sim);
legend("Output Error","Simulation Function");
xlabel("time")


% Play the animations
animate_lip_dp(t_sim, dp_state, lip_state(:,1:2), lip_control, h);

