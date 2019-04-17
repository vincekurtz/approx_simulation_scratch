clear;
clc;

%% Setup

% State variables
syms theta1 theta2 real
syms theta1_dot theta2_dot real
q_sym =[theta1 theta2]';
qd_sym = [theta1_dot theta2_dot]';
x_sym = [q_sym; qd_sym];

% Control Variables
syms tau1 tau2 real
tau_sym = [tau1 tau2]';
u_sym = tau_sym;

double_pendulum_model = autoTree(2);  % produces an unbranched kinematic tree with 2 bodies 
double_pendulum_model.gravity = [0;-9.81;0];   % point gravity in the -y direction

% Camera settings for showmotion
double_pendulum_model.camera.direction = [0 0 1];
double_pendulum_model.camera.up = [0 1 0];
double_pendulum_model.camera.zoom = 0.5;
double_pendulum_model.camera.locus = [0 0];

% State-space representation
qdd_sym = FDab( double_pendulum_model, q_sym, qd_sym, tau_sym);
f_sym = [qd_sym; qdd_sym];

% Generates and saves BalancerDynamics.m
f_func = matlabFunction(f_sym, 'file', 'BalancerDynamics','vars',{x_sym,u_sym});

%% Feedback Linearization

% End effector position
xe_sym = cos(theta1)+cos(theta1+theta2);
ye_sym = sin(theta1)+sin(theta1+theta2);
p_sym = [xe_sym;ye_sym];

% End effector jacobian
J_sym = jacobian(p_sym, q_sym);

% System Equations of motion in the form H(q)qdd+C(q,qd,f_ext)=tau,
% where H is the joint-space inertia matrix and C is the vector of gravity,
% external force, and correolis terms. 
[H_sym,C_sym] = HandC(double_pendulum_model, q_sym, qd_sym);


%% Simulation
x0 = [pi/2 0.01 0 0]';

control_params.x_nom = 0;   % nominal x and y that we'll control to
control_params.y_nom = -1;

[t_store, q_sim , u_store] = Sim(x0, control_params);
q_store = q_sim(:,1:2)';

% Animate the results
showmotion(double_pendulum_model, t_store, q_store)

function [t_store, state_store , u_store] = Sim(x0, control_params)
    x = x0;

    dt = 1e-3;
    
    umin = [-50;-50];
    umax = [50;50];

    state_store = [];
    t_store = [];
    u_store = [];

    for t = 0:dt:10
        u = control_law(x, control_params);
        
        % Set control limits
        u = max(umin,u);
        u = min(umax,u);

        t_store(end+1) = t;
        state_store(end+1,:) = x';
        u_store(end+1,:) = u;

        % Simple euler integration
        dx = BalancerDynamics(x, u);
        x = x+dx*dt;

    end
end

function u = control_law(x, control_params)
    % Note: state is x = [q;qd] = [theta1;theta2;theta1';theta2']
    theta1 = x(1);
    theta2 = x(2);
    theta1_dot = x(3);
    theta2_dot = x(4);
    qd = [theta1_dot;theta2_dot];
    
    % Calculate Jacobian and associated time derivative
    J = [-sin(theta1)-sin(theta1+theta2) , -sin(theta1+theta2);
         cos(theta1)+cos(theta1+theta2) ,  cos(theta1+theta2)];
    
    Jdot = [-theta1_dot*cos(theta1)-(theta1_dot+theta2_dot)*cos(theta1+theta2) , -(theta1_dot+theta2_dot)*cos(theta1+theta2);
            -theta1_dot*sin(theta1)-(theta1_dot+theta2_dot)*sin(theta1+theta2) , -(theta1_dot+theta2_dot)*sin(theta1+theta2)];
    
    % Calculate nominal control for linearized system
    xe_actual = cos(theta1)+cos(theta1+theta2);  % actual end effector position
    ye_actual = sin(theta1)+sin(theta1+theta2);
    
    p = [xe_actual;ye_actual];   % end effector position
    pdot = J*x(3:4);             % end effector velocity
    
    x_lin = [p;pdot];   % linearized state is [x;y;x';y']
    x_lin_err = [p;pdot] - [control_params.x_nom;control_params.y_nom;0;0];
    
    K = [1.0 0 1.7321 0; 0 1.0 0 1.7321];   % from LQR of linearized system
    v = -K*x_lin_err;    % nominal control input
    
    % Feedback linearized control input: we could also use subs() here to
    % compute H, C, etc, but that turns out to be super slow
    H = [cos(theta2)/2 + sin(theta2)^2 + cos(theta2)*(cos(theta2) + 1/2) + 803/1200, cos(theta2)/2 + 803/2400;
                                                           cos(theta2)/2 + 803/2400,                 803/2400];
                                                       
    C = [(981*cos(theta1 + theta2))/200 + (2943*cos(theta1))/200 - (theta2_dot^2*sin(theta2))/2 - theta1_dot*theta2_dot*sin(theta2);
                                                                      (sin(theta2)*theta1_dot^2)/2 + (981*cos(theta1 + theta2))/200];
         
    tau = H*inv(J)*(v-Jdot*qd)+C;
    
    u = tau;
end
