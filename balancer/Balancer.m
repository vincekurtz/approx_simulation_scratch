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
f_func = matlabFunction(f_sym, 'file', 'BalancerDynamics','vars',{x_sym,u_sym})


%% Simulation
x0 = [pi/2 0 0 0]';

anim_params.slow_mo = 1;
anim_params.frame_rate = 30;

control_params.test = 1;

[t_store, x_store , u_store] = Sim(x0, control_params, anim_params);
q_store = x_store(:,1:2)';

% Animate the results
showmotion(double_pendulum_model, t_store, q_store)

function [t_store, state_store , u_store] = Sim(x0, control_params, anim_params )
    x = x0;

    dt = 1e-3;

    state_store = [];
    t_store = [];
    u_store = [];

    for t = 0:dt:10
        u = [0 ; 0];

        t_store(end+1) = t;
        state_store(end+1,:) = x';
        u_store(end+1,:) = u;

        % Simple euler integration
        dx = BalancerDynamics(x, u);
        x = x+dx*dt;

    end

end


