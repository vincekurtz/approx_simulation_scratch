%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use the spatial_v2 library and the associated visualization
% tools to simulate the multilink balancer. 
%
% Assumes that create_multilink_spatial_model has already been
% run, generating a saved model and functions representing the
% balancer dynamics and relevant quantities (H,C,J,Jdot,etc).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Load the model
load('balancer_model')

% Simulation parameters
T = 10;     % simulation time
dt = 1e-2;  % timestep

q0 = [0.2;0.2;0.2;0.4];
qd0 = zeros(4,1);
x = [q0;qd0];   % initial state

% Simulation
state_trajectory = [];
com_trajectory = [];

for t = 1:dt:T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % feedback controller
    tic
    q = x(1:model.NB);    % unpack joint state
    qd = x(model.NB+1:end);  

    x_lin = [p_com(q);J(q)*qd];
    x_des = [0;1.5;0;0];

    K = [1 0 1.7 0; 0 1 0 1.7];
    v = -K*(x_lin-x_des);       % virtual control for linearized system

    % Feedback linearization to generate torques as actual control
    u = J(q)'*(Lambda(q)*v - Lambda(q)*Jdot(q,qd)*qd + Lambda(q)*J(q)*inv(H(q))*C(q,qd));
    toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xdot = BalancerDynamics(x,u);
    x = x + xdot*dt;

    state_trajectory(:,end+1) = x;
    com_trajectory(:,end+1) = p_com(q);
end

% Plot an animation of the motion
q_trajectory = state_trajectory(1:model.NB,:);
showmotion(model, 1:dt:T, q_trajectory);

plot(1:dt:T,com_trajectory)
