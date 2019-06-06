%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use the spatial_v2 library and the associated visualization
% tools to simulate balencing the multilink balancer via
% planning with the LIP model.
%
% Assumes that create_multilink_spatial_model has already been
% run, generating a saved model and functions representing the
% balancer dynamics and relevant quantities (H,C,J,Jdot,etc).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Load the spatial_v2 model of the balancer
load('balancer_model')

% Linearized Centroid Dynamics
m = model.NB;  % total mass
A_lin = [zeros(2,3) 1/m*eye(2); zeros(3,5)];
B_lin = [zeros(2,3);eye(3)];
K = lqr(A_lin,B_lin,diag([10;10;1;1;1]),0.1*eye(3))

% Simulation parameters
T = 10;     % simulation time
dt = 1e-2;  % timestep

q0 = [0.2;0.2;0.2;0.4];   % initial state of full model
qd0 = [0.1;0;0;0];
x = [q0;qd0];

% Quantities to save for plotting
state_trajectory = [];
com_trajectory = [];

for t = 1:dt:T
    tic
    q = x(1:model.NB);    % unpack joint state
    qd = x(model.NB+1:end);  

    % Linearization
    A = [zeros(3,2) eye(3) zeros(3,1)]*A_com(q);    % select only nonzero terms
    Ad_qd = [zeros(3,2) eye(3) zeros(3,1)]*Ad_com_qd(q,qd);

    % Compute the state of the linearized system
    com_pos = p_com(q);
    com_h = A*qd;
    x_lin = [com_pos;com_h];

    % Compute a control to balance the feedback linearized system
    x_des = [0;1.5;0;0;0];
    u_com = -K*(x_lin-x_des);

    % Feedback linearization to generate torques as actual control
    Lambda = inv(A*inv(H(q))*A');
    u = A'*(Lambda*u_com - Lambda*Ad_qd + Lambda*A*inv(H(q))*C(q,qd));

    % Simulate the full model forward
    xdot = BalancerDynamics(x,u);
    x = x + xdot*dt;

    state_trajectory(:,end+1) = x;
    com_trajectory(:,end+1) = x_lin;
    
    toc
end

% Plot an animation of the balancer
q_trajectory = state_trajectory(1:model.NB,:);
showmotion(model, 1:dt:T, q_trajectory);

plot(com_trajectory');

% Show an animation of the LIP
addpath("../balancer")
%figure;
%animate_lip(1:dt:T, lip_trajectory', lip_control',h)
