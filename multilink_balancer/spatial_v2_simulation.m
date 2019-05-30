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

% Derive LIP and linearized CoM dynamics along with an interface
setup_interface

% Simulation parameters
T = 10;     % simulation time
dt = 1e-2;  % timestep

q0 = [0.9;0.2;0.2;0.4];   % initial state of full model
qd0 = zeros(4,1);
x = [q0;qd0];

% initial state of LIP model
com_pos = p_com(q0);
com_vel = pd_com(q0,qd0);
x_lip = [com_pos(1);com_vel(1);h;0];

% Quantities to save for plotting
state_trajectory = [];
com_trajectory = [];
lip_trajectory = [];
lip_control = [];

for t = 1:dt:T
    tic
    q = x(1:model.NB);    % unpack joint state
    qd = x(model.NB+1:end);  

    % Compute LIP controller that will let us balance
    u_lip = -K_lip*x_lip(1:2);

    % Compute the state of the linearized system (CoM x and y position and velocity)
    com_pos = p_com(q);
    com_vel = pd_com(q,qd);
    x_com = [com_pos(1);com_vel(1);com_pos(2);com_vel(2)];

    % Compute virtual control for the feedback linearized system
    u_com = R*u_lip + Q*x_lip + K_joint*(x_com-P*x_lip);

    % Feedback linearization to generate torques as actual control
    u = J(q)'*(Lambda(q)*u_com - Lambda(q)*Jdot(q,qd)*qd + Lambda(q)*J(q)*inv(H(q))*C(q,qd));

    % Compute null-space projector
    Jbar = inv(H(q))*J(q)'*Lambda(q);
    N = (eye(4) - J(q)'*Jbar')';

    % Apply secondary control
    kp = diag([0;0;0;1]);
    kd = diag([1;1;1;0]);
    q_des = [pi/2-0.1;0;0;0.5*sin(t);];
    qd_des = [0;0;0;0];
    u0 = -kp*(q-q_des) - kd*(qd-qd_des);  % secondary priority torques
    u = u + N'*u0;

    % Simulate the full model forward
    xdot = BalancerDynamics(x,u);
    x = x + xdot*dt;

    % Simulate the lip model forward
    xdot_lip = A2*x_lip + B2*u_lip;
    x_lip = x_lip + xdot_lip*dt;

    state_trajectory(:,end+1) = x;
    com_trajectory(:,end+1) = p_com(q);
    lip_trajectory(:,end+1) = x_lip;
    lip_control(:,end+1) = u_lip;
    
    toc
end

% Plot an animation of the balancer
q_trajectory = state_trajectory(1:model.NB,:);
showmotion(model, 1:dt:T, q_trajectory);

% Plot the CoM trajectory
plot(1:dt:T,com_trajectory)

% Show an animation of the LIP
addpath("../balancer")
%figure;
%animate_lip(1:dt:T, lip_trajectory', lip_control',h)
