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

% Derive LIP and Linearized CoM dynamics along with an interface that ensures
% approximate simulation.
setup_interface

% Simulation parameters
T = 5;     % simulation time
dt = 5e-2;  % timestep

% Initial state of full model
q0 = [pi/4;pi/2;-pi/4;-pi/2];   
qd0 = [0.1;0.0;0.0;0.0];
x = [q0;qd0];

% Initial state of LIP model and linearized system
com_pos = p_com(q0);
com_h = A_com(q0)*qd0;
x_com = [com_pos;com_h];
x_lip = [com_pos(1);h;0;com_h(2);0];  % y position fixed at h, angular momentum fixed at 0.

% Define tracking method: AS (ours) or QP (traditional)
tracking_method = "QP";

if tracking_method == "AS"
    % Set up contact constraints
    params.A_lip = A2;
    params.B_lip = B2;
    params.A_com = A1;
    params.B_com = B1;
    params.N = 5;     % MPC horizon
    params.dt = dt;
    params.R = R;
    params.Q = Q;
    params.K = K_joint;
    params.m = m;
    constraint_params = LIPConstraints(params);
end

% Quantities to save for plotting
state_trajectory = [];
com_trajectory = [];
lip_trajectory = [];
lip_control = [];
err_sim = [];
V_sim = [];

for t = 1:dt:T
    tic
    q = x(1:model.NB);    % unpack joint state
    qd = x(model.NB+1:end);  

    % Compute the state of the linearized system
    com_pos = p_com(q);
    com_h = A_com(q)*qd;
    x_com = [com_pos;com_h];

    if tracking_method == "AS"
        % Compute a lip control
        u_lip_trajectory = GenerateLIPTrajectory(x_lip, x_com, constraint_params);
        u_lip = u_lip_trajectory(1);

        % Compute a control to track the LIP trajectory
        u_com = R*u_lip + Q*x_lip + K_joint*(x_com-P*x_lip);

        % Feedback linearization to generate torques as actual control
        Lambda = inv(A_com(q)*inv(H(q))*A_com(q)');
        u = A_com(q)'*(Lambda*u_com - Lambda*Ad_com_qd(q,qd) + Lambda*A_com(q)*inv(H(q))*C(q,qd));

        % Secondary control via null space projector
        Abar = inv(H(q))*A_com(q)'*Lambda;
        N = (eye(4) - A_com(q)'*Abar')';   % null space projector

        kp = diag([0;0;1;0.1]);    % we'll use PD control of joints to
        kd = diag([2;2;1;1]);    % apply commands in the null space
        q_des = [0;pi/4;-pi/2;-pi/4];
        qd_des = [0;0;0;0];
        u0 = -kp*(q-q_des)-kd*(qd-qd_des);
        u = u + N'*u0;

    elseif tracking_method == "QP"
        % Compute a lip control with LQR
        u_lip = -K_lip*[x_lip(1);x_lip(4)];
        
        % Compute a control to track the LIP trajectory
        u = QPTracker(u_lip, x_lip, q, qd, omega);
        
        %u = J(q)'*(Lambda(q)*u_com - Lambda(q)*Jdot(q,qd)*qd + Lambda(q)*J(q)*inv(H(q))*C(q,qd));

    end

    % Simulate the full model forward
    xdot = BalancerDynamics(x,u);
    x = x + xdot*dt;

    % Simulate the lip model forward
    xdot_lip = A2*x_lip + B2*u_lip;
    x_lip = x_lip + xdot_lip*dt;

    state_trajectory(:,end+1) = x;
    com_trajectory(:,end+1) = x_com;
    lip_trajectory(:,end+1) = x_lip;
    lip_control(:,end+1) = u_lip;
    err_sim(:,end+1) = sqrt((x_com-x_lip)'*(x_com-x_lip));
    V_sim(:,end+1) = sqrt((x_com-P*x_lip)'*M*(x_com-P*x_lip));
    
    toc
end

% Plot an animation of the balancer
q_trajectory = state_trajectory(1:model.NB,:);
showmotion(model, 1:dt:T, q_trajectory);

% Plot the error and simulation function
figure;
hold on
plot(1:dt:T, err_sim)
plot(1:dt:T, V_sim)
legend("error","simulation function")
xlabel("time (s)")

% Show an animation of the LIP
addpath("../balancer")
figure;
lip_traj = [lip_trajectory(1,:);1/m*lip_trajectory(4,:)];
animate_lip(1:dt:T, lip_traj', lip_control',h)
