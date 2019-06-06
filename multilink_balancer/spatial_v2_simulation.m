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
T = 10;     % simulation time
dt = 1e-2;  % timestep

% Initial state of full model
q0 = [0.5;0.2;0.2;0.4];   
qd0 = [0.1;0;0;0];
x = [q0;qd0];

% Initial state of LIP model
com_pos = p_com(q0);
com_vel = pd_com(q0,qd0);
x_lip = [com_pos(1);h;0;m*com_vel(1);0];

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

    % Compute LIP control that will let us balance
    u_lip = -K_lip*[x_lip(1);1/m*x_lip(4)];

    % Linearization based on centroid momentum
    A = [zeros(3,2) eye(3) zeros(3,1)]*A_com(q);    % select only nonzero terms
    Ad_qd = [zeros(3,2) eye(3) zeros(3,1)]*Ad_com_qd(q,qd);

    % Compute the state of the linearized system
    com_pos = p_com(q);
    com_h = A*qd;
    x_lin = [com_pos;com_h];

    % Compute a control to track the LIP trajectory
    u_com = R*u_lip + Q*x_lip + K_joint*(x_lin-P*x_lip);

    % Feedback linearization to generate torques as actual control
    Lambda = inv(A*inv(H(q))*A');
    u = A'*(Lambda*u_com - Lambda*Ad_qd + Lambda*A*inv(H(q))*C(q,qd));

    % Secondary control via null space projector
    Abar = inv(H(q))*A'*Lambda;
    N = (eye(4) - A'*Abar')';   % null space projector

    kp = diag([1;0;0;0]);    % we'll use PD control of joints to
    kd = diag([1;1;1;1]);    % apply commands in the null space
    q_des = [pi/2-0.1;0;0;0];
    qd_des = [0;0;0;0];
    u0 = -kp*(q-q_des)-kd*(qd-qd_des);
    u = u + N'*u0;

    % Simulate the full model forward
    xdot = BalancerDynamics(x,u);
    x = x + xdot*dt;

    % Simulate the lip model forward
    xdot_lip = A2*x_lip + B2*u_lip;
    x_lip = x_lip + xdot_lip*dt;

    state_trajectory(:,end+1) = x;
    com_trajectory(:,end+1) = x_lin;
    lip_trajectory(:,end+1) = x_lip;
    lip_control(:,end+1) = u_lip;
    err_sim(:,end+1) = sqrt((x_lin-x_lip)'*(x_lin-x_lip));
    V_sim(:,end+1) = sqrt((x_lin-P*x_lip)'*M*(x_lin-P*x_lip));
    
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
