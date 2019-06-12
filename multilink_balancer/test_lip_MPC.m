%% Quick script to do some debugging with mpc by just simulating the linearized systems

clear;close;clc;

load('balancer_model')
setup_interface

% Simulation parameters
T = 10;     % simulation time
dt = 1e-1;  % timestep

% Initial state of full model
x_com = [0.2011    1.6808    0.1594    0.1326   -0.4532]';
x_lip = [0.2044    1.5000         0    0.1660         0]';

% Quantities to save for plotting
com_trajectory = [];
lip_trajectory = [];

for t = 1:dt:T
    tic
    % Compute LIP control that will let us balance
    params.A_lip = A2;
    params.B_lip = B2;
    params.A_com = A1;
    params.B_com = B1;
    params.N = 2;
    params.dt = dt;
    params.R = R;
    params.Q = Q;
    params.K = K_joint;
    params.m = m;

    u_lip_trajectory = GenerateLIPTrajectory(x_lip, x_com, params);

    u_lip = u_lip_trajectory(1);

    % Compute a control to track the LIP trajectory
    u_com = R*u_lip + Q*x_lip + K_joint*(x_com-P*x_lip);

    % Apply both controls
    x_com = x_com + (A1*x_com + B1*u_com)*dt;
    x_lip = x_lip + (A2*x_lip + B2*u_lip)*dt;

    % Record trajectories
    com_trajectory(:,end+1) = x_com;
    lip_trajectory(:,end+1) = x_lip;
    toc
end

plot(1:dt:T,vecnorm(com_trajectory-lip_trajectory))
xlabel("time")
ylabel("error")
