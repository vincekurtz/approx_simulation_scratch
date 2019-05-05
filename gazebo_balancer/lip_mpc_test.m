% Script to test our LIP MPC controller purely on the lip model

clear all;clc;

% System definition
h = 0.75;
g = 9.812;
omega = sqrt(g/h);

A = [0 1; omega^2 0];
B = [0 ; -omega^2];

% Simulation
dt = 3e-2;
T = 5;

state_trajectory = [];
control_trajectory = [];
x = [0.6517;0.3131];
for t=0:dt:T
    disp(t)

    u = LIPController(x, omega, dt);

    xdot = A*x + B*u;
    x = x + xdot*dt;

    state_trajectory(end+1,:) = x;
    control_trajectory(end+1,:) = u;
end

addpath("../balancer")
animate_lip(0:dt:T,state_trajectory, control_trajectory, h);
