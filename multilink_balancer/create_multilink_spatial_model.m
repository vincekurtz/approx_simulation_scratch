%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to create a model of a multi-link balancer using
% the spatial_v2 library
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Construct the spatial_v2 model
model = autoTree(3);

% Add an "arm"
%model.NB = 4;
%model.jtype = {'R', 'R', 'R', 'R'};
%model.parent = [0 1 2 2];
%model.Xtree{end+1} = model.Xtree{end};  % just copy existing spatial transforms and inertias
%model.I{end+1} = model.I{end};
%model.appearance.body{end+1} = model.appearance.body{end};

% Add gravity in the -y direction
model.gravity = [0;-9.81;0];

% Camera settings for showmotion
model.camera.direction = [0 0 1];
model.camera.up = [0 1 0];
model.camera.zoom = 0.5;
model.camera.locus = [0 0];

% Derive symbolic equations of motion
disp("===> Deriving Equations of Motion")
tau_sym = sym('tau', [model.NB,1], 'real');    % we'll assume all joints are actuated
q_sym = sym('theta', [model.NB,1], 'real');
qd_sym = sym('theta_dot', [model.NB,1], 'real');

u_sym = tau_sym;
x_sym = [q_sym;qd_sym];

qdd_sym = FDab(model, q_sym, qd_sym, tau_sym);
f_sym = [qd_sym; qdd_sym];
f_func = matlabFunction(f_sym, 'file', 'BalancerDynamics', 'vars', {x_sym, u_sym});

% Feedback linearization
disp("===> Computing Feedback Linearization")
ret = EnerMo(model, q_sym, qd_sym); % compute some relevant quantities

% Positions of center of mass and CoM jacobian
xc_sym = simplify(ret.cm(1));  
yc_sym = simplify(ret.cm(2));
pc_sym = [xc_sym;yc_sym];
J_sym = jacobian(pc_sym, q_sym);

vxc_sym = simplify(ret.vcm(1));
vyc_sym = simplify(ret.vcm(2));
vpc_sym = [vxc_sym;vyc_sym]; 
Jdot_sym = jacobian(vpc_sym, q_sym);

% System Equations of motion in the form H(q)qdd+C(q,qd,f_ext)=tau,
% where H is the joint-space inertia matrix and C is the vector of gravity,
% external force, and correolis terms
[H_sym, C_sym] = HandC(model, q_sym, qd_sym);

% Task-space inertia matrix
Lambda_sym = inv(J_sym*inv(H_sym)*J_sym');

% Generate matlab functions for relevant quantities
p_com_func = matlabFunction(pc_sym, 'file', 'p_com', 'vars', {q_sym});
pd_com_func = matlabFunction(vpc_sym, 'file', 'pd_com', 'vars', {q_sym, qd_sym});

J_func = matlabFunction(J_sym, 'file', 'J', 'vars', {q_sym});
Jdot_func = matlabFunction(Jdot_sym, 'file', 'Jdot', 'vars', {q_sym, qd_sym});

H_func = matlabFunction(H_sym, 'file', 'H', 'vars', {q_sym});
C_func = matlabFunction(C_sym, 'file', 'C', 'vars', {q_sym, qd_sym});
Lambda_func = matlabFunction(Lambda_sym, 'file', 'Lambda', 'vars', {q_sym});

% Simulate the systems
disp("===> Simulating System")
T = 10;     % simulation time
dt = 1e-2;  % timestep

% initial state
x = 0.5*ones(2*model.NB,1);   % start away from signularities
state_trajectory = [];
com_trajectory = [];

for t = 1:dt:T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % feedback controller
    tic
    q = x(1:model.NB);    % unpack joint state
    qd = x(model.NB+1:end);  

    x_lin = [p_com(q);J(q)*qd];
    x_des = [0;1.0;0;0];

    K = [1 0 1.7 0; 0 1 0 1.7];
    v = -K*(x_lin-x_des);       % virtual control for linearized system

    %u = H(q)*inv(J(q))*(v-Jdot(q,qd)*qd)+C(q,qd);  % feedback linearization
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

