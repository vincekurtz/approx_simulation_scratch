%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to create a model of a multi-link balancer using
% the spatial_v2 library.
%
% Generates and saves the following files:
%   - BalancerDynamics.m : equations of motion of the balancer
%   - p_com.m            : position of the center of mass as
%                           a function of joint angles q
%   - pd_com.m           : velocity of the center of mass as
%                           a function of q and qdot
%   - J.m                : Center of mass jacobian
%   - Jdot.m             : time derivative of CoM jacobian
%   - H.m                : joint-space intertia matrix
%   - C.m                : correolis and gravitational terms
%   - Lambda.m           : task-space (CoM) inertial matrix
%   - balancer_model.mat : spatial_v2 model of the balancer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Construct the spatial_v2 model
model = autoTree(3);

% Add an "arm"
model.NB = 4;
model.jtype = {'R', 'R', 'R', 'R'};
model.parent = [0 1 2 2];
model.Xtree{end+1} = model.Xtree{end};  % just copy existing spatial transforms and inertias
model.I{end+1} = model.I{end};
model.appearance.body{end+1} = model.appearance.body{end};

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

% Spatial force transforms from contacts to CoM
foot_width = 0.5;
c1 = [foot_width/2 0 0];
c2 = [-foot_width/2 0 0];

o_X_c1 = xlt(c1);   % spatial transforms from contact to base frame
o_X_c2 = xlt(c2);

com_X_o = xlt([xc_sym yc_sym 0]); % transform from base frame to CoM

com_X_c1 = com_X_o*o_X_c1;  % transform from contacts to CoM
com_X_c2 = com_X_o*o_X_c2;

com_Xf_c1 = inv(com_X_c1)'; % force transforms
com_Xf_c2 = inv(com_X_c2)';

% Generate matlab functions for relevant quantities
disp("===> Saving Functions and Model File")
f_func = matlabFunction(f_sym, 'file', 'BalancerDynamics', 'vars', {x_sym, u_sym});

p_com_func = matlabFunction(pc_sym, 'file', 'p_com', 'vars', {q_sym});
pd_com_func = matlabFunction(vpc_sym, 'file', 'pd_com', 'vars', {q_sym, qd_sym});

J_func = matlabFunction(J_sym, 'file', 'J', 'vars', {q_sym});
Jdot_func = matlabFunction(Jdot_sym, 'file', 'Jdot', 'vars', {q_sym, qd_sym});

H_func = matlabFunction(H_sym, 'file', 'H', 'vars', {q_sym});
C_func = matlabFunction(C_sym, 'file', 'C', 'vars', {q_sym, qd_sym});
Lambda_func = matlabFunction(Lambda_sym, 'file', 'Lambda', 'vars', {q_sym});

Xf_1_func = matlabFunction(com_Xf_c1, 'file', 'Xf_1', 'vars', {q_sym});
Xf_2_func = matlabFunction(com_Xf_c2, 'file', 'Xf_2', 'vars', {q_sym});

% Save the model
save('balancer_model','model')
