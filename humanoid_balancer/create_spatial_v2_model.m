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
%   - h_0.m              : centroid momentum
%   - A_com.m            : Centroid dynamics matrix
%   - Ad_com_qd.m        : Centroid dynamics matrix (time derivative)
%   - Xf_1.m             : Force transform from contact 1 to CoM
%   - Xf_2.m             : Force transform from contact 2 to CoM
%   - Xf_0.m             : Force transform from ground frame to CoM
%   - A_cwc.m            : Contact constraint on hd_com as a function of p_com
%   - b_cwc.m            : Contact constraint on hd_com as a function of p_com
%   - balancer_model.mat : spatial_v2 model of the balancer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Construct the spatial_v2 model
model = autoTree(4);

% Extend the "torso" link
model.appearance.body{3}{2} = [0 0 0; 2 0 0];
model.I{3} = mcI(2,[1;0;0], [0.0025       0       0;
                             0       0.0846       0;
                             0            0  0.0846]);
                           
% Add gravity in the -y direction
model.gravity = [0;-9.81;0];

% Camera settings for showmotion
model.camera.direction = [0 0 1];
model.camera.up = [0 1 0];
model.camera.zoom = 0.5;
model.camera.locus = [0 0];

% Determine total mass of the balancer
m = 0;
for i=1:model.NB
    m = m + mcI(model.I{i});
end

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
com_Xf_o = inv(com_X_o)';         % force transform from base frame to CoM

com_X_c1 = com_X_o*o_X_c1;  % transform from contacts to CoM
com_X_c2 = com_X_o*o_X_c2;

com_Xf_c1 = inv(com_X_c1)'; % force transforms
com_Xf_c2 = inv(com_X_c2)';

% Spatial Momentum
h_o = ret.htot;   % spatial momentum in the base frame (for double checking)

fb_model = floatbase(model);   % Centroid dynamics following Wensing and Orin 2015
[HH_sym, CC_sym] = HandC(fb_model, [zeros(5,1);q_sym], [zeros(5,1);qd_sym]);
Psi = [zeros(3) eye(3); eye(3) zeros(3)];
U = [eye(6) zeros(6,3)];
A_com_sym = com_Xf_o*Psi'*U*HH_sym;
A_com_sym = A_com_sym*[zeros(5,model.NB);eye(model.NB)];  % Assume floating base joints are fixed at zero
A_com_sym = [zeros(3,2) eye(3) zeros(3,1)]*A_com_sym;     % select only nonzero momentum terms
mg = [0;0;0;m*model.gravity];
Ad_com_qd_sym = com_Xf_o*Psi'*U*CC_sym+mg;
Ad_com_qd_sym = [zeros(3,2) eye(3) zeros(3,1)]*Ad_com_qd_sym;     % select only nonzero momentum terms

% Contact Constraints 
mu = 0.2;     % friction coefficient
l = 0.5;      % foot width

p_com = sym('p_com', [2,1], 'real');  % center of mass
o_Xf_com_p = [eye(3)  [0         0        p_com(2) ;   % Force transform from the CoM frame to the 0 frame
                       0         0        -p_com(1);   % in terms of position of CoM
                       -p_com(2) p_com(1) 0       ];
              zeros(3)   eye(3)                     ];

% Constraints on the contact wrench (note that we only consider motion in the x-y plane)
A = [0 0  0  0  -1 0;   % positive normal force
     0 0  0  1 -mu 0;   % Coulomb friction
     0 0  0 -1 -mu 0;
     0 0  1  0  -l 0;   % Center of pressure constraint
     0 0 -1  0  -l 0];

% Constraints on u_com
A_cwc_sym = A*o_Xf_com_p;
A_cwc_sym = A_cwc_sym*[zeros(2,3);eye(3);zeros(1,3)]; % map u_com to [0;0;u_com;0] implicitly
b_cwc_sym = A*o_Xf_com_p*mg;


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

h_func = matlabFunction(h_o, 'file', 'h_o', 'vars', {q_sym, qd_sym});
A_com_func = matlabFunction(A_com_sym, 'file', 'A_com', 'vars', {q_sym});
Ad_com_qd_func = matlabFunction(Ad_com_qd_sym, 'file', 'Ad_com_qd', 'vars', {q_sym, qd_sym});

Xf_1_func = matlabFunction(com_Xf_c1, 'file', 'Xf_1', 'vars', {q_sym});
Xf_2_func = matlabFunction(com_Xf_c2, 'file', 'Xf_2', 'vars', {q_sym});

com_Xf_0_func = matlabFunction(com_Xf_o, 'file', 'Xf_0', 'vars', {q_sym});

A_cwc_func = matlabFunction(A_cwc_sym, 'file', 'A_cwc', 'vars', {p_com});
b_cwc_func = matlabFunction(b_cwc_sym, 'file', 'b_cwc', 'vars', {p_com});

% Save the model
save('balancer_model','model')
