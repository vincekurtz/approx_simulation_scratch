%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate an interface that will guarantee that the full robot dynamics approximately 
% simulate the LIP dynamics.
%
% This script will define two functions:
%   interface_fcn : u_lip x_lip x_robot --> u_robot
%
%   simulation_fcn : x_lip, x_robt --> V
%
%   where 
%     u_lip is the x-position of the CoP of the LIP model
%     x_lip is the position of the CoM of the LIP model
%     x_robot is the double pendulum's joint angles and velocities
%     u_robot is the double pendulum's joint torques
%     V is the simulation function 
%             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Definitions                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear Inverted Pendulum
h = 0.75;
g = 9.81;
omega = sqrt(g/h);

A_lip = [0 1; omega^2 0];
B_lip = [0; -omega^2];

% Stabilizing controller to balance the LIP
Q_lip = diag([10;1]);
R_lip = 1.0;
K_lip = lqr(A_lip, B_lip, Q_lip, R_lip);  % Stabilizing gain to balance the LIP

% Double Pendulum
syms theta1 theta2 theta1_dot theta2_dot real;
syms tau1 tau2 real
q_sym = [theta1;theta2];
qd_sym = [theta1_dot;theta2_dot];

x_sym = [q_sym;qd_sym];
u_sym = [tau1;tau2];

dp_model = autoTree(2);   % an unbranched kinematic tree with 2 bodies
m = 2;  % total mass
dp_model.gravity = [0;-g;0];
qdd_sym = FDab( dp_model, q_sym, qd_sym, u_sym);
f_sym = [qd_sym; qdd_sym];
f_func = matlabFunction(f_sym,'vars',{x_sym, u_sym});

% Feedback Linearized Inverted Pendulum
xc_sym = 3/4*cos(theta1)+1/4*cos(theta1+theta2);   % COM position
yc_sym = 3/4*sin(theta1)+1/4*sin(theta1+theta2);
pc_sym = [xc_sym;yc_sym];
pc_func = matlabFunction(pc_sym, 'vars',{x_sym});

J_sym = jacobian(pc_sym, q_sym);   % Jacobian maps state to COM velocities
J = matlabFunction(J_sym, 'vars', {q_sym});

[H_sym, C_sym] = HandC(dp_model, q_sym, qd_sym);  % Parameters of the full system dynamics for feedback linearization
H = matlabFunction(H_sym, 'vars', {q_sym});
C = matlabFunction(C_sym, 'vars', {q_sym, qd_sym});

% TODO: figure out how to compute this symbolically
Jdot_sym =  [- (3*cos(theta1)*theta1_dot)/4 - (cos(theta1 + theta2)*(theta1_dot + theta2_dot))/4, -(cos(theta1 + theta2)*(theta1_dot + theta2_dot))/4;
            - (3*sin(theta1)*theta1_dot)/4 - (sin(theta1 + theta2)*(theta1_dot + theta2_dot))/4, -(sin(theta1 + theta2)*(theta1_dot + theta2_dot))/4];
Jdot = matlabFunction(Jdot_sym, 'vars', {q_sym, qd_sym});

% Convinience function that maps joint states to center of mass position
% and velocity, ie [x_c,x_c',y_c,y_c']
pc_dot_sym = J_sym*qd_sym;   % center of mass velocity
x_com_sym = [pc_sym(1);pc_dot_sym(1);pc_sym(2);pc_dot_sym(2)];
x_com = matlabFunction(x_com_sym, 'vars', {x_sym});

% Define transformations from ground contacts to the center of mass
foot_width = 0.75;          % conservative under-approximation of foot width
c1 = [foot_width/2 0 0];    % we'll consider two contacts along the x-axis
c2 = [-foot_width/2 0 0];

% Spatial transforms from contact to base frame
o_X_c1 = xlt(c1);  
o_X_c2 = xlt(c2);

% Spatial transforms from base frame to COM
com_X_o = xlt([xc_sym yc_sym 0]);

% Spatial transforms from contact to COM
com_X_c1_sym = com_X_o*o_X_c1;
com_X_c2_sym = com_X_o*o_X_c2;
com_X_c1 = matlabFunction(com_X_c1_sym, 'vars', {x_sym});  % express as a function of [q,qdot]
com_X_c2 = matlabFunction(com_X_c2_sym, 'vars', {x_sym});

% Force transforms from contact to COM
com_Xf_c1_sym = inv(com_X_c1_sym)';
com_Xf_c2_sym = inv(com_X_c2_sym)';
com_Xf_c1 = matlabFunction(com_Xf_c1_sym, 'vars', {x_sym});  % express as a function of [q,qdot]
com_Xf_c2 = matlabFunction(com_Xf_c2_sym, 'vars', {x_sym});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface Design                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Concrete system is a double integrator in x and y position, which
% is obtained from task-space linearization of the double pendulum
A1 = [0 1 0 0;
      0 0 0 0;
      0 0 0 1;
      0 0 0 0];
B1 = [0 0;
      1 0;
      0 0;
      0 1];
C1 = eye(4);     % the outputs of both systems are the full states

% The abstract system is the linear inverted pendulum
A2 = [0       1 0 0;
      omega^2 0 0 0;
      0       0 0 1;
      0       0 0 0];
B2 = [0        ;
      -omega^2 ;
      0        ;
      0        ];
C2 = eye(4);

% external control is position of CoP for the lip
u = sdpvar(1);   

% State of the concrete system
x1 = sdpvar(4,1);

% State of the abstract system
x2 = sdpvar(4,1);

% Interface definition
P = eye(4);  % P and Q chosen such that PA2 = A1P+B1*Q and C2 = C1*P
Q = [omega^2 0 0 0;
     0       0 0 0];
R = [-omega^2;0];
K_joint = -lqr(A1,B1,eye(4),eye(2));    % A control gain that stabilizes the concrete system
                                            % Could derive from a PD controller as well

lambda = 0.1;
Mbar = sdpvar(4,4);
Kbar = K_joint*Mbar;

F = [ [Mbar Mbar*C1'; C1*Mbar eye(4)] >= 0];
F = F + [Mbar*A1' + A1*Mbar + Kbar'*B1'+B1*Kbar + 2*lambda*Mbar <= 0 ];
optimize(F);

M = inv(value(Mbar));

% Double check the original conditions
check1 = all(eig(M-C1'*C1) >= 0);   % M >= C'C
check2 = all(eig((A1+B1*K_joint)'*M + M*(A1+B1*K_joint) + 2*lambda*M) <= 0);  % (A+BK)'M+M(A+BK) <= -2lambdaM

if check1 & check2
    disp("Checks passed!")
else
    disp("Checks Failed!")
    return;
end


% Simulation function
disp("V = (x1-x2)'*M*(x1-x2)")
disp("M: ")
M
disp("")

% Class K function of control input
gmma = norm(sqrt(M)*(B1*R-P*B2))/lambda;
disp("Gamma: ")
sdisplay(gmma)
disp("")

% Interface in terms of joint states and lip state
u_lip_sym = sym('u_lip');
x_lip_sym = sym('x_lip',[4,1]);
% x_com_sym defined above in terms of joint angles and velocities

% Linear interface u1 = Ru2 + Qx2 + K(x1-Px2)
u_com = R*u_lip_sym + Q*x_lip_sym + K_joint*(x_com_sym-P*x_lip_sym); 

% Feedback linearizing map to joint torques
tau_sym = H_sym*inv(J_sym)*(u_com-Jdot_sym*qd_sym)+C_sym;

interface_fcn = matlabFunction(tau_sym, 'file', 'InterfaceFcn', 'vars',{u_lip_sym, x_lip_sym, x_sym});

% Simulation function in terms of joint angles and velocities
V = (P*x_lip_sym-x_com_sym)'*M*(P*x_lip_sym-x_com_sym);
simulation_fcn = matlabFunction(V, 'file', 'SimulationFcn', 'vars', {x_lip_sym, x_sym});

