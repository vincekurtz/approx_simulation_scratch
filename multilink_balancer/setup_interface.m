%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define parameters of the LIP model and design an associated 
% interface.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIP model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1.5;  % height of lip model
g = 9.81; % acceleration due to gravity
omega = sqrt(g/h);
A_lip = [0 1; omega^2 0];
B_lip = [0; -omega^2];

% Stabilizing controller for LIP model
Q_lip = diag([1,1]);
R_lip = 1.0;
K_lip = lqr(A_lip, B_lip, Q_lip, R_lip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface Design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Concrete model is a linearization of the centroid dynamics, where
% the input is the hd_com.
m = 1*model.NB;  % total mass of the balancer
A1 = [zeros(2,3) 1/m*eye(2); 
      zeros(3,5)           ];
B1 = [zeros(2,3);
      eye(3)    ];
C1 = eye(5);  % outputs of both systems are the full states

% Abstract model is the LIP, including (dummy) y position and velocity
A2 = [zeros(2,3) 1/m*eye(2);
      zeros(1,5)           ;
      omega^2 zeros(1,4)   ;
      zeros(1,5)           ];
B2 = [zeros(3,1);
      -omega^2  ;
      0         ];
C2 = eye(5);

% Symbolic variables for interface design
u = sdpvar(1);   % external control is position of CoP of the LIP
x1 = sdpvar(5,1);
x2 = sdpvar(5,1);

% Interface Derivation
P = eye(5);  % P and Q chosen such that PA2 = A1P+B1*Q and C2 = C1*P
Q = [0       0 0 0 0;
     omega^2 0 0 0 0;
     0       0 0 0 0];
R = [0;-omega^2;0];
K_joint = -lqr(A1,B1,diag([1 1 1 1 1]),0.01*eye(3));    % A control gain that stabilizes the concrete system
                                        % Could derive from a PD controller as well

lambda = 0.1;
Mbar = sdpvar(5,5);
Kbar = K_joint*Mbar;

F = [ [Mbar Mbar*C1'; C1*Mbar eye(5)] >= 0];
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

% Now we have an interface given by 
% u_com = R*u + Q*x2 + K_joint*(x1-P*x2);

% Simulation function
V = (P*x2-x1)'*M*(P*x2-x1);
disp("V = (x1-x2)'*M*(x1-x2)")
disp("M: ")
M
disp("")

% Class K function of control input
gmma = norm(sqrt(M)*(B1*R-P*B2))/lambda;
disp("Gamma: ")
sdisplay(gmma)
disp("")
