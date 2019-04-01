% http://ftp.feq.ufu.br/Luis_Claudio/HYBRID/mpt_ver262/mpt/yalmip/htmldata/sos.htm

% Requires PENBMI

clear all
yalmip('clear');

% States...
sdpvar x1 x2
x = [x1;x2];

% Non-quadratic Lyapunov z'Pz 
z = [x1;x2;x1^2];
P = sdpvar(3,3);
V = z'*P*z;

% Non-linear state feedback
v = [x1;x2;x1^2];
K = sdpvar(1,3);
u = K*v;

% System x' = f(x)+Bu
f = [1.5*x1^2-0.5*x1^3-x2; 3*x1-x2];
B = [0;1];

% Closed loop system, u = Kv
fc = f+B*K*v;

% Stability and performance constraint dVdt < -x'x-u'u
% NOTE : This polynomial is bilinear in P and K
F = [sos(-jacobian(V,x)*fc-x'*x-u'*u)];

% P is positive definite, bound P and K for numerical reasons
F = F + [sos(V)]
F = F +[25>=P(:)>=-25]+[25>=K>=-25]

% Minimize trace(P)
% Parametric variables P and K automatically detected
% by YALMIP since they are both constrained
solvesos(F,trace(P))
value(K)