%% Finding bisimulation function for linear systems
% Following Girard & Pappas, Hierarchical control system design using
% approximate simulation. Automatica 2009.

clear all;
clc;

% Systems in question are a double integrator
A1 = [0 1;
      0 0];
B1 = [0;1];
C1 = [1 0];

% And a single integrator
A2 = 0;
B2 = 1;
C2 = 1;

%% Finding M, lambda
% First we'll find matrix M to satisfy the matrix inequalities of Lemma 1. 

% First we'll solve equation (7) to get Mbar = M^-1

K = [-1 -1];  % K and lambda are specified a priori
lambda = 0.1;
Mbar = sdpvar(2,2);
Kbar = K*Mbar;

F = [ [Mbar  Mbar*C1';
        C1*Mbar eye(1)] >= 0];
F = [F, Mbar*A1' + A1*Mbar + Kbar'*B1'+B1*Kbar + 2*lambda*Mbar <= 0 ]
optimize(F)

M = inv(value(Mbar));

% Now we'll double check to make sure the original conditions hold
eig(M-C1'*C1) >= 0   % M >= C'C
eig((A1+B1*K)'*M + M*(A1+B1*K) + 2*lambda*M) <= 0  % (A+BK)'M+M(A+BK) <= -2lambdaM

%% Finding the bisimulation function and interface

% First, we need to find P, Q such that (8) and (9) hold

% We'll do this by hand
P = [1;0];
Q = 0;

% And double check of course
P*A2 == A1*P+B1*Q;
C2 == C1*P;

% Now we can find the Bisimulation function
x1 = sdpvar(2,1);
x2 = sdpvar(1);

V = (P*x2-x1)'*M*(P*x2-x1);     % really should be sqrt(V). Equivalent?
disp("Simulation Function: ")
sdisplay(V,[3])
disp(' ')

% And the interface
R = 1;   % user defined
u2 = sdpvar(1);
u2 = 0.1; % DEBUG
interface = u2 + K*(x1-P*x2);
disp("Interface: ")
sdisplay(interface,[3])
disp(" ")

% We can also find the class K function gamma
gamma = norm(sqrt(M)*(B1*R-P*B2))/lambda*abs(u2);
disp("Gamma function: ")
sdisplay(gamma,[3])

%% Check if this known solution satisfies the SOS constraints

% Define input and output equations
f1 = A1*x1+B1*interface;  
f2 = A2*x2+B2*u2;
g1 = [1 0]*x1;
g2 = [1]*x2;

% Now we define the partial of V w.r.t x1 and x2
nablaV_x1 = jacobian(V,x1);
nablaV_x2 = jacobian(V,x2); 

% Simple Simulation
% % S-procedure
% [s,cs] = polynomial([u2;x1;x2],3);
% 
% % The constraints
% F = sos(V-(g1-g2)'*(g1-g2));
% F = F + [sos(s)];
% F = F + [sos( -nablaV_x1*f1-nablaV_x2*f2 + s*(gamma-V))];
% [sol,u,Q,res] = solvesos(F,[],[],[cs]);


% Bi-simulation
l = sdpvar(1);
g = sdpvar(1);
F = sos(V-(g1-g2)'*(g1-g2));
F = F + [sos( -nablaV_x1*f1-nablaV_x2*f2 - l*V + g*(interface-u2)^2)];
[sol,u,Q,res] = solvesos(F,[],[],[l;g]);

if ~sol.problem
    disp(" ")
    disp("SOS Constraints Satisfied!")
else
    disp(" ")
    disp("Constraints Violated :(")
end

%% To double check our results, we can simulate our system  

% Simulation parameters
dt = 0.1;
Ns = 500;

% Place to record the results
x1 = zeros(2,Ns);
x2 = zeros(1,Ns);
u2 = zeros(1,Ns);

% Initial conditions
x1(:,1) = [0.2;0];
x2(:,1) = 0;
for i = 1:(Ns-1)
    u2(:,i) = 0.1*sin(0.03*i);  % control pre-specified
    
    if i>300
        u2(:,i) = 0;
    end
    
    % Simulate x1 forward in time
    interface = u2(:,i) + K*(x1(:,i)-P*x2(:,i));
    x1_dot = A1*x1(:,i)+B1*interface;
    x1(:,i+1) = x1(:,i) + x1_dot*dt;
    
    % Simulate x2 forward in time
    x2_dot = A2*x2(:,i)+B2*u2(:,i);
    x2(:,i+1) = x2(:,i) + x2_dot*dt;
    
end

% Plot the resulting trajectories: x1 should follow x2
hold on
plot(x1(1,:))
plot(x2(1,:))
legend("x_1","x_2")

% Plot the resulting error
figure;
hold on
error = (x1(1,:)-x2(1,:)).^2;
plot(error)

% And the bisimulation function, which should bound the error
V = zeros(1,Ns);
for i = 1:(Ns-1)
    % Copied the bisimulation function by hand for now
    V(:,i) = 1.929*x1(1,i)^2+1.458*x1(1,i)*x1(2,i)-3.858*x1(1,i)*x2(i)+1.884*x1(2,i)^2-1.458*x1(2,i)*x2(i)+1.929*x2(i)^2;
end
plot(V)

legend("||g_1(x_1)-g_2(x_2)||","V(x_1,x_2)")

% Note that the bisimulation function is nondecreasing only when
% gamma(u)<V(x1,x2)

figure;
hold on;
gamma = norm(sqrt(M)*(B1*R-P*B2))/lambda * abs(u2);
plot(gamma)
plot(V)
legend("gamma(||u_2||)","V(x_1,x_2)")

figure;
hold on;
% Plot difference in outputs
plot(error)
% vs bound for bisimulation
u1 = u2(1,:)-x1(1,:)-x1(2,:)+x2(1,:);
plot(norm(sqrt(M)*(B1*R-P*B2))/lambda*(u2(1,:)-u1(1,:)).^2);
legend("error","gamma term")

