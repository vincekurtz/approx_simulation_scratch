%% Finding a quadratic simulation function for linear systems
% Following Girard & Pappas, Hierarchical control system design using
% approximate simulation. Automatica 2009.

close all;
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
x1 = sdpvar(1)
x2 = sdpvar(1)
x3 = sdpvar(1);

V = (P*x3-[x1;x2])'*M*(P*x3-[x1;x2]);     % really should be sqrt(V). Equivalent?
disp("Simulation Function: ")
sdisplay(V,[3])
disp(' ')

% And the interface
R = 1;   % user defined
u2 = sdpvar(1);
u2 = 0.1; % DEBUG
interface = u2 + K*([x1;x2]-P*x3);
disp("Interface: ")
sdisplay(interface,[3])
disp(" ")

% We can also find the class K function gamma
gamma = norm(sqrt(M)*(B1*R-P*B2))/lambda*abs(u2);
disp("Gamma function: ")
sdisplay(gamma,[3])

%% Check if this known solution satisfies the SOS constraints

% Define input and output equations of the joint system
f = [A1*[x1;x2] + B1*interface;
     A2*x3 + B2*u2];
 
% System outputs
g = x1-x3;

%% To double check our results, we can simulate our system  

% Simulation parameters
dt = 0.1;
Ns = 500;

% Place to record the results
x0 = [0.2;0;0];   % initial conditions

% Place to record the results
x1_sim = zeros(2,Ns);
x2_sim = zeros(1,Ns);
u2 = zeros(1,Ns);

% Initial conditions
x1_sim(:,1) = x0(1:2);
x2_sim(:,1) = x0(3);
for i = 1:(Ns-1)
    u2(:,i) = 0.1*sin(0.03*i);  % control pre-specified
    
    if i>300
        u2(:,i) = 0;
    end
    
    % Simulate x1 forward in time
    interface = u2(:,i) + K*(x1_sim(:,i)-P*x2_sim(:,i));
    x1_dot = A1*x1_sim(:,i)+B1*interface;
    x1_sim(:,i+1) = x1_sim(:,i) + x1_dot*dt;
    
    % Simulate x2 forward in time
    x2_dot = A2*x2_sim(:,i)+B2*u2(:,i);
    x2_sim(:,i+1) = x2_sim(:,i) + x2_dot*dt;
    
end

% Plot x1 vs x2
figure;
hold on;
plot(x1_sim(1,:))
plot(x2_sim(1,:))
xlabel('time')
ylabel('position')
legend("x_1: concrete","x_2: abstract")
title("System Trajectories")

% Plot the simulation function
V_sim = zeros(1,Ns);
for i = 1:Ns
    V_sim(:,i)= replace(V,[x1;x2;x3],[x1_sim(:,i);x2_sim(:,i)]);
end

figure;
hold on;
plot(V_sim)
legend("V(x)")
title("Simulation Function")
xlabel('time')
