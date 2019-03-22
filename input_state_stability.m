% We can formulate the simulation problem as Input-Output stability of the
% joint system, so here we'll consider starting from an ISS perspective

close all;
clear all;
clc;

%% Original two systems and interface
% Systems in question are a double integrator
A1 = [0 1;
      0 0];
B1 = [0;1];
C1 = [1 0];

% And a single integrator
A2 = 0;
B2 = 1;
C2 = 1;

% An interface that maps the input of system 2 to an input of system 1
% is u2 + K*(x1-P*x2), where
K = [-1, -1];
P = [1;0];

%% ISS of joint system system

u = sdpvar(1);

x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

x = [x1;x2;x3];
x0 = [0.2;0;0];   % initial conditions

% System definition
f = [A1*[x1;x2] + B1*[u+K*([x1;x2]-P*x3)];
     A2*x3 + B2*u];
 
% System outputs
g = x1-x3;

% Candidate ISS Lyapunov function
[V,vc] = polynomial(x,2);
Vdot = jacobian(V,x)*f;

% Class K_inf functions should be even polynomials
a1c = sdpvar(2,1);
a1poly = monolist(g,2).^2;
a1poly = a1poly(2:end);
alpha1 = a1c'*a1poly;

[alpha1,a1c,da1] = even_polynomial(g,2);
[alpha2,a2c,da2] = even_polynomial(x,2);
[alpha3,a3c,da3] = even_polynomial(x,2);
[alpha4,a4c,da4] = even_polynomial(u,2);

% To help with numerical issues
epsilon = 1e-8;

% SOS constraints
F = [sos(V-alpha1-epsilon*x'*x)];
F = F + [sos(alpha2-V-epsilon*x'*x)];
F = F + [sos(-alpha3+alpha4-Vdot)];
F = F + [sos(da1),sos(da2),sos(da3),sos(da4)];

[sol] = solvesos(F,[],[],[vc;a1c;a2c;a3c;a4c]);

if ~sol.problem
    disp("Solution Found!")
    disp("V: ")
    sdisplay(replace(V,vc,value(vc)),[3])
else
    disp("No solution found :(")
end

%% Simulating both systems to check the results

% Simulation parameters
dt = 0.1;
Ns = 500;

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
legend("x_1: concrete","x_2: abstract")
title("System Trajectories")

% Plot some of the quantities we've calculated to certify simulation
figure;
hold on;

V_sim = zeros(1,Ns);
Vdot_sim = zeros(1,Ns);
alpha1_sim = zeros(1,Ns);
alpha2_sim = zeros(1,Ns);
alpha3_sim = zeros(1,Ns);
alpha4_sim = zeros(1,Ns);
for i = 1:Ns
    V_sim(:,i)= replace(V,[vc;x1;x2;x3],[value(vc);x1_sim(:,i);x2_sim(:,i)]);
    Vdot_sim(:,i)= replace(Vdot,[vc;x1;x2;x3;u],[value(vc);x1_sim(:,i);x2_sim(:,i);u2(:,i)]);
    alpha1_sim(:,i) = replace(alpha1,[a1c;x1;x2;x3],[value(a1c);x1_sim(:,i);x2_sim(:,i)]);
    alpha2_sim(:,i) = replace(alpha2,[a2c;x1;x2;x3],[value(a2c);x1_sim(:,i);x2_sim(:,i)]);
    alpha3_sim(:,i) = replace(alpha3,[a3c;x1;x2;x3],[value(a3c);x1_sim(:,i);x2_sim(:,i)]);
    alpha4_sim(:,i) = replace(alpha4,[a4c;x1;x2;x3;u],[value(a4c);x1_sim(:,i);x2_sim(:,i);u2(:,i)]);
end

title("\alpha_1 \leq V \leq \alpha_2")
plot(alpha2_sim)
plot(V_sim)
plot(alpha1_sim)
legend("\alpha_2","V","\alpha_1")

figure;
hold on

title("Vdot \leq -\alpha_3 + \alpha_4")
plot(Vdot_sim)
plot(-alpha3_sim+alpha4_sim)
legend("Vdot","-\alpha_3+\alpha_4")
error = x1_sim(1,:)-x2_sim(1,:);

figure;
hold on
title("Class K_{\infty} functions")
s = linspace(0,10);
a1 = value(a1c)'*[s.^2;s.^4];
a2 = value(a2c)'*[s.^2;s.^4];
a3 = value(a3c)'*[s.^2;s.^4];
a4 = value(a4c)'*[s.^2;s.^4];
plot(s,a1)
plot(s,a2)
plot(s,a3)
plot(s,a4)

legend("\alpha_1","\alpha_2","\alpha_3","\alpha_4")
xlabel("s")
ylabel("\alpha(s)")


function [p, c, sdp_ds] = even_polynomial(s, n)
    % Helper function to construct a symbolic even polynomial
    %
    % Arguments:
    %    s - a YALMIP symbolic variable
    %    n - the number of coefficients
    %
    % Returns:
    %    p      - a symbolic even polynomial of the form p = c1*s^2+c2s^4+...
    %    c      - the associated constants
    %    sdp_ds - the polynomial s*dp(s)/ds, which forms a sos constraint
    %             for p to be a class K_inf function

    c = sdpvar(n,1);
    raw_poly = [s'*s];
    for i=2:n
        raw_poly(i,:) = raw_poly(1,:)^(i);
    end
    p = c'*raw_poly;

    cds = [2:2:2*n]'.*c    % coefficients of s*dp(s)/ds
    sdp_ds = cds'*raw_poly;

    %sdp_ds = jacobian(p,s);
end

