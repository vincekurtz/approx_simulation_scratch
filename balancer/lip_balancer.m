%% Simulation of balancing using the Linear Inverted Pendulum (LIP) model

% Parameters for the LIP model
lip_params.h = 0.7;
lip_params.g = 9.81;
lip_params.omega = sqrt(lip_params.g/lip_params.h);

lip_params.A = [ 0 1; lip_params.omega^2 0];  % x' = Ax+Bu
lip_params.B = [1; -lip_params.omega^2];

[K,P,e] = lqr(lip_params.A,lip_params.B,eye(2),1);  % We'll use an LQR control scheme to regulate about the upright
lip_params.K = K;

% Initial position of th CoM
x0 = [0.1;0.1];

[t_sim, x_sim, u_sim] = Sim(x0, lip_params);

animate_lip(t_sim, x_sim, u_sim, lip_params.h)

function [t_store, state_store, u_store] = Sim(x0, params)
    % Simulate the LIP model with the given initial condition, and 
    % control law as defined below.

    x = x0;
    dt = 1e-2;

    state_store = [];
    t_store = [];
    u_store = [];

    for t = 0:dt:5
        % A brief disturbance
        if t == 3
            x = [1;0.1];
        end
        % Compute the control input
        u = control_law(x, params);

        % Save the state and control
        t_store(end+1) = t;
        state_store(end+1,:) = x';
        u_store(end+1,:) = u;

        % Euler integration
        dx = params.A*x + params.B*u;
        x = x+dx*dt;
    end
end

function u = control_law(x, params)
    % We'll choose to place the CoP at the "caputre point"
    u = -params.K*x;
end

