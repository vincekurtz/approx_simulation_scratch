%% Simulation of balancing using the Linear Inverted Pendulum (LIP) model

% Parameters for the simulation
params.h = 0.7;
params.g = 9.81;
params.omega = sqrt(params.g/params.h);

params.A = [ 0 1; params.omega^2 0];  % x' = Ax+Bu
params.B = [1; -params.omega^2];

% Parameters for the controller: we'll use an LQR scheme
[K,P,e] = lqr(params.A,params.B,eye(2),1);
params.K = K;

% Initial position of th CoM
x0 = [0.1;0.1];

[t_sim, x_sim, u_sim] = Sim(x0, params);

animate(t_sim, x_sim, u_sim, params)

function [] = animate(t_sim, x_sim, u_sim, params)
    % Animate the behavior of the LIP model

    % Y height is fixed
    y = params.h;

    % Ground frame
    ground = plot(gca, [-15 15],[0 0],'k','LineWidth',2);
    frame0 = hgtransform(gca);

    % Center of pressure
    frame1 = hgtransform(frame0);
    r = 0.1;
    cop = rectangle('Curvature',[1,1],'Parent',frame1);  % Center of pressure
    cop.Position = [-r/2, -r/2, r, r];
    cop.EdgeColor = 'none';
    cop.FaceColor = 'blue';

    % Center of mass
    frame2 = hgtransform(frame1);
    cop = rectangle('Curvature',[1,1],'Parent',frame2);  % Center of pressure
    cop.Position = [-r/2, -r/2, r, r];
    cop.EdgeColor = 'none';
    cop.FaceColor = 'red';

    % Line between COP and COM
    l = line('Parent',frame1);


    axis equal
    xlim([-2.5,2.5])

    for i=1:length(t_sim)-1
        tic
        dt = t_sim(i+1) - t_sim(i);

        % Move the center of pressure
        u = u_sim(i);
        frame1.Matrix = makehgtform('translate',[u 0 0]);

        % move the center of mass
        x = x_sim(i,1);
        frame2.Matrix = makehgtform('translate',[x-u y 0]);

        l.XData = [0,x-u];
        l.YData = [0,y];

        pause(dt-toc)
    end


end

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


