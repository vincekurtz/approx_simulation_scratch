
% Define an interface that ensures approximate simulation
compute_interface;

% Simulation parameters
T = 10;   % length of simulation
dt = 1e-2;

% For storing data to plot later
t_sim = [];
lip_state = [];
lip_control = [];
dp_state = [];

% Torque limits
tau_min = [-100;-100];
tau_max = [100;100];

% Initial state of the full system
x = [pi/2-0.8;0.2;0;0];         % initial joint states and velocities

% Initial state of the LIP model
com_init = x_com(x); 
x_lip = [com_init(1:2);h;0];  % we start the lip with same x value and velocity as the full system

for t = 1:dt:T

    % Compute the LIP control that will let us balance.
    % Only the x position and velocity are considered
    u_lip = -K_lip*x_lip(1:2);

    % Translate the LIP control to the full system control
    tau = InterfaceFcn(u_lip, x_lip, x);

    % Apply the controls to the full system
    dx = f_func(x, tau);
    x = x + dx*dt;
    
    % Apply the LIP control to the LIP model
    dx_lip = A2*x_lip + B2*u_lip;
    x_lip = x_lip + dx_lip*dt;

    % Record the corresponding values
    t_sim(end+1,:) = t;
    lip_state(end+1,:) = x_lip;
    lip_control(end+1,:) = u_lip;
    dp_state(end+1,:) = x;

end

% Animate the double pendulum's motion
addpath("../balancer")
animate_lip_dp(1:dt:T,dp_state, lip_state, lip_control, h)

% Calculate the actual CoM position for the recorded trajectory
com_state = x_com(dp_state')';

% Plot the errors along the trajectories
subplot(2,2,1);
hold on
plot(t_sim, lip_state(:,1)-com_state(:,1));
title("CoM X position");
xlabel("time")
ylabel("error")

subplot(2,2,2);
hold on
plot(t_sim, lip_state(:,3)-com_state(:,3));
title("CoM Y position");
xlabel("time")
ylabel("error")

subplot(2,2,3);
hold on
plot(t_sim, lip_state(:,2)-com_state(:,2));
title("CoM X velocity");
xlabel("time")
ylabel("error")

subplot(2,2,4);
hold on
plot(t_sim, lip_state(:,4)-com_state(:,4));
title("CoM Y velocity");
xlabel("time")
ylabel("error")

% Plot the net error and Bisimulation function
err = zeros(size(t_sim));
V_sim = zeros(size(t_sim));
for t = 1:length(t_sim)
    err(t) = (lip_state(t,:)-com_state(t,:))*(lip_state(t,:)-com_state(t,:))';
    V_sim(t) = SimulationFcn(lip_state(t,:)',dp_state(t,:)');
end

figure;
hold on
plot(t_sim, err);
plot(t_sim, V_sim);
legend("Output Error","Simulation Function");
xlabel("time")

