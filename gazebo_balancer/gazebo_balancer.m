%% Use the LIP model to balance a double pendulum which is simulated in Gazebo via ROS

clear;
clc;

% Tools for animating the results

% Launch the gazebo simulation (note that matlab overwrites the
% LD_LIBRARY_PATH, so we need to set it manually first
[status,cmdout] = system(["./start_gazebo_sim.sh"])

pause(5)  % wait a few seconds so that the gazebo simulation can start roscore

try
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Controller Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute an interface that guarantees approximate simulation.
    compute_interface;

    % Torque limits
    tau_min = -1000;
    tau_max = 1000;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Start a matlab ros node
    ROS_MASTER_URI="http://ubuntu-p50:11311/";
    rosinit(ROS_MASTER_URI)

    % Connect to the physics un-pausing service
    pause_client = rossvcclient('/gazebo/unpause_physics');
    pause_msg = rosmessage(pause_client);  % an empty message
    
    % create a subscriber to read joint angles and velocities
    global t1 t2 t1_dot t2_dot;
    joint_sub = rossubscriber('/double_pendulum/joint_states', @joint_state_callback);
    
    % Create publishers for the joint torques
    joint1_pub = rospublisher('/double_pendulum/joint1_torque_controller/command', 'std_msgs/Float64');
    joint2_pub = rospublisher('/double_pendulum/joint2_torque_controller/command', 'std_msgs/Float64');
    joint1_msg = rosmessage(joint1_pub);  % empty messages of the correct type
    joint2_msg = rosmessage(joint2_pub);
    
    % Number of timesteps and time discritization
    T = 5;  % simulation time in seconds
    dt = 3e-2;
    
    pause(2) % Wait 2s
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % unpause the physics engine
    pause_resp = call(pause_client, pause_msg);

    pause(0.1)

    % Get the initial state of the robot: note that there is a mismatch between the
    % definitions of theta1, theta2 in our matlab model and in gazebo
    x = [pi/2-t1;-t2;-t1_dot;-t2_dot];

    % Set the initial state of the LIP model
    com_init = x_com(x);
    x_lip = [com_init(1:2);h;0];
    x_lip = [0.6517;0.3131;h;0];

    % Record trajectories
    joint_trajectory = [];
    com_trajectory = [];
    lip_trajectory = [];
    lip_control = [];
        
    for timestep=1:dt:T
        tic
        % Compute the LIP control that will let us balance.
        u_lip = LIPController(x_lip, omega, dt);
        
        % Apply the LIP control to the LIP model
        dx_lip = A2*x_lip + B2*u_lip;
        x_lip = x_lip + dx_lip*dt;

        % Get full system state
        x = [pi/2-t1;-t2;-t1_dot;-t2_dot];

        % Compute torques to apply to the full system
        tau = InterfaceFcn(u_lip, x_lip, x);

        % Correct for different angle definitions
        tau1 = -tau(1);
        tau2 = -tau(2);
        
        % Apply the torques to the full system
        joint1_msg.Data = min(tau_max, max(tau_min, tau1));
        joint2_msg.Data = min(tau_max, max(tau_min, tau2));

        send(joint1_pub, joint1_msg)
        send(joint2_pub, joint2_msg)

        % Record the resulting trajectories
        joint_trajectory(end+1,:) = x;
        com_trajectory(end+1,:) = x_com(x);
        lip_trajectory(end+1,:) = x_lip;
        lip_control(end+1,:) = u_lip;
        
        pause(dt-toc)
    end

    disp("Simulation Finished")

    % Shut down the gazebo simulation
    cleanupFcn(cmdout);
catch E
    disp(getReport(E))
    cleanupFcn(cmdout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Playback and analysis of the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the net error and Bisimulation function
t_sim = 1:dt:T;
err = zeros(size(t_sim));
V_sim = zeros(size(t_sim));
for t = 1:length(t_sim)
    err(t) = (lip_trajectory(t,:)-com_trajectory(t,:))*(lip_trajectory(t,:)-com_trajectory(t,:))';
    V_sim(t) = SimulationFcn(lip_trajectory(t,:)',joint_trajectory(t,:)');
end
figure;
hold on
plot(t_sim, err);
plot(t_sim, V_sim);
legend("Output Error","Simulation Function");
xlabel("time")

% Play back an animation of both the lip and the full model
addpath("../balancer")
animate_lip_dp(1:dt:T, joint_trajectory, lip_trajectory, lip_control, h)


function joint_state_callback(~, data)
    % a callback function that sets global variables that represent the
    % angles of both joints
    global t1 t2 t1_dot t2_dot;
    
    t1 = data.Position(1);
    t2 = data.Position(2);
    t1_dot = data.Velocity(1);
    t2_dot = data.Velocity(2);
end

function cleanupFcn(pid)
    % A helper function to do all the things we need to do on exit
    disp("Running Cleanup Function")
    
    % Close the gazebo simulation
    system(['kill ', pid]);
    
    % shutdown the ros node
    rosshutdown;
end
