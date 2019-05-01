%% Use the LIP model to balance a double pendulum which is simulated in Gazebo via ROS

clear;
clc;

% Launch the gazebo simulation (note that matlab overwrites the
% LD_LIBRARY_PATH, so we need to set it manually first
[status,cmdout] = system(["./start_gazebo_sim.sh"])

pause(5)  % wait a few seconds so that the gazebo simulation can start roscore

try
    %% Simulation Setup
    % Start a matlab ros node
    ROS_MASTER_URI="http://ubuntu-p50:11311/";
    rosinit(ROS_MASTER_URI)

    % Connect to the physics un-pausing service
    pause_client = rossvcclient('/gazebo/unpause_physics');
    pause_msg = rosmessage(pause_client);  % an empty message
    
    % create a subscriber to read joint angles and velocities
    global theta1 theta2 theta1_dot theta2_dot;
    joint_sub = rossubscriber('/double_pendulum/joint_states', @joint_state_callback);
    
    % Create publishers for the joint torques
    joint1_pub = rospublisher('/double_pendulum/joint1_torque_controller/command', 'std_msgs/Float64');
    joint2_pub = rospublisher('/double_pendulum/joint2_torque_controller/command', 'std_msgs/Float64');
    joint1_msg = rosmessage(joint1_pub);  % empty messages of the correct type
    joint2_msg = rosmessage(joint2_pub);
    
    % Number of timesteps and time discritization
    sim_time = 5;  % simulation time in seconds
    dt = 1e-3;
    Ns = sim_time/dt;
    
    pause(2) % Wait 2s
    
    %% Controller Setup
    kp_1 = 50;  % proportional controller for each torque
    kp_2 = 30;
    
    %% Simulation
    % unpause the physics engine
    pause_resp = call(pause_client, pause_msg);
    
    for i=1:Ns        
        % Compute torques to apply
        tau1 = -kp_1*theta1;
        tau2 = -kp_2*theta2;
        
        % Apply the torques
        joint1_msg.Data = tau1;
        joint2_msg.Data = tau2;
        send(joint1_pub, joint1_msg)
        send(joint2_pub, joint2_msg)
        
        pause(dt)
    end

    % Shut down the gazebo simulation
    cleanupFcn(cmdout);
catch E
    disp(E)
    cleanupFcn(cmdout);
end

function joint_state_callback(~, data)
    % a callback function that sets global variables that represent the
    % angles of both joints
    global theta1 theta2 theta1_dot theta2_dot;
    
    theta1 = data.Position(1);
    theta2 = data.Position(2);
    theta1_dot = data.Velocity(1);
    theta2_dot = data.Velocity(2);
end

function cleanupFcn(pid)
    % A helper function to do all the things we need to do on exit
    disp("Running Cleanup Function")
    
    % Close the gazebo simulation
    system(['kill ', pid]);
    
    % shutdown the ros node
    rosshutdown;
end