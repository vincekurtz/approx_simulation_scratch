%% Use the LIP model to balance a double pendulum which is simulated in Gazebo via ROS

clear;
clc;

% Launch the gazebo simulation (note that matlab overwrites the
% LD_LIBRARY_PATH, so we need to set it manually first
[status,cmdout] = system(["./start_gazebo_sim.sh"]);

pause(5)  % wait a few seconds so that the gazebo simulation can start roscore

% Start a matlab ros node
ROS_MASTER_URI="http://ubuntu-p50:11311/";
rosinit(ROS_MASTER_URI)

try
    % Connect to the physics un-pausing service
    pause_client = rossvcclient('/gazebo/unpause_physics');
    pause_msg = rosmessage(pause_client);  % an empty message
    
    % create a subscriber to the joint_states topic
    %global theta1 theta2 theta1_dot theta2_dot;
    %joint_sub = rossubscriber('/double_pendulum/joint_states', @joint_state_callback);
    joint_sub = rossubscriber('/double_pendulum/joint_states');
    
    pause(2) % Wait 2s
    
    % unpause the physics engine
    pause_resp = call(pause_client, pause_msg);
    
    % Read the current joint states
    data = receive(joint_sub, 10);  % timeout of 10s
    data.showdetails

    pause

    % Shut down the gazebo simulation
    cleanupFcn(cmdout);
catch E
    disp(E)
    cleanupFcn(cmdout);
end

function joint_state_callback(data)
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