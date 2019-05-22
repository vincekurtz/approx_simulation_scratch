%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use Gazebo/Ros to imulate balencing the multilink balancer via
% planning with the LIP model.
%
% Assumes that create_multilink_spatial_model has already been
% run, generating a saved model and functions representing the
% balancer dynamics and relevant quantities (H,C,J,Jdot,etc).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

% Launch the gazebo simulation
[status,PID] = system(["./start_gazebo_sim.sh"])

pause(5)  % wait a few seconds so that the gazebo simulation can start roscore

try
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Controller Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Derive LIP and linearized CoM dynamics along with an interface
    setup_interface;

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
    
    % Define a signal to indicate that we have recieved joint angle data
    global ready
    ready = false;
    
    % create a subscriber to read joint angles and velocities
    global t1 t2 t3 t4 t1_dot t2_dot t3_dot t4_dot;
    joint_sub = rossubscriber('/multilink_balancer/joint_states', @joint_state_callback);

    % Create publishers for the joint torques
    joint1_pub = rospublisher('/multilink_balancer/joint1_torque_controller/command', 'std_msgs/Float64');
    joint2_pub = rospublisher('/multilink_balancer/joint2_torque_controller/command', 'std_msgs/Float64');
    joint3_pub = rospublisher('/multilink_balancer/joint3_torque_controller/command', 'std_msgs/Float64');
    joint4_pub = rospublisher('/multilink_balancer/joint4_torque_controller/command', 'std_msgs/Float64');
    joint1_msg = rosmessage(joint1_pub);  % empty messages of the correct type
    joint2_msg = rosmessage(joint2_pub);
    joint3_msg = rosmessage(joint3_pub);
    joint4_msg = rosmessage(joint4_pub);
    
    % Number of timesteps and time discritization
    T = 5;  % simulation time in seconds
    dt = 3e-2;  % note that we get joint angles from ROS at ~50Hz
  
    pause(2) % Wait 2s to initialize ROS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % unpause the physics engine
    pause_resp = call(pause_client, pause_msg);

    % Wait to get joint angle readings
    while ~ready
        disp(ready)
        pause(0.001)
    end

    % Get the initial state of the robot: note that there is a mismatch between the
    % definitions of theta1, theta2 in our matlab model and in gazebo
    q0 = [pi/2-t1; -t2; -t3; -t4];
    qd0 = [-t1_dot; -t2_dot; -t3_dot; -t4_dot];
    x = [q0;qd0];

    % Set the initial state of the LIP model
    com_pos = p_com(q0);
    com_vel = pd_com(q0,qd0);
    x_lip = [com_pos(1);com_vel(1);h;0];

    % Initial value of the simulation function: we can treat as an error bound
    % between the LIP model state and the true model state
    %V = SimulationFcn(x_lip, x)

    % Record trajectories
    joint_trajectory = [];
    com_trajectory = [];
    lip_trajectory = [];
    lip_control = [];

    for timestep=1:dt:T
        tic

        %tic
        % Compute the LIP control that will let us balance.
        u_lip = -K_lip*x_lip(1:2);
   
        % Compute double pendulum state (joint angle definitions differ from Gazebo)
        q = [pi/2-t1; -t2; -t3; -t4];
        qd = [-t1_dot; -t2_dot; -t3_dot; -t4_dot];
        x = [q;qd];

        % Compute the state of the linearized system (CoM x and y position and velocity)
        com_pos = p_com(q);
        com_vel = pd_com(q,qd);
        x_com = [com_pos(1);com_vel(1);com_pos(2);com_vel(2)];

        % Compute virtual control for the feedback linearized system
        u_com = R*u_lip + Q*x_lip + K_joint*(x_com-P*x_lip);

        % Compute torques to apply to the full system
        tau = J(q)'*(Lambda(q)*u_com - Lambda(q)*Jdot(q,qd)*qd + Lambda(q)*J(q)*inv(H(q))*C(q,qd));

        % Apply the torques to the full system  
        joint1_msg.Data = min(tau_max, max(tau_min, -tau(1)));   % torque limits + correct for angle definitions
        joint2_msg.Data = min(tau_max, max(tau_min, -tau(2)));
        joint3_msg.Data = min(tau_max, max(tau_min, -tau(3)));
        joint4_msg.Data = min(tau_max, max(tau_min, -tau(4)));
        send(joint1_pub, joint1_msg)
        send(joint2_pub, joint2_msg)
        send(joint3_pub, joint3_msg)
        send(joint4_pub, joint4_msg)
        
        % Apply the LIP control to the LIP model
        dx_lip = A2*x_lip + B2*u_lip;
        x_lip = x_lip + dx_lip*dt;

        % Record the resulting trajectories
        joint_trajectory(end+1,:) = x;
        com_trajectory(end+1,:) = x_com;
        lip_trajectory(end+1,:) = x_lip;
        lip_control(end+1,:) = u_lip;
        %toc
        
        pause(dt-toc)
    end

    disp("Simulation Finished")

    % Shut down the gazebo simulation
    cleanupFcn(PID);
catch E
    disp(getReport(E))
    cleanupFcn(PID);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Playback and analysis of the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('balancer_model')  % import spatial_v2 model
q_trajectory = joint_trajectory(:,1:model.NB)';
showmotion(model, 1:dt:T, q_trajectory);

plot(com_trajectory)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function joint_state_callback(~, data)
    % a callback function that sets global variables that represent the
    % angles and velocities of the goints
    global t1 t2 t3 t4 t1_dot t2_dot t3_dot t4_dot;

    % Indicate that we are ready to start the simulation
    global ready;
    ready = true;
    
    t1 = data.Position(1);
    t2 = data.Position(2);
    t3 = data.Position(3);
    t4 = data.Position(4);
    t1_dot = data.Velocity(1);
    t2_dot = data.Velocity(2);
    t3_dot = data.Velocity(3);
    t4_dot = data.Velocity(4);
end

function cleanupFcn(pid)
    % A helper function to do all the things we need to do on exit
    disp("Running Cleanup Function")
    
    % Close the gazebo simulation
    system(['kill ', pid]);
    
    % shutdown the ros node
    rosshutdown;
end