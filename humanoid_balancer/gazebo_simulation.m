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

    % Indicate the approach to tracking the LIP model we take
    %   'QP' : Traditional approach using a Quadratic Program
    %   'AS' : Our approach using Approximate Simulation
    tracking_method = 'AS';

    % Derive LIP and linearized CoM dynamics along with an interface
    setup_interface;

    % Torque limits
    tau_min = -100;
    tau_max = 100;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Start a matlab ros node
    ROS_MASTER_URI="http://ubuntu-p50:11311/";
    rosinit(ROS_MASTER_URI)

    % Connect to the physics un-pausing service
    pause_client = rossvcclient('/gazebo/unpause_physics');
    pause_msg = rosmessage(pause_client);  % an empty message

    % Connect to a service that will allow us to apply a "push" force to the model
    push_client = rossvcclient('/gazebo/apply_body_wrench');
    push_msg = rosmessage(push_client);
    
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
    dt = 5e-2;  % note that we get joint angles from ROS at ~50Hz
 
    % Compute MPC constraints for the LIP model
    params.A_lip = A2;
    params.B_lip = B2;
    params.A_com = A1;
    params.B_com = B1;
    params.N = 5;     % MPC horizon
    params.dt = dt;
    params.R = R;
    params.Q = Q;
    params.K = K_joint;
    params.m = m;
    constraint_params = LIPConstraints(params);
   
    pause(2) % Wait 2s to initialize ROS

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % unpause the physics engine
    pause_resp = call(pause_client, pause_msg);

    % Wait to get joint angle readings
    while ~ready
        pause(0.001)
    end

    % Move the balancer to an initial position
    for timestep=1:dt:T
        tic
        % Get joint states (joint angle definitions differ from Gazebo)
        q = [pi/2-t1; -t2; -t3; -t4];
        qd = [-t1_dot; -t2_dot; -t3_dot; -t4_dot];
        x_com = [p_com(q);A_com(q)*qd];

        % Move to desired initial condition
        x_com_des = [0.0 h 0 0 0]';
        u_com = -K_joint*(x_com_des-x_com);
        u_com = constrain_ucom(u_com, x_com);

        % Compute and apply torques
        Lambda = inv(A_com(q)*inv(H(q))*A_com(q)');
        tau = A_com(q)'*(Lambda*u_com - Lambda*Ad_com_qd(q,qd) + Lambda*A_com(q)*inv(H(q))*C(q,qd));
        
        % Secondary control objective via null space projector
        Abar = inv(H(q))*A_com(q)'*Lambda;
        N = (eye(4) - A_com(q)'*Abar')';   % null space projector
        kp = diag([0;0;1;1]);    % we'll use PD control of joints to
        kd = diag([1;1;1;1]);    % apply commands in the null space
        q_des = [0;pi/4;-pi/2;-pi/4];
        qd_des = [0;0;0;0];
        tau0 = -kp*(q-q_des)-kd*(qd-qd_des);
        tau = tau + N'*tau0;

        joint1_msg.Data = min(tau_max, max(tau_min, -tau(1)));   % torque limits + correct for angle definitions
        joint2_msg.Data = min(tau_max, max(tau_min, -tau(2)));
        joint3_msg.Data = min(tau_max, max(tau_min, -tau(3)));
        joint4_msg.Data = min(tau_max, max(tau_min, -tau(4)));
        send(joint1_pub, joint1_msg)
        send(joint2_pub, joint2_msg)
        send(joint3_pub, joint3_msg)
        send(joint4_pub, joint4_msg)

        pause(dt-toc)
    end

    % Give the robot a "push" by applying a force to the torso link
    disp("Applying Push")
    push_msg.BodyName = 'multilink_balancer::link4';
    push_msg.ReferencePoint.Z = 1.0;
    push_msg.Wrench.Force.X = -150;
    push_msg.Duration.Nsec = 1e7;
    push_resp = call(push_client, push_msg);
    pause(0.01);

    % Get the initial state of the robot: note that there is a mismatch between the
    % definitions of theta1, theta2 in our matlab model and in gazebo
    q0 = [pi/2-t1; -t2; -t3; -t4];
    qd0 = [-t1_dot; -t2_dot; -t3_dot; -t4_dot];
    x = [q0;qd0];

    % Set the initial state of the LIP model and the linearized CoM system
    com_pos = p_com(q0);
    com_h = A_com(q0)*qd0;
    x_com = [com_pos;com_h];
    x_lip = [com_pos(1);h;0;com_h(2);0];  % y position fixed at h, angular momentum fixed at 0.
    
    % "practice" generating a lip trajectory: this will help warm-start the MPC solver
    %GenerateLIPTrajectory(x_lip, x_com, constraint_params);

    % Record trajectories
    joint_trajectory = [];
    com_trajectory = [];
    lip_trajectory = [];
    lip_control = [];
    torques = [];

    for timestep=dt:dt:T
        disp("")
        tic

        % Compute double pendulum state (joint angle definitions differ from Gazebo)
        q = [pi/2-t1; -t2; -t3; -t4];
        qd = [-t1_dot; -t2_dot; -t3_dot; -t4_dot];
        x = [q;qd];

        % Compute the state of the linearized system (CoM x and y position and velocity)
        com_pos = p_com(q);
        com_h = A_com(q)*qd;
        x_com = [com_pos;com_h];

        if tracking_method == 'AS'
            % Use our approach, which certifies approximate simulation

            % Constrain the LIP model in MPC according to the CWC criterion
            u_lip_trajectory = GenerateLIPTrajectory(x_lip, x_com, constraint_params);
            u_lip = u_lip_trajectory(1);

            % Apply the LIP control to the LIP model
            dx_lip = A2*x_lip + B2*u_lip;
            x_lip = x_lip + dx_lip*dt;

            % Compute the associated virtual control
            u_com = R*u_lip + Q*x_lip + K_joint*(x_com-P*x_lip);

            % Compute torques to apply to the full system
            Lambda = inv(A_com(q)*inv(H(q))*A_com(q)');
            tau = A_com(q)'*(Lambda*u_com - Lambda*Ad_com_qd(q,qd) + Lambda*A_com(q)*inv(H(q))*C(q,qd));
        
            % Secondary control objective via null space projector
            Abar = inv(H(q))*A_com(q)'*Lambda;
            N = (eye(4) - A_com(q)'*Abar')';   % null space projector

            kp = diag([0;0;1;0.1]);    % we'll use PD control of joints to
            kd = diag([2;2;1;1]);    % apply commands in the null space
            q_des = [0;pi/4;-pi/2;-pi/4];
            qd_des = [0;0;0;0];
            tau0 = -kp*(q-q_des)-kd*(qd-qd_des);
            tau = tau + N'*tau0;

        elseif tracking_method == 'QP'
            % Use a traditional Quadratic Program to generate torques that track
            % the LIP model

            % Use LQR to find a trajectory for the LIP model
            u_lip = -K_lip*[x_lip(1);1/m*x_lip(4)];

            % Apply the LIP control to the LIP model
            dx_lip = A2*x_lip + B2*u_lip;
            x_lip = x_lip + dx_lip*dt;

            % Track this trajectory with a QP that enforces contact constraints
            params.omega = omega;
            tau = QPTracker(u_lip, x_lip, q, qd, params);
            
        end

        % Apply the torques to the full system  
        joint1_msg.Data = min(tau_max, max(tau_min, -tau(1)));   % torque limits + correct for angle definitions
        joint2_msg.Data = min(tau_max, max(tau_min, -tau(2)));
        joint3_msg.Data = min(tau_max, max(tau_min, -tau(3)));
        joint4_msg.Data = min(tau_max, max(tau_min, -tau(4)));
        send(joint1_pub, joint1_msg)
        send(joint2_pub, joint2_msg)
        send(joint3_pub, joint3_msg)
        send(joint4_pub, joint4_msg)

        % Record the resulting trajectories
        joint_trajectory(end+1,:) = x;
        com_trajectory(end+1,:) = x_com;
        lip_trajectory(end+1,:) = x_lip;
        lip_control(end+1,:) = u_lip;
        torques(end+1,:) = tau;
        toc
        
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

% Plot the error vs the simulation function
err_sim = [];
V_sim = [];
for i=1:length(com_trajectory)
    x_com = com_trajectory(i,:)';
    x_lip = lip_trajectory(i,:)';

    err_sim(end+1) = sqrt((x_com-x_lip)'*(x_com-x_lip));
    V_sim(end+1) = sqrt((x_com-P*x_lip)'*M*(x_com-P*x_lip));
end
hold on
plot(dt:dt:T,err_sim)
plot(dt:dt:T,V_sim)

legend("Error","Simulation Function")
xlabel("time")
ylabel("value")

% Plot the CoM trajectory
figure;
plot(dt:dt:T,com_trajectory)
legend("x position", "y position", "angular momentum", "x velocity", "y velocity");

% Animate the LIP trajectory
addpath("../balancer")
figure;
lip_traj = [lip_trajectory(:,1);1/m*lip_trajectory(:,4)];
animate_lip(dt:dt:T, lip_traj, lip_control,h)


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
    t4 = data.Position(4); t1_dot = data.Velocity(1);
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
