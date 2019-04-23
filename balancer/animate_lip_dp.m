function [] = animate_lip_dp(t_sim, x_dp, x_lip, u_lip, h)
    % Animate the motion of our simple double pendulum model and the reduced-order
    % LIP model
    %
    % Parameters:
    %    t_sim : a 1xT*dt vector of timestamps
    %    x_dp  : a 4xT*dt vector of dp states x = [theta1;theta2;theta1_dot;theta2_dot]
    %    x_lip : a 2xT*dt sequence of states representing CoM position and velocity of the LIP model
    %    u_sim : a 1xT*dt vector of control inputs representing the position of the CoP of the LIP
    %    h     : the (fixed) height of the CoM of the LIP

    figure('Position',[100 100 1000 500]);

    x_lim = [-2.5,2.5];
    y_lim = [-1,3];

    % Double Pendulum animation setup
    subplot(1,2,1)

    arm_width = 0.1;
    arm_length = 1;  % our model assumes arms of unit length
   
    ground = plot(gca, [-15 15],[0 0],'k','LineWidth',2);
    frame0 = hgtransform(gca);

    % First arm
    frame1 = hgtransform(frame0);
    arm1 = rectangle('Parent', frame1);
    arm1.Position = [-arm_width/2 -arm_width/2 arm_length+arm_width, arm_width];
    arm1.FaceColor = 'red';

    % Second arm
    frame2 = hgtransform(frame1);
    frame2.Matrix = makehgtform('translate',[1 0 0]);
    arm2 = rectangle('Parent', frame2);
    arm2.Position = [-arm_width/2 -arm_width/2 arm_length+arm_width, arm_width];
    arm2.FaceColor = 'red';
    
    axis equal
    xlim(x_lim)
    ylim(y_lim)
    title("Rigid Body Model")

    % LIP animation setup
    subplot(1,2,2)

    ground_lip = plot(gca, [-15 15],[0 0],'k','LineWidth',2);
    frame0_lip = hgtransform(gca);

    % Center of pressure
    frame1_lip = hgtransform(frame0_lip);
    r = 0.1;
    cop = rectangle('Curvature',[1,1],'Parent',frame1_lip);  % Center of pressure
    cop.Position = [-r/2, -r/2, r, r];
    cop.EdgeColor = 'none';
    cop.FaceColor = 'blue';

    % Center of mass
    frame2_lip = hgtransform(frame1_lip);
    cop = rectangle('Curvature',[1,1],'Parent',frame2_lip);  % Center of pressure
    cop.Position = [-r/2, -r/2, r, r];
    cop.EdgeColor = 'none';
    cop.FaceColor = 'red';

    % Line between COP and COM
    l = line('Parent',frame1_lip);

    % Line that COM is constrained to
    yline(h,'--');

    axis equal
    xlim(x_lim)
    ylim(y_lim)
    title("Reduced Order Model")

    for t = 1:length(t_sim)-1
        tic
        dt = t_sim(t+1) - t_sim(t);

        % Update the Rigid Body model
        theta1 = x_dp(t,1);
        theta2 = x_dp(t,2);
        frame1.Matrix = makehgtform('zrotate',theta1);
        frame2.Matrix = makehgtform('translate',[1 0 0],'zrotate',theta2);

        % Update the LIP model
        u = u_lip(t);
        frame1_lip.Matrix = makehgtform('translate',[u 0 0]);  % Move the CoP
        x = x_lip(t,1);
        frame2_lip.Matrix = makehgtform('translate',[x-u h 0]);    % move the CoM

        l.XData = [0,x-u];   % update the line between them
        l.YData = [0,h];

        pause(dt-toc)      % Keep as close to real time as possible
    end

end

