function [] = animate_double_pendulum(t_sim, x_sim)
    % Animate the motion of our simple double pendulum model
    %
    % Parameters:
    %   t_sim : a 1xT*dt vector of timestamps
    %   x_sim : a 4xT*dt vector of states x = [theta1;theta2;theta1_dot;theta2_dot]

    arm_width = 0.1;
    arm_length = 1;  % our model assumes arms of unit length
    
    % Ground frame
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
    xlim([-2.5,2.5])
    ylim([-1,3])

    for t = 1:length(t_sim)-1
        tic
        dt = t_sim(t+1) - t_sim(t);

        theta1 = x_sim(t,1);
        theta2 = x_sim(t,2);

        frame1.Matrix = makehgtform('zrotate',theta1);
        frame2.Matrix = makehgtform('translate',[1 0 0],'zrotate',theta2);

        pause(dt-toc)
    end

end

