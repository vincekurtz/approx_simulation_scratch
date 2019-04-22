function [] = animate_lip(t_sim, x_sim, u_sim, h)
    % Animate the behavior of the LIP model.
    %
    % Parameters:
    %    t_sim : a 1xT*dt vector of timesteps
    %    x_sim : a 2xT*dt sequence of states representing CoM position and velocity
    %    u_sim : a 1xT*dt vector of control inputs representing the position of the CoP
    %    h     : the (fixed) height of the CoM

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
    ylim([-1,3])

    for i=1:length(t_sim)-1
        tic
        dt = t_sim(i+1) - t_sim(i);

        % Move the center of pressure
        u = u_sim(i);
        frame1.Matrix = makehgtform('translate',[u 0 0]);

        % move the center of mass
        x = x_sim(i,1);
        frame2.Matrix = makehgtform('translate',[x-u h 0]);

        l.XData = [0,x-u];
        l.YData = [0,h];

        pause(dt-toc)
    end
end
