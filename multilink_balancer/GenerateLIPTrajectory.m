function u_lip = GenerateLIPTrajectory(x_lip_init, x_com_init, params)
    % Compute a trajectory for the LIP model that respects contact constraints
    % of the full model.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameterizing the Cost Function %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Running Cost
    Q = diag([100;0;0;10;0]);
    R = 0.01;

    % Terminal cost based on optimal cost-to-go of LQR
    Qf = [32.1091 0 0  0.1544 0;
           0      0 0  0      0;
           0      0 0  0      0;
           0.1544 0 0  0.0491 0;
           0      0 0  0      0];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization via Direct Collocation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    opti = casadi.Opti();
    
    u_lip = opti.variable(1,params.N-1);
    x_lip = opti.variable(5,params.N);
    u_com = opti.variable(3,params.N-1);
    x_com = opti.variable(5,params.N);

    % Cost Function
    QQ = kron(eye(params.N-1),Q);  % matrix with Q on diagonal
    QQ = blkdiag(QQ,2*Qf);           % the last element is the terminal cost
    RR = kron(eye(params.N-1),R);  % matrix with R on diagonal
    J = x_lip(:)'*QQ*x_lip(:) + u_lip(:)'*RR*u_lip(:);

    % Constraints
    opti.subject_to(x_lip(:,1) == x_lip_init);
    opti.subject_to(x_com(:,1) == x_com_init);

    for t=1:params.N-1
        % Interface constraint
        opti.subject_to( u_com(:,t) == params.R*u_lip(:,t) + params.Q*x_lip(:,t) + params.K*(x_com(:,t)-x_lip(:,t)) );

        % Contact constraint
        opti.subject_to( A_cwc(x_com(:,t))*u_com(:,t) <= b_cwc(x_com(:,t)) );

        % Euler integration constraints
        opti.subject_to(x_lip(:,t+1) == x_lip(:,t) + (params.A_lip*x_lip(:,t) + params.B_lip*u_lip(:,t))*params.dt);
        opti.subject_to(x_com(:,t+1) == x_com(:,t) + (params.A_com*x_com(:,t) + params.B_com*u_com(:,t))*params.dt);

    end

    % Run the optimization
    opti.minimize(J);

    options.ipopt.print_level = 0;
    options.print_time = false;
    opti.solver('ipopt', options);
    opti.set_initial(u_lip,1);
    sol = opti.solve();
    u_lip = sol.value(u_lip);

end

