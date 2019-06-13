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
    QQ = blkdiag(QQ,Qf);           % the last element is the terminal cost
    RR = kron(eye(params.N-1),R);  % matrix with R on diagonal
    J = x_lip(:)'*QQ*x_lip(:) + u_lip(:)'*RR*u_lip(:);

    % Constraints
    opti.subject_to(x_lip(:,1) == x_lip_init);
    opti.subject_to(x_com(:,1) == x_com_init);

    % Constraints on the contact wrench (note that we only consider motion in the x-y plane)
    mu = 0.2;
    l = 0.5;
    mg = -39.24;
    A = [0 0  0  0  -1 0;   % positive normal force
         0 0  0  1 -mu 0;   % Coulomb friction
         0 0  0 -1 -mu 0;
         0 0  1  0  -l 0;   % Center of pressure constraint
         0 0 -1  0  -l 0];

    for t=1:params.N-1
        % Interface constraint
        opti.subject_to( u_com(:,t) == params.R*u_lip(:,t) + params.Q*x_lip(:,t) + params.K*(x_com(:,t)-x_lip(:,t)) );

        % Contact constraint
      
        % Attempt to use acceleration bound to linearize contact constraints
        bnd = 3.0;
        opti.subject_to( -bnd <= u_com(2,t) <= bnd );  % bound linear acceleration
        opti.subject_to( -bnd <= u_com(3,t) <= bnd );

        Acwc{t,1} = A*[0;0;u_com(:,t);0] + A*[S([0;mg;0])-S([bnd;bnd;0]);zeros(3)]*[x_com(1:2,t);0];
        Acwc{t,2} = A*[0;0;u_com(:,t);0] + A*[S([0;mg;0])-S([bnd;-bnd;0]);zeros(3)]*[x_com(1:2,t);0];
        Acwc{t,3} = A*[0;0;u_com(:,t);0] + A*[S([0;mg;0])-S([-bnd;bnd;0]);zeros(3)]*[x_com(1:2,t);0];
        Acwc{t,4} = A*[0;0;u_com(:,t);0] + A*[S([0;mg;0])-S([-bnd;-bnd;0]);zeros(3)]*[x_com(1:2,t);0];

        opti.subject_to( Acwc{t,1} <= A*[0;0;0;0;mg;0] );
        opti.subject_to( Acwc{t,2} <= A*[0;0;0;0;mg;0] );
        opti.subject_to( Acwc{t,3} <= A*[0;0;0;0;mg;0] );
        opti.subject_to( Acwc{t,4} <= A*[0;0;0;0;mg;0] );

        %% Original contact constraint
        %Acwc{t} = A*[0;0;u_com(:,t);0] + A*[S([0;mg;0])-S([u_com(2:3);0]);zeros(3)]*[x_com(1:2,t);0];
        %opti.subject_to( Acwc{t} <= A*[0;0;0;0;mg;0] );

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

function y = S(x)
    % Skew symmetric cross product matrix
    y = [0     -x(3)   x(2);
         x(3)   0     -x(1);
         -x(2)  x(1)    0  ];
end
