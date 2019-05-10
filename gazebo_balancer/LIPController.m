function u_lip = LIPController(x_lip, params)
    % Compute a control for the (simplified) LIP model
   
    tic
    import casadi.*
    
    % Horizon length 
    N = 5;
   
    % LIP dynamics
    A_lip = params.A_lip;
    B_lip = params.B_lip;

    % Running Cost
    Q = diag([100;10;0;0]);
    R = 0.01;

    % Terminal cost based on optimal cost-to-go of LQR
    %[K,Qf] = lqr(A_lip,B_lip,Q,R)
    Qf = [31.8676 0.0772 0 0;
           0.0772 0.0244 0 0;
           0      0      0 0;
           0      0      0 0];

    % Direct Collocation
    opti = casadi.Opti();
    u = opti.variable(1,N);
    x = opti.variable(4,N);

    J = 0;
    opti.subject_to( x(:,1) == x_lip );  % initial constraint
    for t=1:N-1
        J = J + x(:,t)'*Q*x(:,t) + u(:,t)'*R*u(:,t);

        % Forward euler integration constraint
        opti.subject_to( x(:,t+1) == x(:,t) + (A_lip*x(:,t) + B_lip*u(:,t))*params.dt );

    end
    J = J + x(:,N)'*Qf*x(:,N);
    opti.minimize(J);

    % Control constraints
    opti.subject_to(u <= params.u_max);

    options.ipopt.print_level = 0;
    options.print_time = false;
    opti.solver('ipopt', options);
    sol = opti.solve();
    u_lip = sol.value(u(1));

    % Check the validity of this control
    check_ulip(u_lip, params)

    toc

end

