function u_lip = LIPController(x_lip, omega, dt)
    % Compute a control for the (simplified) LIP model
   
    import casadi.*
    
    % Horizon length 
    N = 10;
    
    % LIP dynamics
    A_lip = [0 1; omega^2 0]; 
    B_lip = [0 ; -omega^2];

    % Cost function definition
    Q = diag([1;1]);
    R = 10;
    [K,Qf] = lqr(A_lip,B_lip,Q,R); % Terminal cost based on optimal cost-to-go of LQR

    % Single Shooting
    opti = casadi.Opti();
    u = opti.variable(1,N);

    J = 0;
    x = x_lip(1:2);
    for t=1:N
        J = J + x'*Q*x + u(:,t)'*R*u(:,t);

        % Forward euler integration
        x = x + (A_lip*x + B_lip*u(:,t))*dt;
    end
    J = J + x'*Qf*x;
    opti.minimize(J);

    options.ipopt.print_level = 0;
    options.print_time = false;
    opti.solver('ipopt', options);
    sol = opti.solve();
    u_lip = sol.value(u(1));

end

