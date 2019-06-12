function tau = QPTracker(u_lip, x_lip, q, qd, params)
    % Given the state and control of the reduced-order LIP model, 
    % formulate and solve a Quadratic Program that attempts to track
    % the reduced-order model while also enforcing contact constraints.

    % Actual CoM position and velocity
    pp_com = p_com(q);  % avoid name conflict
    pd_com = J(q)*qd;

    % LIP model CoM position and velocity
    p_com_des = x_lip(1:2);
    pd_com_des = x_lip(4:5);

    a_lip = params.omega^2*(p_com_des(1)-u_lip);

    % Desired CoM acceleration
    Kp = 30*eye(2);
    Kd = 20*eye(2);
    a_com_des = [a_lip;0] - Kp*(pp_com - p_com_des) - Kd*(pd_com - pd_com_des);

    % Temporary check: direct feedback linearization
    %tau = J(q)'*(Lambda(q)*a_com_des - Lambda(q)*Jdot(q,qd)*qd + Lambda(q)*J(q)*inv(H(q))*C(q,qd));

    % Quadratic Program to enforce contact constraints
    opti = casadi.Opti();
    qdd = opti.variable(4,1);
    tau = opti.variable(4,1);

    % Cost Function
    cost = J(q)*qdd + Jdot(q,qd)*qd - a_com_des;
    cost = cost'*cost;
    cost = cost + qdd'*diag([0;0;0.3;0.3])*qdd;  % penalty to discourage spinning top links
    opti.minimize(cost);

    % Constraints
    opti.subject_to( H(q)*qdd+C(q,qd) == tau );
    opti.subject_to( -50 <= tau <= 50 );
    hd_com = A_com(q)*qdd + Ad_com_qd(q,qd);
    opti.subject_to( A_cwc(pp_com)*hd_com <= b_cwc(pp_com) );

    options.ipopt.print_level = 0;
    options.print_time = false;
    opti.solver('ipopt', options);
    sol = opti.solve();
    tau = sol.value(tau);


