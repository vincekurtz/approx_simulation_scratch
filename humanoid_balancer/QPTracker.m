function tau = QPTracker(u_lip, x_lip, q, qd, omega)
    % Given the state and control of the reduced-order LIP model, 
    % formulate and solve a Quadratic Program that attempts to track
    % the reduced-order model while also enforcing contact constraints.

    % Actual CoM position and velocity
    pp_com = p_com(q);  % avoid name conflict
    pd_com = J(q)*qd;

    % LIP model CoM position and velocity
    p_com_des = x_lip(1:2);
    pd_com_des = x_lip(4:5);

    a_lip = omega^2*(p_com_des(1)-u_lip);

    % Desired CoM acceleration
    Kp = 30*eye(2);
    Kd = 20*eye(2);
    a_com_des = [a_lip;0] - Kp*(pp_com - p_com_des) - Kd*(pd_com - pd_com_des);

    % Quadratic Program to enforce contact constraints
    opti = casadi.Opti();
    qdd = opti.variable(4,1);
    tau = opti.variable(4,1);

    % Main Cost Function
    cost = J(q)*qdd + Jdot(q,qd)*qd - a_com_des;
    cost = cost'*cost;

    % Secondary cost to try to stay in desired configuration
    q_des = [pi/4;pi/2;-pi/4;-pi/2];   
    qd_des = [0;0;0;0];
    qdd_des = -1*(q - q_des) - 1*(qd-qd_des);
    cost = cost + (qdd-qdd_des)'*0.1*eye(4)*(qdd-qdd_des);

    opti.minimize(cost);

    % Constraints
    opti.subject_to( H(q)*qdd+C(q,qd) == tau );
    opti.subject_to( -100 <= tau <= 100 );
    hd_com = A_com(q)*qdd + Ad_com_qd(q,qd);
    opti.subject_to( A_cwc(pp_com)*hd_com <= b_cwc(pp_com) );

    options.ipopt.print_level = 0;
    options.print_time = false;
    opti.solver('ipopt', options);
    sol = opti.solve();
    tau = sol.value(tau);


